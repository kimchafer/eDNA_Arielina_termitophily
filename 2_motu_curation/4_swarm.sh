#!/usr/bin/env bash
#SBATCH --job-name=4_swarm
#SBATCH --output=4_swarm_%j.out
#SBATCH --error=4_swarm_%j.err
#SBATCH --partition=Node4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=30-00:00:00

# ==============================================================================
# PURPOSE
#   Cluster subsequence-resolved COI ASVs with Swarm and reconstruct MOTU tables.
#   Run one user-selected d per job; default d=2.
#
# INPUT
#   /storage2/jihoonkim/Arielina_eDNA/analyses/1_20261103/3_subsequence/asvs.subsequence_resolved.fasta
#   /storage2/jihoonkim/Arielina_eDNA/analyses/1_20261103/3_subsequence/asv_id_to_sequence_resolved.tsv
#   /storage2/jihoonkim/Arielina_eDNA/analyses/1_20261103/3_subsequence/sample_by_asv_resolved.tsv
#
# OUTPUT (example: d=2)
#   4_swarm_motu/ALL_d2/ALL_d2.swarm             Swarm clusters
#   4_swarm_motu/ALL_d2/ALL_d2.stats             Cluster statistics
#   4_swarm_motu/ALL_d2/ALL_d2_seeds_raw.fasta   Raw seed FASTA
#   4_swarm_motu/ALL_d2/ALL_d2_members.tsv       ASV-to-MOTU membership
#   4_swarm_motu/ALL_d2/ALL_d2_MOTU.fasta        MOTU representatives
#   4_swarm_motu/ALL_d2/ALL_d2_MOTU_table.tsv    Sample abundances
#   4_swarm_motu/ALL_d2/ALL_d2_cluster_summary.tsv
#   4_swarm_motu/ALL_d2/ALL_d2_sample_summary.tsv
#   4_swarm_motu/ALL_d2/ALL_d2_summary.tsv
#   4_swarm_motu/ALL_d2/logs/                    Parameters and checksum
#
# SETTINGS
#   d=2 by default; override with first sbatch argument.
#   d=1 uses --fastidious (-f); d>=2 uses ordinary clustering.
#   --usearch-abundance (-z); default Swarm scoring and OTU-breaking.
#   No read removal or taxonomic filtering; preserve sample read totals.
#   Existing output directories are protected against overwriting.
# ==============================================================================
set -euo pipefail

WORK_DIR="/storage2/jihoonkim/Arielina_eDNA/analyses/1_20261103"
D_VAL="${1:-2}"
THREADS="${SLURM_CPUS_PER_TASK:-4}"
SCRIPT="${WORK_DIR}/4_swarm.sh"
INPUT_DIR="${WORK_DIR}/3_subsequence"
INPUT_FASTA="${INPUT_DIR}/asvs.subsequence_resolved.fasta"
INPUT_MAP="${INPUT_DIR}/asv_id_to_sequence_resolved.tsv"
INPUT_COUNTS="${INPUT_DIR}/sample_by_asv_resolved.tsv"
POOL="ALL_d${D_VAL}"
OUT_DIR="${WORK_DIR}/4_swarm_motu/${POOL}"
LOG_DIR="${OUT_DIR}/logs"
PREFIX="${OUT_DIR}/${POOL}"
SWARM_FASTA="${PREFIX}_swarm_input.fasta"
CLUSTERS="${PREFIX}.swarm"
STATS="${PREFIX}.stats"
SEEDS="${PREFIX}_seeds_raw.fasta"

# ==============================================================================
# PREFLIGHT
# ==============================================================================

[[ "$D_VAL" =~ ^[1-9][0-9]*$ ]] || {
    echo "[ERROR] d must be a positive integer: $D_VAL" >&2
    exit 1
}
[[ -d "$WORK_DIR" && -d "$INPUT_DIR" ]] || {
    echo "[ERROR] Missing work or subsequence directory" >&2
    exit 1
}
for file in "$INPUT_FASTA" "$INPUT_MAP" "$INPUT_COUNTS"; do
    [[ -s "$file" ]] || {
        echo "[ERROR] Missing or empty input: $file" >&2
        exit 1
    }
done
[[ ! -e "$OUT_DIR" ]] || {
    echo "[ERROR] Existing output protected: $OUT_DIR" >&2
    exit 1
}

source /home/JihoonKim/miniconda3/etc/profile.d/conda.sh
conda activate dada2_env
for exe in swarm python3; do
    command -v "$exe" >/dev/null 2>&1 || {
        echo "[ERROR] Missing executable: $exe" >&2
        exit 1
    }
done

mkdir -p "$LOG_DIR"
{
    printf 'parameter\tvalue\n'
    printf 'date\t%s\n' "$(date -Is)"
    printf 'work_dir\t%s\n' "$WORK_DIR"
    printf 'input_fasta\t%s\n' "$INPUT_FASTA"
    printf 'input_map\t%s\n' "$INPUT_MAP"
    printf 'input_counts\t%s\n' "$INPUT_COUNTS"
    printf 'swarm_version\t%s\n' "$(swarm -v 2>&1 | head -n 1)"
    printf 'python_version\t%s\n' "$(python3 --version 2>&1)"
    printf 'distance_d\t%s\n' "$D_VAL"
    printf 'fastidious\t%s\n' "$([[ "$D_VAL" -eq 1 ]] && printf true || printf false)"
    printf 'usearch_abundance\ttrue\n'
    printf 'default_scoring\ttrue\n'
    printf 'default_otu_breaking\ttrue\n'
    printf 'threads\t%s\n' "$THREADS"
} > "${LOG_DIR}/parameters.tsv"

echo "[INFO] Start: $(date)"
echo "[INFO] Working directory: $WORK_DIR"
echo "[INFO] Swarm: $(swarm -v 2>&1 | head -n 1)"
echo "[INFO] Distance d: $D_VAL"
echo "[INFO] Threads: $THREADS"
echo "[INFO] Output directory: $OUT_DIR"

# ==============================================================================
# STEP 1. VALIDATE ASVS AND PREPARE SWARM INPUT FASTA
# ==============================================================================

python3 - "$INPUT_FASTA" "$INPUT_MAP" "$INPUT_COUNTS" "$SWARM_FASTA" <<'PY'
import csv
import sys

fasta_path, map_path, counts_path, output_path = sys.argv[1:]


def records(path):
    with open(path, newline='') as fh:
        r = csv.DictReader(fh, delimiter='\t')
        if not r.fieldnames or len(r.fieldnames) != len(set(r.fieldnames)):
            raise ValueError(f'Invalid TSV header: {path}')
        return list(r), r.fieldnames


def fasta(path):
    seqs, name, parts = {}, None, []
    with open(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith('>'):
                if name is not None:
                    if name in seqs:
                        raise ValueError(f'Duplicate FASTA ID: {name}')
                    seqs[name] = ''.join(parts).upper()
                name, parts = line[1:].split()[0], []
            else:
                if name is None:
                    raise ValueError('FASTA lacks initial header')
                parts.append(line)
    if name is not None:
        if name in seqs:
            raise ValueError(f'Duplicate FASTA ID: {name}')
        seqs[name] = ''.join(parts).upper()
    return seqs


def count(raw, label):
    try:
        value = int(raw)
    except (ValueError, TypeError) as exc:
        raise ValueError(f'Invalid count for {label}: {raw}') from exc
    if value < 0:
        raise ValueError(f'Negative abundance for {label}')
    return value


seqs = fasta(fasta_path)
meta, meta_cols = records(map_path)
rows, cols = records(counts_path)
if not seqs or not meta or not rows:
    raise ValueError('Empty ASV input')
if not {'ASV_ID', 'sequence', 'total_reads'} <= set(meta_cols):
    raise ValueError('Missing ASV mapping columns')
ids = [r['ASV_ID'] for r in meta]
if len(ids) != len(set(ids)) or set(ids) != set(seqs):
    raise ValueError('ASV ID sets differ across FASTA and mapping')
if set(cols) != set(ids) | {'Sample'} or len(cols) != len(ids) + 1:
    raise ValueError('Abundance matrix must have Sample and all FASTA ASVs')
samples = [row['Sample'] for row in rows]
if any(not s.strip() for s in samples) or len(samples) != len(set(samples)):
    raise ValueError('Invalid or duplicate sample ID')

input_reads = 0
for row in meta:
    aid = row['ASV_ID']
    seq = row['sequence'].upper()
    total = count(row['total_reads'], aid)
    if not seq or seq != seqs[aid] or set(seq) - set('ACGT'):
        raise ValueError(f'Invalid or mismatched sequence: {aid}')
    if 'length_bp' in row and count(row['length_bp'], aid) != len(seq):
        raise ValueError(f'Sequence length mismatch: {aid}')
    if total < 1 or total != sum(
        count(r[aid], r['Sample'] + '/' + aid) for r in rows
    ):
        raise ValueError(f'Abundance mismatch: {aid}')
    input_reads += total
if len(set(seqs.values())) != len(seqs):
    raise ValueError('Duplicate nucleotide sequences in Swarm input')

with open(output_path, 'w') as fh:
    for row in meta:
        aid = row['ASV_ID']
        fh.write(f'>{aid};size={int(row["total_reads"])};\n{seqs[aid]}\n')
print(f'[QC PASS] Input: {len(ids)} ASVs, {input_reads:,} reads, {len(samples)} samples')
PY

# ==============================================================================
# STEP 2. RUN SWARM CLUSTERING
# ==============================================================================

echo "[STEP 2] Swarm clustering d=$D_VAL"
if (( D_VAL == 1 )); then
    swarm -d "$D_VAL" -f -t "$THREADS" -z \
        -o "$CLUSTERS" -s "$STATS" -w "$SEEDS" "$SWARM_FASTA"
else
    swarm -d "$D_VAL" -t "$THREADS" -z \
        -o "$CLUSTERS" -s "$STATS" -w "$SEEDS" "$SWARM_FASTA"
fi
for file in "$CLUSTERS" "$STATS" "$SEEDS"; do
    [[ -s "$file" ]] || {
        echo "[ERROR] Missing Swarm output: $file" >&2
        exit 1
    }
done
echo "[QC PASS] Swarm clusters, statistics and seeds created"

# ==============================================================================
# STEP 3. RECONSTRUCT MOTU TABLES AND VALIDATE CLUSTERS
# ==============================================================================

python3 - "$INPUT_FASTA" "$INPUT_MAP" "$INPUT_COUNTS" \
    "$CLUSTERS" "$STATS" "$SEEDS" "$PREFIX" "$D_VAL" <<'PY'
import csv
import re
import sys

(fa_path, map_path, counts_path, clusters_path,
 stats_path, seeds_path, prefix, d_str) = sys.argv[1:]
d = int(d_str)


def read_tsv(path):
    with open(path, newline='') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def write(suffix, cols, rows):
    with open(prefix + suffix, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def read_fasta(path):
    seqs, head, parts = {}, None, []
    with open(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith('>'):
                if head is not None:
                    if head in seqs:
                        raise ValueError(f'Duplicate representative header: {head}')
                    seqs[head] = ''.join(parts).upper()
                head, parts = line[1:].split()[0], []
            else:
                parts.append(line)
    if head is not None:
        if head in seqs:
            raise ValueError(f'Duplicate representative header: {head}')
        seqs[head] = ''.join(parts).upper()
    return seqs


def parse_label(value):
    match = re.fullmatch(r'([^;\s]+)(?:;size=([0-9]+);?)?', value)
    if match is None:
        raise ValueError(f'Unexpected Swarm label: {value!r}')
    return match.group(1), int(match.group(2)) if match.group(2) is not None else None


meta = read_tsv(map_path)
sample_rows = read_tsv(counts_path)
samples = [row['Sample'] for row in sample_rows]
seqs = read_fasta(fa_path)
all_ids = [row['ASV_ID'] for row in meta]
totals = {row['ASV_ID']: int(row['total_reads']) for row in meta}
if set(seqs) != set(totals) or len(all_ids) != len(totals):
    raise ValueError('ASV input consistency failed')
with open(clusters_path) as fh:
    clusters = [line.split() for line in fh if line.strip()]
with open(stats_path) as fh:
    stats = [line.rstrip('\r\n').split('\t') for line in fh if line.strip()]
seeds = read_fasta(seeds_path)
if not clusters or not (len(clusters) == len(stats) == len(seeds)):
    raise ValueError('Swarm cluster, statistics, and seed counts disagree')

seen, seen_seeds = set(), set()
members, cluster_qc, table, long_table, final_seqs = [], [], [], [], []
for number, (tokens, stat) in enumerate(zip(clusters, stats), 1):
    if len(stat) < 4:
        raise ValueError(f'Malformed Swarm statistics row: {number}')
    ids = [parse_label(token)[0] for token in tokens]
    if (not ids or len(set(ids)) != len(ids) or
            any(aid not in totals or aid in seen for aid in ids)):
        raise ValueError(f'Duplicate or unknown ASV in MOTU {number}')
    for token in tokens:
        aid, size = parse_label(token)
        if size is not None and size != totals[aid]:
            raise ValueError(f'Swarm membership abundance differs for {aid}')
    seen.update(ids)
    seed, _ = parse_label(stat[2])
    if seed not in ids or seed in seen_seeds:
        raise ValueError(f'Invalid or reused seed in MOTU {number}: {seed}')
    seen_seeds.add(seed)
    total = sum(totals[aid] for aid in ids)
    if (int(stat[0]), int(stat[1]), int(stat[3])) != (len(ids), total, totals[seed]):
        raise ValueError(f'Swarm statistics mismatch for MOTU {number}')
    motu = f'd{d}_MOTU{number:06d}'
    cluster_qc.append(dict(MOTU_ID=motu, seed_ASV=seed, n_ASVs=len(ids),
                           total_reads=total, members=','.join(ids)))
    for aid in ids:
        members.append(dict(ASV_ID=aid, MOTU_ID=motu, seed_ASV=seed,
                            ASV_reads=totals[aid], MOTU_reads=total))
    per_sample = {
        sample: sum(int(row[aid]) for aid in ids)
        for sample, row in zip(samples, sample_rows)
    }
    if sum(per_sample.values()) != total:
        raise ValueError(f'Sample counts disagree for {motu}')
    table.append(dict(MOTU_ID=motu, **per_sample))
    for sample in samples:
        if per_sample[sample]:
            long_table.append(dict(MOTU_ID=motu, sample=sample, reads=per_sample[sample]))
    final_seqs.append((motu, seqs[seed]))

# Validate every ASV, seed identity, and cluster-level abundance.
if seen != set(totals):
    raise ValueError(f'ASVs omitted from Swarm: {sorted(set(totals)-seen)}')
seed_labels = {}
for header, sequence in seeds.items():
    aid, size = parse_label(header)
    if aid in seed_labels:
        raise ValueError(f'Duplicate representative ID: {aid}')
    seed_labels[aid] = (size, sequence)
if set(seed_labels) != seen_seeds:
    raise ValueError('Representative FASTA seeds disagree with Swarm statistics')
for c in cluster_qc:
    size, sequence = seed_labels[c['seed_ASV']]
    if sequence != seqs[c['seed_ASV']] or (size is not None and size != c['total_reads']):
        raise ValueError(f'Representative sequence or abundance mismatch: {c["MOTU_ID"]}')

# Export original memberships, representative sequences, and MOTU abundances.
write('_members.tsv', ['ASV_ID', 'MOTU_ID', 'seed_ASV', 'ASV_reads', 'MOTU_reads'], members)
write('_cluster_summary.tsv', ['MOTU_ID', 'seed_ASV', 'n_ASVs', 'total_reads', 'members'], cluster_qc)
write('_MOTU_table.tsv', ['MOTU_ID'] + samples, table)
write('_MOTU_abundance_long.tsv', ['MOTU_ID', 'sample', 'reads'], long_table)
with open(prefix + '_MOTU.fasta', 'w') as fh:
    for motu, sequence in final_seqs:
        fh.write(f'>{motu}\n{sequence}\n')

# Validate read conservation within each sample and globally.
sample_qc = []
for sample, row in zip(samples, sample_rows):
    incoming = sum(int(row[aid]) for aid in all_ids)
    outgoing = sum(motu[sample] for motu in table)
    if incoming != outgoing:
        raise ValueError(f'Read conservation failed for sample {sample}')
    sample_qc.append(dict(
        sample=sample, input_reads=incoming, output_reads=outgoing,
        input_detected_ASVs=sum(int(row[aid]) > 0 for aid in all_ids),
        output_detected_MOTUs=sum(motu[sample] > 0 for motu in table)
    ))
write('_sample_summary.tsv', ['sample', 'input_reads', 'output_reads',
                              'input_detected_ASVs', 'output_detected_MOTUs'], sample_qc)
input_reads = sum(totals.values())
if sum(c['total_reads'] for c in cluster_qc) != input_reads:
    raise ValueError('Global read conservation failed')
summary = [
    ('d', d), ('input_ASVs', len(all_ids)), ('input_reads', input_reads),
    ('MOTUs', len(clusters)),
    ('singleton_ASV_MOTUs', sum(c['n_ASVs'] == 1 for c in cluster_qc)),
    ('multi_ASV_MOTUs', sum(c['n_ASVs'] > 1 for c in cluster_qc)),
    ('ASVs_in_multi_ASV_MOTUs', sum(c['n_ASVs'] for c in cluster_qc if c['n_ASVs'] > 1)),
    ('maximum_ASVs_in_MOTU', max(c['n_ASVs'] for c in cluster_qc)),
    ('output_reads', sum(c['total_reads'] for c in cluster_qc))
]
write('_summary.tsv', ['metric', 'value'], [dict(metric=k, value=v) for k, v in summary])

print('=' * 66)
print(f'SWARM CLUSTERING SUMMARY (d={d})')
print('=' * 66)
for metric, value in summary:
    print(f'{metric:34s} {value:>12,d}')
print('\nSAMPLE-LEVEL READ CONSERVATION')
for row in sample_qc:
    print(f'{row["sample"]}: {row["input_reads"]:,} -> {row["output_reads"]:,}; '
          f'{row["input_detected_ASVs"]} ASVs -> {row["output_detected_MOTUs"]} MOTUs')
print('\nLARGEST MULTI-ASV MOTUS')
for row in sorted((c for c in cluster_qc if c['n_ASVs'] > 1),
                  key=lambda c: (-c['n_ASVs'], -c['total_reads']))[:15]:
    print(f'{row["MOTU_ID"]}: {row["n_ASVs"]} ASVs, '
          f'{row["total_reads"]:,} reads: {row["members"]}')
print('\n[QC PASS] All ASVs assigned once; every sample and all reads conserved')
PY

# ==============================================================================
# STEP 4. VALIDATE FINAL OUTPUTS AND RECORD SCRIPT CHECKSUM
# ==============================================================================

for name in \
    "${POOL}.swarm" \
    "${POOL}.stats" \
    "${POOL}_seeds_raw.fasta" \
    "${POOL}_members.tsv" \
    "${POOL}_cluster_summary.tsv" \
    "${POOL}_MOTU_table.tsv" \
    "${POOL}_MOTU_abundance_long.tsv" \
    "${POOL}_MOTU.fasta" \
    "${POOL}_sample_summary.tsv" \
    "${POOL}_summary.tsv"; do
    [[ -s "${OUT_DIR}/${name}" ]] || {
        echo "[ERROR] Missing or empty output: ${OUT_DIR}/${name}" >&2
        exit 1
    }
done
if [[ -f "$SCRIPT" ]]; then
    sha256sum "$SCRIPT" > "${LOG_DIR}/script_sha256.txt"
fi

echo "[QC PASS] Swarm output files validated"
echo "[DONE] Swarm clustering completed: $(date)"
echo "[OUTPUT] ${PREFIX}_MOTU.fasta"
echo "[OUTPUT] ${PREFIX}_MOTU_table.tsv"
echo "[OUTPUT] ${PREFIX}_summary.tsv"

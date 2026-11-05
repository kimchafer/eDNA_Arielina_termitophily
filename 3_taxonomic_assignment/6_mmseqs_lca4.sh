#!/usr/bin/env bash
#SBATCH --job-name=6_mmseqs_lca4
#SBATCH --output=6_mmseqs_lca4_%j.out
#SBATCH --error=6_mmseqs_lca4_%j.err
#SBATCH --partition=Node7
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=400G
#SBATCH --time=30-00:00:00

# ==============================================================================
# PURPOSE
#   Assign taxonomy to LULU-curated MOTU representatives using MMseqs2.
#   Compare LCA mode 4 at minimum sequence identities of 80% and 97%.
#   Only --min-seq-id differs between the two conditions.
#
# INPUT
#   /storage2/jihoonkim/Arielina_eDNA/analyses/1_20261103/5_lulu/ALL_d2/LULU_curated_MOTU.fasta
#   /storage2/jihoonkim/Arielina_eDNA/analyses/1_20261103/5_lulu/ALL_d2/LULU_curated_MOTU_table.tsv
#   /storage/woojunbang/mmseq2_for_taxonomy/nt.fnaDB (taxonomy database)
#
# OUTPUT
#   6_mmseqs_lca4/input/lulu_motu_representatives.fasta  Validated MOTU FASTA
#   6_mmseqs_lca4/input/motu_ids.tsv                 Expected MOTU IDs
#   6_mmseqs_lca4/db/lulu_motu                       MMseqs2 query database
#   6_mmseqs_lca4/results/taxonomy_lca4_a20_80.tsv  Taxonomy at 80% identity
#   6_mmseqs_lca4/results/taxonomy_lca4_a20_97.tsv  Taxonomy at 97% identity
#   6_mmseqs_lca4/logs/input_summary.tsv             Sample-level input QC
#   6_mmseqs_lca4/logs/input_metrics.tsv             Input MOTU/read counts
#   6_mmseqs_lca4/logs/parameters.tsv                Analysis parameters
#   6_mmseqs_lca4/logs/taxonomy_qc.tsv               Output assignment QC
#   6_mmseqs_lca4/logs/script_sha256.txt             Script checksum
#
# SETTINGS
#   Swarm d=2, LULU-curated representatives; no hard-coded MOTU/read counts.
#   LCA mode 4, min-aln-len 300, max-accept 20.
#   Search type 3, alignment mode 3, strand 2, ORF filter off.
#   Sensitivity 5, e-value 1e-5, query coverage 0.80, cov-mode 2.
#   Minimum sequence identities: 0.80 and 0.97 (independent runs).
#   No additional taxonomic filtering or modification of LULU abundances.
#   Existing output directory is protected against overwriting.
# ==============================================================================
set -euo pipefail

WORK_DIR="/storage2/jihoonkim/Arielina_eDNA/analyses/1_20261103"
REF_DB="/storage/woojunbang/mmseq2_for_taxonomy/nt.fnaDB"
THREADS="${SLURM_CPUS_PER_TASK:-40}"
SCRIPT="${WORK_DIR}/6_mmseqs_lca4.sh"

LULU_DIR="${WORK_DIR}/5_lulu/ALL_d2"
INPUT_FASTA="${LULU_DIR}/LULU_curated_MOTU.fasta"
INPUT_TABLE="${LULU_DIR}/LULU_curated_MOTU_table.tsv"

OUT_DIR="${WORK_DIR}/6_mmseqs_lca4"
INPUT_DIR="${OUT_DIR}/input"
QUERY_FASTA="${INPUT_DIR}/lulu_motu_representatives.fasta"
QUERY_IDS="${INPUT_DIR}/motu_ids.tsv"
QUERY_DB="${OUT_DIR}/db/lulu_motu"
RESULTS_DIR="${OUT_DIR}/results"
TMP_DIR="${OUT_DIR}/tmp"
LOG_DIR="${OUT_DIR}/logs"
INPUT_SUMMARY="${LOG_DIR}/input_summary.tsv"
INPUT_METRICS="${LOG_DIR}/input_metrics.tsv"
TAXONOMY_QC="${LOG_DIR}/taxonomy_qc.tsv"

# ==============================================================================
# PREFLIGHT
# ==============================================================================

source /home/JihoonKim/miniconda3/etc/profile.d/conda.sh
conda activate mm

for exe in mmseqs python3 sha256sum; do
    command -v "$exe" >/dev/null 2>&1 || {
        echo "[ERROR] Missing executable: $exe" >&2
        exit 1
    }
done

[[ -d "$WORK_DIR" && -d "$LULU_DIR" ]] || {
    echo "[ERROR] Work or LULU directory does not exist" >&2
    exit 1
}
for file in "$INPUT_FASTA" "$INPUT_TABLE" "${REF_DB}.dbtype"; do
    [[ -s "$file" ]] || {
        echo "[ERROR] Missing or empty input: $file" >&2
        exit 1
    }
done
[[ ! -e "$OUT_DIR" ]] || {
    echo "[ERROR] Existing output protected: $OUT_DIR" >&2
    exit 1
}

mkdir -p "$INPUT_DIR" "${OUT_DIR}/db" "$RESULTS_DIR" "$TMP_DIR" "$LOG_DIR"

echo "[INFO] Start: $(date)"
echo "[INFO] MMseqs2: $(mmseqs version)"
echo "[INFO] Threads: $THREADS"
echo "[INFO] Input: $INPUT_FASTA"
echo "[INFO] Reference DB: $REF_DB"

# ==============================================================================
# STEP 1. VALIDATE LULU MOTUS AND PREPARE QUERY FASTA
# ==============================================================================

echo "[STEP 1] Validate LULU-curated MOTU sequences and abundances"

python3 - "$INPUT_FASTA" "$INPUT_TABLE" "$QUERY_FASTA" "$QUERY_IDS" \
    "$INPUT_SUMMARY" "$INPUT_METRICS" <<'PY'
import csv
import sys
from pathlib import Path

fasta_path, table_path, query_path, ids_path, summary_path, metrics_path = map(Path, sys.argv[1:])

sequences = {}
order = []
identifier = None
parts = []

with fasta_path.open() as handle:
    for raw in handle:
        line = raw.strip()
        if not line:
            continue
        if line.startswith('>'):
            if identifier is not None:
                if identifier in sequences:
                    raise RuntimeError(f'Duplicate FASTA ID: {identifier}')
                sequences[identifier] = ''.join(parts).upper()
                order.append(identifier)
            header = line[1:].split()
            if not header:
                raise RuntimeError('Empty FASTA header')
            identifier, parts = header[0], []
        else:
            if identifier is None:
                raise RuntimeError('FASTA sequence without a header')
            parts.append(line)
    if identifier is not None:
        if identifier in sequences:
            raise RuntimeError(f'Duplicate FASTA ID: {identifier}')
        sequences[identifier] = ''.join(parts).upper()
        order.append(identifier)

if not sequences or any(not s or set(s) - set('ACGT') for s in sequences.values()):
    raise RuntimeError('Empty or invalid nucleotide sequences in curated FASTA')

with table_path.open(newline='') as handle:
    reader = csv.DictReader(handle, delimiter='\t')
    fields = reader.fieldnames
    if not fields or fields[0] != 'MOTU_ID' or len(fields) < 2 or len(fields) != len(set(fields)):
        raise RuntimeError('Curated table requires unique MOTU_ID and sample columns')
    samples = fields[1:]
    rows = list(reader)

if not rows:
    raise RuntimeError('Curated MOTU table is empty')
ids = [row['MOTU_ID'] for row in rows]
if len(ids) != len(set(ids)) or set(ids) != set(sequences):
    raise RuntimeError('Curated FASTA and table MOTU identifiers disagree')
if any(not sample.strip() for sample in samples):
    raise RuntimeError('Blank sample identifier')

sample_totals = {sample: 0 for sample in samples}
sample_detected = {sample: 0 for sample in samples}
for row in rows:
    mid = row['MOTU_ID']
    motu_total = 0
    for sample in samples:
        value = row[sample]
        try:
            n = int(value)
        except (ValueError, TypeError) as exc:
            raise RuntimeError(f'Invalid abundance for {mid}/{sample}: {value!r}') from exc
        if n < 0:
            raise RuntimeError(f'Negative abundance for {mid}/{sample}')
        sample_totals[sample] += n
        sample_detected[sample] += int(n > 0)
        motu_total += n
    if motu_total == 0:
        raise RuntimeError(f'Zero-abundance curated MOTU: {mid}')

with query_path.open('w') as handle:
    for mid in order:
        handle.write(f'>{mid}\n{sequences[mid]}\n')
with ids_path.open('w', newline='') as handle:
    writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
    writer.writerow(['MOTU_ID'])
    writer.writerows((mid,) for mid in order)
with summary_path.open('w', newline='') as handle:
    writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
    writer.writerow(['sample', 'input_reads', 'detected_MOTUs'])
    writer.writerows((s, sample_totals[s], sample_detected[s]) for s in samples)
with metrics_path.open('w', newline='') as handle:
    writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
    writer.writerow(['metric', 'value'])
    writer.writerows([
        ('query_MOTUs', len(sequences)),
        ('input_samples', len(samples)),
        ('input_reads', sum(sample_totals.values())),
    ])

print(f'[QC PASS] Curated input: {len(sequences)} MOTUs, {sum(sample_totals.values()):,} reads')
for sample in samples:
    print(f'[QC] {sample}: {sample_totals[sample]:,} reads; {sample_detected[sample]} MOTUs')
PY

QUERY_COUNT="$(awk -F '\t' '$1 == "query_MOTUs" {print $2}' "$INPUT_METRICS")"
INPUT_READS="$(awk -F '\t' '$1 == "input_reads" {print $2}' "$INPUT_METRICS")"
[[ "$QUERY_COUNT" =~ ^[1-9][0-9]*$ && "$INPUT_READS" =~ ^[1-9][0-9]*$ ]] || {
    echo "[ERROR] Invalid dynamic input counts" >&2
    exit 1
}

{
    printf 'parameter\tvalue\n'
    printf 'date\t%s\n' "$(date -Is)"
    printf 'work_dir\t%s\n' "$WORK_DIR"
    printf 'input_fasta\t%s\n' "$INPUT_FASTA"
    printf 'input_table\t%s\n' "$INPUT_TABLE"
    printf 'input_fasta_sha256\t%s\n' "$(sha256sum "$INPUT_FASTA" | cut -d' ' -f1)"
    printf 'reference_db\t%s\n' "$REF_DB"
    printf 'mmseqs_version\t%s\n' "$(mmseqs version)"
    printf 'representative_fasta_sha256\t%s\n' "$(sha256sum "$QUERY_FASTA" | cut -d' ' -f1)"
    printf 'query_MOTUs\t%s\n' "$QUERY_COUNT"
    printf 'input_reads\t%s\n' "$INPUT_READS"
    printf 'lca_mode\t4\n'
    printf 'min_seq_id\t0.80,0.97\n'
    printf 'min_aln_len\t300\n'
    printf 'search_type\t3\n'
    printf 'orf_filter\t0\n'
    printf 'alignment_mode\t3\n'
    printf 'strand\t2\n'
    printf 'sensitivity\t5\n'
    printf 'evalue\t1e-5\n'
    printf 'coverage\t0.80\n'
    printf 'cov_mode\t2\n'
    printf 'max_accept\t20\n'
    printf 'max_rejected\t50\n'
    printf 'tax_lineage\t1\n'
    printf 'threads\t%s\n' "$THREADS"
} > "${LOG_DIR}/parameters.tsv"

# ==============================================================================
# STEP 2. CREATE MMSEQS2 QUERY DATABASE
# ==============================================================================

echo "[STEP 2] Create query database ($QUERY_COUNT MOTUs)"
mmseqs createdb "$QUERY_FASTA" "$QUERY_DB"
[[ -s "${QUERY_DB}.dbtype" ]] || {
    echo "[ERROR] Missing MMseqs2 query database" >&2
    exit 1
}

echo "[QC PASS] Query database created"

# ==============================================================================
# STEP 3. RUN TAXONOMY WITH LCA MODE 4 AT 80% AND 97% IDENTITY
# ==============================================================================

run_taxonomy() {
    local label="$1"
    local identity="$2"
    local result="${RESULTS_DIR}/taxonomy_lca4_a20_${label}"
    local out_tsv="${result}.tsv"
    local qc_row="${LOG_DIR}/taxonomy_lca4_a20_${label}_qc.tsv"

    echo "[STEP 3] LCA mode 4 / minimum identity ${identity}: $(date)"

    mmseqs taxonomy \
        "$QUERY_DB" \
        "$REF_DB" \
        "$result" \
        "${TMP_DIR}/lca4_a20_${label}" \
        --tax-lineage 1 \
        --threads "$THREADS" \
        --search-type 3 \
        --orf-filter 0 \
        --lca-mode 4 \
        --alignment-mode 3 \
        --strand 2 \
        -s 5 \
        -e 1e-5 \
        --min-seq-id "$identity" \
        --min-aln-len 300 \
        -c 0.80 \
        --cov-mode 2 \
        --max-accept 20 \
        --max-rejected 50

    mmseqs createtsv "$QUERY_DB" "$result" "$out_tsv" --threads "$THREADS"

    [[ -s "$out_tsv" ]] || {
        echo "[ERROR] Missing or empty taxonomy TSV: $out_tsv" >&2
        exit 1
    }

    python3 - "$QUERY_IDS" "$out_tsv" "$qc_row" "$label" "$identity" <<'PY'
import csv
import sys
from pathlib import Path

ids_path, result_path, qc_path = map(Path, sys.argv[1:4])
label, identity = sys.argv[4:6]
with ids_path.open(newline='') as handle:
    reader = csv.DictReader(handle, delimiter='\t')
    expected = [row['MOTU_ID'] for row in reader]

observed = []
with result_path.open(newline='') as handle:
    for row_number, record in enumerate(csv.reader(handle, delimiter='\t'), 1):
        if len(record) != 5:
            raise RuntimeError(f'Unexpected taxonomy columns ({len(record)}) in row {row_number}')
        observed.append(record[0])

if len(observed) != len(expected) or len(observed) != len(set(observed)) or set(observed) != set(expected):
    missing = sorted(set(expected) - set(observed))
    unexpected = sorted(set(observed) - set(expected))
    raise RuntimeError(f'Taxonomy ID/count mismatch: expected={len(expected)}, '
                       f'observed={len(observed)}, missing={missing[:10]}, '
                       f'unexpected={unexpected[:10]}')

with qc_path.open('w', newline='') as handle:
    writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
    writer.writerow(['condition', 'min_seq_id', 'expected_MOTUs', 'taxonomy_rows', 'unique_MOTU_IDs', 'id_check'])
    writer.writerow([label, identity, len(expected), len(observed), len(set(observed)), 'PASS'])
print(f'[QC PASS] {label}: {len(observed)} unique taxonomy rows; all MOTUs represented')
PY

    echo "[OUTPUT] $out_tsv"
}

run_taxonomy 80 0.80
run_taxonomy 97 0.97

# ==============================================================================
# STEP 4. SUMMARIZE QC AND RECORD SCRIPT CHECKSUM
# ==============================================================================

{
    head -n 1 "${LOG_DIR}/taxonomy_lca4_a20_80_qc.tsv"
    tail -n +2 "${LOG_DIR}/taxonomy_lca4_a20_80_qc.tsv"
    tail -n +2 "${LOG_DIR}/taxonomy_lca4_a20_97_qc.tsv"
} > "$TAXONOMY_QC"

if [[ -f "$SCRIPT" ]]; then
    sha256sum "$SCRIPT" > "${LOG_DIR}/script_sha256.txt"
fi

echo "[QC PASS] Both taxonomy assignments validated for ${QUERY_COUNT} MOTUs"
echo "[DONE] MMseqs2 LCA mode 4 comparisons completed: $(date)"
echo "[OUTPUT] ${RESULTS_DIR}/taxonomy_lca4_a20_80.tsv"
echo "[OUTPUT] ${RESULTS_DIR}/taxonomy_lca4_a20_97.tsv"
echo "[OUTPUT] $TAXONOMY_QC"

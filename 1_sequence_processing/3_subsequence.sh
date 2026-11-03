#!/usr/bin/env bash
#SBATCH --job-name=3_subsequence
#SBATCH --output=3_subsequence_%j.out
#SBATCH --error=3_subsequence_%j.err
#SBATCH --partition=Node4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=30-00:00:00

# ==============================================================================
# PURPOSE
#   Identify exact nucleotide containment among DADA2 COI ASVs and resolve
#   unambiguous terminally truncated sequences without losing sample-level reads.
#
# INPUT
#   /storage2/jihoonkim/Arielina_eDNA/analyses/1_20261103/2_dada2/ASVs.nochim.fasta
#   /storage2/jihoonkim/Arielina_eDNA/analyses/1_20261103/2_dada2/ASV_ID_to_sequence.tsv
#   /storage2/jihoonkim/Arielina_eDNA/analyses/1_20261103/2_dada2/sample_by_ASV.tsv
#
# OUTPUT
#   3_subsequence/asvs.subsequence_matches.tsv
#   3_subsequence/asvs.subsequence_resolution.tsv
#   3_subsequence/asvs.subsequence_manual_review.tsv
#   3_subsequence/asvs.subsequence_resolved.fasta
#   3_subsequence/asv_id_to_sequence_resolved.tsv
#   3_subsequence/sample_by_asv_resolved.tsv
#   3_subsequence/asv_by_sample_resolved.tsv
#   3_subsequence/asvs.subsequence_by_sample.tsv
#   3_subsequence/asvs.subsequence_summary.tsv
#   3_subsequence/logs/parameters.tsv
#   3_subsequence/logs/script_sha256.txt
#
# SETTINGS
#   Reference length = 313 bp; compare every strictly shorter/longer ASV pair.
#   Matches require 100% nucleotide identity, without mismatches or gaps.
#   Report terminal and internal matches independently of length class.
#   Child (length - 313) divisible by 3: retain and flag matches for review.
#   Other children: auto-merge only with one unique longer parent and exclusively
#   terminal occurrence(s). Ambiguous/internal matches require manual review.
#   Preserve ASV IDs and sample-specific counts; verify read conservation.
#   Existing output directory is protected against overwriting.
#   MACSE filtering is deferred until after taxonomic assignment.
# ==============================================================================
set -euo pipefail

WORK_DIR="/storage2/jihoonkim/Arielina_eDNA/analyses/1_20261103"
INPUT_DIR="${WORK_DIR}/2_dada2"
INPUT_FASTA="${INPUT_DIR}/ASVs.nochim.fasta"
INPUT_MAP="${INPUT_DIR}/ASV_ID_to_sequence.tsv"
INPUT_ABUND="${INPUT_DIR}/sample_by_ASV.tsv"
OUT_DIR="${WORK_DIR}/3_subsequence"
LOG_DIR="${OUT_DIR}/logs"
SCRIPT="${WORK_DIR}/3_subsequence.sh"
REFERENCE_LENGTH=313

# ==============================================================================
# PREFLIGHT
# ==============================================================================

command -v python3 >/dev/null 2>&1 || {
    echo "[ERROR] Missing executable: python3" >&2
    exit 1
}

[[ -d "$WORK_DIR" && -d "$INPUT_DIR" ]] || {
    echo "[ERROR] Work or DADA2 directory does not exist" >&2
    exit 1
}

for file in "$INPUT_FASTA" "$INPUT_MAP" "$INPUT_ABUND"; do
    [[ -s "$file" ]] || {
        echo "[ERROR] Missing or empty input: $file" >&2
        exit 1
    }
done

[[ ! -e "$OUT_DIR" ]] || {
    echo "[ERROR] Existing output protected: $OUT_DIR" >&2
    exit 1
}

mkdir -p "$LOG_DIR"

{
    printf 'parameter\tvalue\n'
    printf 'date\t%s\n' "$(date -Is)"
    printf 'work_dir\t%s\n' "$WORK_DIR"
    printf 'input_fasta\t%s\n' "$INPUT_FASTA"
    printf 'input_map\t%s\n' "$INPUT_MAP"
    printf 'input_abundance\t%s\n' "$INPUT_ABUND"
    printf 'python_version\t%s\n' "$(python3 --version 2>&1)"
    printf 'reference_length_bp\t%s\n' "$REFERENCE_LENGTH"
    printf 'minimum_identity_pct\t100\n'
    printf 'require_unique_parent\ttrue\n'
    printf 'require_exclusively_terminal_match\ttrue\n'
    printf 'auto_merge_non_313_mod3_children_only\ttrue\n'
    printf 'preserve_original_asv_ids\ttrue\n'
    printf 'per_sample_read_conservation\ttrue\n'
} > "${LOG_DIR}/parameters.tsv"

echo "[INFO] Start: $(date)"
echo "[INFO] Work directory: $WORK_DIR"
echo "[INFO] Reference length: ${REFERENCE_LENGTH} bp"
echo "[INFO] Search: all shorter/longer ASV pairs; exact identity"
echo "[INFO] Automatic merge: unique terminal match, non-3n child length class"

# ==============================================================================
# STEP 1. VALIDATE, MATCH, RESOLVE AND SUMMARIZE ASVs
# ==============================================================================

python3 - "$INPUT_FASTA" "$INPUT_MAP" "$INPUT_ABUND" "$OUT_DIR" <<'PY'
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

fasta_path, map_path, abund_path, out_path = sys.argv[1:]
out = Path(out_path)
REFERENCE_LENGTH = 313


def read_tsv(path, required):
    with open(path, newline='') as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        fields = reader.fieldnames
        if not fields or len(set(fields)) != len(fields):
            raise ValueError('Invalid/duplicate TSV headers: ' + path)
        if not set(required).issubset(fields):
            raise ValueError('Missing columns in {}: {}'.format(path, set(required) - set(fields)))
        return list(reader), fields


def write_tsv(filename, columns, records):
    with (out / filename).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(records)


def read_fasta(path):
    sequences = {}
    name, chunks = None, []
    with open(path) as stream:
        for raw in stream:
            line = raw.strip()
            if not line:
                continue
            if line.startswith('>'):
                if name is not None:
                    if name in sequences:
                        raise ValueError('Duplicate FASTA identifier: ' + name)
                    sequences[name] = ''.join(chunks).upper()
                name, chunks = line[1:].split()[0], []
            else:
                if name is None:
                    raise ValueError('Sequence line before FASTA header')
                chunks.append(line)
        if name is not None:
            if name in sequences:
                raise ValueError('Duplicate FASTA identifier: ' + name)
            sequences[name] = ''.join(chunks).upper()
    return sequences


def as_count(value, description):
    try:
        number = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError('Invalid count for {}: {}'.format(description, value)) from exc
    if number < 0:
        raise ValueError('Negative count: ' + description)
    return number


# ==============================================================================
# STEP 1. Validate current input sequences and sample-level abundances.
# ==============================================================================

fasta_seqs = read_fasta(fasta_path)
metadata, _ = read_tsv(map_path, {'ASV_ID', 'sequence', 'total_reads'})
sample_rows, sample_fields = read_tsv(abund_path, set())
if not fasta_seqs or not metadata or not sample_rows:
    raise ValueError('Input file is empty')

all_ids = [row['ASV_ID'] for row in metadata]
if len(all_ids) != len(set(all_ids)):
    raise ValueError('Duplicate ASV identifier in metadata')
meta = {row['ASV_ID']: row for row in metadata}
if not set(fasta_seqs).issubset(meta):
    raise ValueError('Input FASTA contains an ASV absent from sequence metadata')

# DADA2 sample column can be named Sample, sample, or sample_ID.
# Its identity is determined by exclusion of all ASV identifier columns.
non_asv_columns = [col for col in sample_fields if col not in meta]
if len(non_asv_columns) != 1 or set(sample_fields) != set(meta) | set(non_asv_columns):
    raise ValueError('Expected every original ASV column and exactly one sample-ID column')
sample_column = non_asv_columns[0]
samples = [row[sample_column] for row in sample_rows]
if any(not name or not name.strip() for name in samples) or len(samples) != len(set(samples)):
    raise ValueError('Blank or duplicate sample name')

seq = {}
original_totals = {}
for aid in all_ids:
    row = meta[aid]
    sequence = row['sequence'].upper()
    if not sequence or set(sequence) - set('ACGT'):
        raise ValueError('Invalid nucleotide sequence: ' + aid)
    if 'length_bp' in row and as_count(row['length_bp'], aid) != len(sequence):
        raise ValueError('Sequence length mismatch: ' + aid)
    seq[aid] = sequence
    original_totals[aid] = as_count(row['total_reads'], aid)

input_ids = [aid for aid in all_ids if aid in fasta_seqs]
if len({seq[aid] for aid in input_ids}) != len(input_ids):
    raise ValueError('Input FASTA has duplicate nucleotide sequences under different IDs')
for aid in input_ids:
    if fasta_seqs[aid] != seq[aid]:
        raise ValueError('Input FASTA and sequence metadata disagree: ' + aid)

counts = {}
for row in sample_rows:
    sample = row[sample_column]
    counts[sample] = {aid: as_count(row[aid], sample + '/' + aid) for aid in all_ids}
for aid in all_ids:
    if sum(counts[sample][aid] for sample in samples) != original_totals[aid]:
        raise ValueError('Sample abundance total mismatch: ' + aid)

input_reads = sum(original_totals[aid] for aid in input_ids)


# ==============================================================================
# STEP 2. Exhaustively examine exact containment for every shorter/longer pair.
# ==============================================================================

matches = []
parents_by_child = defaultdict(list)
compared_pairs = 0
for child in input_ids:
    child_seq = seq[child]
    for parent in input_ids:
        parent_seq = seq[parent]
        if len(parent_seq) <= len(child_seq):
            continue
        compared_pairs += 1
        positions = []
        index = parent_seq.find(child_seq)
        while index != -1:
            positions.append(index)
            index = parent_seq.find(child_seq, index + 1)
        if not positions:
            continue
        location_types = []
        for index in positions:
            match_type = ('5prime_terminal' if index == 0 else
                          '3prime_terminal' if index + len(child_seq) == len(parent_seq) else
                          'internal')
            location_types.append(match_type)
            matches.append(dict(
                child_ASV=child, child_length_bp=len(child_seq),
                child_reads=original_totals[child],
                child_313_mod3=(len(child_seq) - REFERENCE_LENGTH) % 3,
                parent_ASV=parent, parent_length_bp=len(parent_seq),
                parent_original_reads=original_totals[parent],
                length_difference_bp=len(parent_seq) - len(child_seq),
                match_type=match_type, start_1based=index + 1,
                end_1based=index + len(child_seq), identity_pct='100.0'
            ))
        parents_by_child[child].append((parent, location_types))

match_cols = [
    'child_ASV', 'child_length_bp', 'child_reads', 'child_313_mod3',
    'parent_ASV', 'parent_length_bp', 'parent_original_reads',
    'length_difference_bp', 'match_type', 'start_1based', 'end_1based', 'identity_pct'
]
write_tsv('asvs.subsequence_matches.tsv', match_cols, matches)


# ==============================================================================
# STEP 3. Apply conservative, reference-length-based merge decisions.
# ==============================================================================

resolution = {}
merge_to = {}
for child in input_ids:
    parents = parents_by_child[child]
    n_terminal = sum(any(x.endswith('terminal') for x in kinds) for _, kinds in parents)
    n_internal = sum('internal' in kinds for _, kinds in parents)
    if not parents:
        status, parent = 'retained_no_exact_parent', ''
    elif (len(seq[child]) - REFERENCE_LENGTH) % 3 == 0:
        status, parent = 'manual_review_313_frame_class', ''
    elif len(parents) > 1:
        status, parent = 'manual_review_multiple_parents', ''
    elif n_internal or not n_terminal:
        status, parent = 'manual_review_internal_match', ''
    else:
        status, parent = 'merged_unique_terminal', parents[0][0]
        merge_to[child] = parent
    resolution[child] = dict(
        status=status, direct_parent_ASV=parent,
        matching_parents=len(parents), terminal_parents=n_terminal,
        internal_parents=n_internal,
        matching_parent_ASVs=','.join(p for p, _ in parents)
    )


def destination(aid):
    visited = set()
    while aid in merge_to:
        if aid in visited:
            raise ValueError('Cyclic parent assignments')
        visited.add(aid)
        aid = merge_to[aid]
    return aid


final_destination = {aid: destination(aid) for aid in input_ids}
survivors = [aid for aid in input_ids if aid not in merge_to]
if set(final_destination.values()) != set(survivors):
    raise ValueError('ASV destination mapping is inconsistent')

resolution_rows = []
for aid in input_ids:
    resolution_rows.append(dict(
        ASV_ID=aid, length_bp=len(seq[aid]), original_reads=original_totals[aid],
        child_313_mod3=(len(seq[aid]) - REFERENCE_LENGTH) % 3,
        **resolution[aid], final_destination_ASV=final_destination[aid]
    ))
resolution_cols = [
    'ASV_ID', 'length_bp', 'original_reads', 'child_313_mod3', 'status',
    'direct_parent_ASV', 'final_destination_ASV', 'matching_parents',
    'terminal_parents', 'internal_parents', 'matching_parent_ASVs'
]
write_tsv('asvs.subsequence_resolution.tsv', resolution_cols, resolution_rows)
manual_rows = [row for row in resolution_rows if row['status'].startswith('manual_review_')]
write_tsv('asvs.subsequence_manual_review.tsv', resolution_cols, manual_rows)


# ==============================================================================
# STEP 4. Aggregate reads into surviving ASVs, separately for each sample.
# ==============================================================================

resolved = {sample: {aid: 0 for aid in survivors} for sample in samples}
for sample in samples:
    for aid in input_ids:
        resolved[sample][final_destination[aid]] += counts[sample][aid]

final_totals = {aid: sum(resolved[sample][aid] for sample in samples) for aid in survivors}
output_reads = sum(final_totals.values())
if output_reads != input_reads:
    raise ValueError('Global read conservation FAILED')

sample_qc = []
for sample in samples:
    before = sum(counts[sample][aid] for aid in input_ids)
    after = sum(resolved[sample].values())
    if before != after:
        raise ValueError('Per-sample read conservation FAILED: ' + sample)
    sample_qc.append(dict(
        sample=sample, input_reads=before,
        transferred_child_reads=sum(counts[sample][aid] for aid in merge_to),
        output_reads=after,
        input_detected_ASVs=sum(counts[sample][aid] > 0 for aid in input_ids),
        output_detected_ASVs=sum(resolved[sample][aid] > 0 for aid in survivors)
    ))
write_tsv('asvs.subsequence_by_sample.tsv', [
    'sample', 'input_reads', 'transferred_child_reads', 'output_reads',
    'input_detected_ASVs', 'output_detected_ASVs'
], sample_qc)

with (out / 'asvs.subsequence_resolved.fasta').open('w') as stream:
    for aid in survivors:
        stream.write('>{}\n{}\n'.format(aid, seq[aid]))

resolved_map = []
for aid in survivors:
    children = [c for c in merge_to if final_destination[c] == aid]
    resolved_map.append(dict(
        ASV_ID=aid, sequence=seq[aid], length_bp=len(seq[aid]),
        total_reads=final_totals[aid], original_reads=original_totals[aid],
        transferred_reads=final_totals[aid] - original_totals[aid],
        merged_children=','.join(children)
    ))
write_tsv('asv_id_to_sequence_resolved.tsv', [
    'ASV_ID', 'sequence', 'length_bp', 'total_reads',
    'original_reads', 'transferred_reads', 'merged_children'
], resolved_map)
write_tsv('sample_by_asv_resolved.tsv', ['Sample'] + survivors, [
    dict(Sample=sample, **resolved[sample]) for sample in samples
])
write_tsv('asv_by_sample_resolved.tsv', ['ASV_ID'] + samples + ['Total_reads'], [
    dict(ASV_ID=aid, **{sample: resolved[sample][aid] for sample in samples},
         Total_reads=final_totals[aid]) for aid in survivors
])


# ==============================================================================
# STEP 5. Summarize this stage alone and verify complete accounting.
# ==============================================================================

statuses = Counter(row['status'] for row in resolution_rows)
summary = [
    ('input_asvs', len(input_ids)),
    ('input_reads', input_reads),
    ('examined_child_asvs', len(input_ids)),
    ('examined_shorter_longer_pairs', compared_pairs),
    ('exact_match_rows', len(matches)),
    ('exact_match_child_asvs', sum(bool(parents_by_child[c]) for c in input_ids)),
    ('manual_review_313_frame_class', statuses['manual_review_313_frame_class']),
    ('manual_review_multiple_parents', statuses['manual_review_multiple_parents']),
    ('manual_review_internal_match', statuses['manual_review_internal_match']),
    ('manual_review_total_asvs', len(manual_rows)),
    ('auto_merged_child_asvs', len(merge_to)),
    ('auto_merged_child_reads', sum(original_totals[c] for c in merge_to)),
    ('output_asvs', len(survivors)),
    ('output_reads', output_reads),
]
if len(resolution_rows) != len(input_ids) or len(merge_to) + len(survivors) != len(input_ids):
    raise ValueError('ASV accounting FAILED')
write_tsv('asvs.subsequence_summary.tsv', ['metric', 'value'], [
    dict(metric=metric, value=value) for metric, value in summary
])

print('=' * 68)
print('ASV SUBSEQUENCE SEARCH AND RESOLUTION')
print('=' * 68)
for metric, value in summary:
    print('{:<39s} {:>12,d}'.format(metric, value))
print('\nAUTOMATIC TERMINAL MERGES')
if merge_to:
    for child, parent in merge_to.items():
        print('  {} ({} bp, {} reads) -> {} ({} bp), final={}'.format(
            child, len(seq[child]), original_totals[child],
            parent, len(seq[parent]), final_destination[child]
        ))
else:
    print('  None')
print('\nMANUAL REVIEW ASVS')
if manual_rows:
    for row in manual_rows:
        print('  {} ({} bp): {}; parents={}'.format(
            row['ASV_ID'], row['length_bp'],
            row['status'], row['matching_parent_ASVs']
        ))
else:
    print('  None')
print('\nSAMPLE-LEVEL READ CONSERVATION')
for row in sample_qc:
    print('  {}: {:,} -> {:,}; transferred={:,}'.format(
        row['sample'], row['input_reads'], row['output_reads'],
        row['transferred_child_reads']
    ))
print('\n[QC PASS] All ASVs accounted for; all sample read totals conserved')
PY

# ==============================================================================
# STEP 2. VALIDATE OUTPUT FILES AND RECORD SCRIPT CHECKSUM
# ==============================================================================

for name in \
    asvs.subsequence_matches.tsv \
    asvs.subsequence_resolution.tsv \
    asvs.subsequence_manual_review.tsv \
    asvs.subsequence_resolved.fasta \
    asv_id_to_sequence_resolved.tsv \
    sample_by_asv_resolved.tsv \
    asv_by_sample_resolved.tsv \
    asvs.subsequence_by_sample.tsv \
    asvs.subsequence_summary.tsv; do
    [[ -s "${OUT_DIR}/${name}" ]] || {
        echo "[ERROR] Missing or empty output: ${OUT_DIR}/${name}" >&2
        exit 1
    }
done

if [[ -f "$SCRIPT" ]]; then
    sha256sum "$SCRIPT" > "${LOG_DIR}/script_sha256.txt"
fi

echo "[QC PASS] Subsequence outputs validated"
echo "[DONE] Subsequence resolution completed: $(date)"
echo "[OUTPUT] ${OUT_DIR}/asvs.subsequence_resolved.fasta"
echo "[OUTPUT] ${OUT_DIR}/asvs.subsequence_resolution.tsv"
echo "[OUTPUT] ${OUT_DIR}/asvs.subsequence_summary.tsv"

#!/usr/bin/env bash
#SBATCH --job-name=1_trimming
#SBATCH --output=1_trimming_%j.out
#SBATCH --error=1_trimming_%j.err
#SBATCH --partition=Node6
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G
#SBATCH --time=30-00:00:00

# ==============================================================================
# PURPOSE
#   Trim adapters and COI primers from paired-end reads for DADA2.
#   Require at least one 5' primer per pair, but apply length/N filters to
#   both reads independently (reject a pair if either read fails).
#
# INPUT
#   /storage2/jihoonkim/Arielina_eDNA/raw_data/*_R[12]_001.fastq.gz
#
# OUTPUT
#   1_trimming/1_adapter_trim/   Adapter-trimmed FASTQ.gz (intermediate)
#   1_trimming/2_5prime_trim/    5'-primer-trimmed FASTQ.gz (intermediate)
#   1_trimming/3_primer_trim/    Final paired FASTQ for DADA2
#   1_trimming/logs/             Reports, parameters and QC tables
#
# SETTINGS
#   Adapter: overlap 15 bp, error rate 0.20, no indels.
#   5' primer: overlap 17 bp, error rate 0.20, no indels, at least one match.
#   3' primer: overlap 15 bp, error rate 0.20, no indels, optional match.
#   Final filters: minimum length 50 bp and zero Ns in EACH read.
#   Original primer sequences and wildcard handling are unchanged.
# ==============================================================================
set -euo pipefail

WORK_DIR="/storage2/jihoonkim/Arielina_eDNA/analyses/1_20261103"
RAW_DIR="/storage2/jihoonkim/Arielina_eDNA/raw_data"
THREADS="${SLURM_CPUS_PER_TASK:-12}"

OUT_DIR="${WORK_DIR}/1_trimming"
ADAPTER_DIR="${OUT_DIR}/1_adapter_trim"
FIVEP_DIR="${OUT_DIR}/2_5prime_trim"
FINAL_DIR="${OUT_DIR}/3_primer_trim"
LOG_DIR="${OUT_DIR}/logs"
MAP="${LOG_DIR}/sample_name_map.tsv"
SCRIPT="${WORK_DIR}/1_coi313_trimming.sh"

ADAPT_R1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
ADAPT_R2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
FWD="GGWACWGGWTGAACWGTWTAYCCYCC"
REV_I="TAIACYTCIGGRTGICCRAARAAYCA"
REV_N="${REV_I//I/N}"

revcomp_iupac() {
    printf '%s' "$1" | rev | tr \
        'ACGTRYSWKMBDHVN' 'TGCAYRSWMKVHDBN'
}
REV_RC_N="$(revcomp_iupac "$REV_N")"
RC_FWD="$(revcomp_iupac "$FWD")"

# ==============================================================================
# PREFLIGHT
# ==============================================================================

source /home/JihoonKim/miniconda3/etc/profile.d/conda.sh
conda activate /home/JihoonKim/miniconda3/envs/cutadapt

for exe in cutadapt python3; do
    command -v "$exe" >/dev/null 2>&1 || {
        echo "[ERROR] Missing executable: $exe" >&2
        exit 1
    }
done

[[ -d "$WORK_DIR" && -d "$RAW_DIR" ]] || {
    echo "[ERROR] Work or raw-data directory does not exist" >&2
    exit 1
}
[[ ! -e "$OUT_DIR" ]] || {
    echo "[ERROR] Existing output protected: $OUT_DIR" >&2
    exit 1
}

mkdir -p "$ADAPTER_DIR" "$FIVEP_DIR" "$FINAL_DIR" "$LOG_DIR"

# ==============================================================================
# STEP 1. VALIDATE INPUT PAIRS AND RECORD SETTINGS
# ==============================================================================

echo "[STEP 1] Discover and validate paired FASTQ files"

python3 - "$RAW_DIR" "$MAP" <<'PY'
import csv
import re
import sys
from pathlib import Path

raw_dir, map_path = map(Path, sys.argv[1:])
r1_files = sorted(raw_dir.glob('*_R1_001.fastq.gz'))
r2_files = sorted(raw_dir.glob('*_R2_001.fastq.gz'))

if not r1_files or len(r1_files) != len(r2_files):
    raise SystemExit('[ERROR] Missing R1 files or unequal R1/R2 counts')

rows = []
seen = set()
for r1 in r1_files:
    original = r1.name.removesuffix('_R1_001.fastq.gz')
    r2 = raw_dir / f'{original}_R2_001.fastq.gz'
    if not r1.is_file() or not r2.is_file() or not r1.stat().st_size or not r2.stat().st_size:
        raise SystemExit(f'[ERROR] Missing/empty pair: {original}')
    match = re.fullmatch(r'(.+)_E_Lib_S\d+_L\d{3}', original)
    if match is None:
        match = re.fullmatch(r'(.+)_S\d+_L\d{3}', original)
    sample = match.group(1) if match else original
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', sample) or sample in seen:
        raise SystemExit(f'[ERROR] Invalid/duplicate sample ID: {sample}')
    seen.add(sample)
    rows.append((original, sample))

with map_path.open('w', newline='') as f:
    writer = csv.writer(f, delimiter='\t', lineterminator='\n')
    writer.writerow(['original_sample_ID', 'sample_ID'])
    writer.writerows(rows)

print(f'[QC PASS] {len(rows)} unique paired-end samples')
for original, sample in rows:
    print(f'[QC] {original} -> {sample}')
PY

{
    printf 'parameter\tvalue\n'
    printf 'date\t%s\n' "$(date -Is)"
    printf 'work_dir\t%s\nraw_dir\t%s\n' "$WORK_DIR" "$RAW_DIR"
    printf 'cutadapt_version\t%s\n' "$(cutadapt --version)"
    printf 'adapter_R1\t%s\nadapter_R2\t%s\n' "$ADAPT_R1" "$ADAPT_R2"
    printf 'forward_primer\t%s\nreverse_primer_inosine\t%s\n' "$FWD" "$REV_I"
    printf 'reverse_primer_N_mask\t%s\n' "$REV_N"
    printf 'adapter_min_overlap\t15\nadapter_error_rate\t0.20\n'
    printf 'five_prime_type\tnon-internal_X\n'
    printf 'five_prime_min_overlap\t17\nfive_prime_error_rate\t0.20\n'
    printf 'three_prime_min_overlap\t15\nthree_prime_error_rate\t0.20\n'
    printf 'allow_indels\tfalse\nmatch_read_wildcards\ttrue\n'
    printf 'five_prime_discard_untrimmed\ttrue\n'
    printf 'five_prime_pair_filter\tboth\n'
    printf 'three_prime_discard_untrimmed\tfalse\n'
    printf 'final_pair_filter\tany\n'
    printf 'minimum_length\t50\nmax_n\t0\n'
    printf 'trim_n\ttrue\n'
    printf 'trunc_len_diagnostic_R1\t240\ntrunc_len_diagnostic_R2\t200\n'
    printf 'threads\t%s\n' "$THREADS"
} > "${LOG_DIR}/parameters.tsv"

echo "[INFO] Cutadapt $(cutadapt --version); threads=$THREADS"

# ==============================================================================
# STEP 2. ADAPTER -> 5' PRIMER -> 3' PRIMER + FINAL QUALITY FILTER
# ==============================================================================

{
    IFS= read -r header
    while IFS=$'\t' read -r original sample; do
        echo "[SAMPLE] $sample"

        raw_r1="${RAW_DIR}/${original}_R1_001.fastq.gz"
        raw_r2="${RAW_DIR}/${original}_R2_001.fastq.gz"
        adapter_r1="${ADAPTER_DIR}/${sample}_R1_001.fastq.gz"
        adapter_r2="${ADAPTER_DIR}/${sample}_R2_001.fastq.gz"
        fivep_r1="${FIVEP_DIR}/${sample}_R1_001.fastq.gz"
        fivep_r2="${FIVEP_DIR}/${sample}_R2_001.fastq.gz"
        final_r1="${FINAL_DIR}/${sample}_R1_001.fastq"
        final_r2="${FINAL_DIR}/${sample}_R2_001.fastq"

        echo "[TRIM 1/3] $sample: sequencing adapters"
        cutadapt -j "$THREADS" \
            -a "$ADAPT_R1" -A "$ADAPT_R2" \
            -O 15 -e 0.20 --no-indels \
            --trim-n --minimum-length 50 --pair-filter=any \
            -o "$adapter_r1" -p "$adapter_r2" \
            "$raw_r1" "$raw_r2" \
            > "${LOG_DIR}/${sample}_adapter.txt"

        echo "[TRIM 2/3] $sample: 5' primers (one match sufficient)"
        cutadapt -j "$THREADS" \
            -g "X${FWD};e=0.20" -G "X${REV_I};e=0.20" \
            -O 17 --no-indels --match-read-wildcards \
            --discard-untrimmed --pair-filter=both --trim-n \
            -o "$fivep_r1" -p "$fivep_r2" \
            "$adapter_r1" "$adapter_r2" \
            > "${LOG_DIR}/${sample}_5prime.txt"

        echo "[TRIM 3/3] $sample: 3' primers and per-read filters"
        cutadapt -j "$THREADS" \
            -a "$REV_RC_N" -A "$RC_FWD" \
            -O 15 -e 0.20 --no-indels --match-read-wildcards \
            --trim-n --max-n 0 --minimum-length 50 --pair-filter=any \
            -o "$final_r1" -p "$final_r2" \
            "$fivep_r1" "$fivep_r2" \
            > "${LOG_DIR}/${sample}_3prime.txt"

        echo "[QC] $sample: all three Cutadapt passes completed"
    done
} < "$MAP"

# ==============================================================================
# STEP 3. VALIDATE FINAL FASTQ PAIRS AND SUMMARIZE READ RETENTION
# ==============================================================================

echo "[STEP 3] Final FASTQ integrity, pairing and length QC"

python3 - "$MAP" "$FINAL_DIR" "$LOG_DIR" <<'PY'
import csv
import re
import sys
from collections import Counter
from itertools import zip_longest
from pathlib import Path

mapping_file, final_dir, log_dir = map(Path, sys.argv[1:])


def get_count(path, label):
    text = path.read_text()
    match = re.search(r'^\s*' + re.escape(label) + r':\s*([\d,]+)', text, re.M)
    if not match:
        raise RuntimeError(f'Missing Cutadapt counter {label} in {path}')
    return int(match.group(1).replace(',', ''))


def reads(path):
    with path.open(encoding='ascii') as f:
        while (header := f.readline()):
            sequence = f.readline().rstrip('\r\n')
            plus = f.readline()
            quality = f.readline().rstrip('\r\n')
            if not header.startswith('@') or not plus.startswith('+') or len(sequence) != len(quality):
                raise RuntimeError(f'Invalid FASTQ record: {path}')
            name = header.split()[0][1:]
            if name.endswith(('/1', '/2')):
                name = name[:-2]
            yield name, sequence


with mapping_file.open(newline='') as f:
    samples = [row['sample_ID'] for row in csv.DictReader(f, delimiter='\t')]

count_rows = []
length_rows = []
for sample in samples:
    raw = get_count(log_dir / f'{sample}_adapter.txt', 'Total read pairs processed')
    adapter = get_count(log_dir / f'{sample}_adapter.txt', 'Pairs written (passing filters)')
    fiveprime = get_count(log_dir / f'{sample}_5prime.txt', 'Pairs written (passing filters)')
    final = get_count(log_dir / f'{sample}_3prime.txt', 'Pairs written (passing filters)')
    if not raw >= adapter >= fiveprime >= final > 0:
        raise RuntimeError(f'Unexpected read retention: {sample}')

    r1 = final_dir / f'{sample}_R1_001.fastq'
    r2 = final_dir / f'{sample}_R2_001.fastq'
    h1, h2 = Counter(), Counter()
    observed = too_short = 0

    for mate1, mate2 in zip_longest(reads(r1), reads(r2)):
        if mate1 is None or mate2 is None or mate1[0] != mate2[0]:
            raise RuntimeError(f'Unsynchronized reads: {sample}, pair {observed + 1}')
        n1, n2 = len(mate1[1]), len(mate2[1])
        if n1 < 50 or n2 < 50 or 'N' in mate1[1].upper() or 'N' in mate2[1].upper():
            raise RuntimeError(f'Final length/N filter failed: {sample}, pair {observed + 1}')
        observed += 1
        h1[n1] += 1
        h2[n2] += 1
        too_short += int(n1 < 240 or n2 < 200)

    if observed != final:
        raise RuntimeError(f'Cutadapt/FASTQ count mismatch: {sample}')
    count_rows.append([sample, raw, adapter, fiveprime, final,
                       round(100 * final / raw, 4), too_short])
    for direction, histogram in (('R1', h1), ('R2', h2)):
        length_rows.extend([sample, direction, length, count]
                           for length, count in sorted(histogram.items()))
    print(f'[QC PASS] {sample}: {final:,} pairs; below DADA2 truncLen: {too_short:,}')

with (log_dir / 'read_counts_summary.tsv').open('x', newline='') as f:
    writer = csv.writer(f, delimiter='\t', lineterminator='\n')
    writer.writerow(['Sample_ID', 'Step1_Raw', 'Step2_Adapter', 'Step3_5prime',
                     'Step4_Primer_Final', 'retained_pct_raw', 'pairs_below_truncLen'])
    writer.writerows(count_rows)

with (log_dir / 'final_length_distribution.tsv').open('x', newline='') as f:
    writer = csv.writer(f, delimiter='\t', lineterminator='\n')
    writer.writerow(['sample', 'read', 'length_bp', 'count'])
    writer.writerows(length_rows)

print(f'[QC PASS] Total final pairs: {sum(row[4] for row in count_rows):,}')
PY

if [[ -f "$SCRIPT" ]]; then
    sha256sum "$SCRIPT" > "${LOG_DIR}/script_sha256.txt"
fi

echo "[DONE] COI313 trimming completed: $(date)"
echo "[OUTPUT] ${FINAL_DIR}"
echo "[OUTPUT] ${LOG_DIR}/read_counts_summary.tsv"

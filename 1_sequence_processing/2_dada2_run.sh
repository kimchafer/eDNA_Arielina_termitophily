#!/usr/bin/env bash
#SBATCH --job-name=2_dada2
#SBATCH --output=2_dada2_%j.out
#SBATCH --error=2_dada2_%j.err
#SBATCH --partition=Node4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --time=30-00:00:00

# ==============================================================================
# PURPOSE
#   Run paired-end COI ASV inference using DADA2.
#   Execute quality filtering, error learning, pseudo-pooling, paired-end
#   merging, chimera removal, sequence-length selection and read tracking.
#
# INPUT
#   /storage2/jihoonkim/Arielina_eDNA/analyses/1_20261103/1_trimming/3_primer_trim/*_R[12]_001.fastq
#   /storage2/jihoonkim/Arielina_eDNA/analyses/1_20261103/1_trimming/logs/read_counts_summary.tsv
#   /storage2/jihoonkim/Arielina_eDNA/analyses/1_20261103/2_dada2.r
#
# OUTPUT
#   2_dada2/                         DADA2 ASV tables, FASTA and intermediate RDS
#   2_dada2/filtered/                Filtered paired-end reads
#   2_dada2/plots/                   Read-quality and error-model plots
#   2_dada2/logs/                    Parameters, QC, tracking and script checksum
#
# SETTINGS
#   R environment: dada2_coi313 (existing Conda environment).
#   CPU: 12; memory: 64 GB; walltime: 30 days; partition: Node4.
#   DADA2 parameters and ASV selection are defined in 2_dada2.r.
# ==============================================================================
set -euo pipefail

WORK_DIR="/storage2/jihoonkim/Arielina_eDNA/analyses/1_20261103"
IN_DIR="${WORK_DIR}/1_trimming/3_primer_trim"
OUT_DIR="${WORK_DIR}/2_dada2"
LOG_DIR="${OUT_DIR}/logs"
R_SCRIPT="${WORK_DIR}/2_dada2.r"
SCRIPT="${WORK_DIR}/2_dada2.sh"
THREADS="${SLURM_CPUS_PER_TASK:-12}"

# ==============================================================================
# PREFLIGHT
# ==============================================================================

echo "[INFO] DADA2 preflight: $(date)"

source /home/JihoonKim/miniconda3/etc/profile.d/conda.sh
conda activate dada2_coi313

command -v Rscript >/dev/null 2>&1 || {
    echo "[ERROR] Missing executable: Rscript" >&2
    exit 1
}

[[ -d "$WORK_DIR" && -d "$IN_DIR" ]] || {
    echo "[ERROR] Work or trimmed-read directory does not exist" >&2
    exit 1
}
[[ -s "$R_SCRIPT" ]] || {
    echo "[ERROR] R script missing or empty: $R_SCRIPT" >&2
    exit 1
}
[[ -s "${WORK_DIR}/1_trimming/logs/read_counts_summary.tsv" ]] || {
    echo "[ERROR] Step 1 QC summary missing or empty" >&2
    exit 1
}
[[ ! -e "$OUT_DIR" ]] || {
    echo "[ERROR] Existing output protected: $OUT_DIR" >&2
    exit 1
}

echo "[INFO] Working directory: $WORK_DIR"
echo "[INFO] Rscript: $(command -v Rscript)"
echo "[INFO] DADA2: $(Rscript -e 'cat(as.character(packageVersion("dada2")))')"
echo "[INFO] Threads: $THREADS"

# ==============================================================================
# STEP 1. EXECUTE DADA2 ANALYSIS
# ==============================================================================

export COI_WORKDIR="$WORK_DIR"

Rscript --vanilla "$R_SCRIPT"

# ==============================================================================
# STEP 2. RECORD SCRIPT CHECKSUMS AND FINAL OUTPUT
# ==============================================================================

if [[ -d "$LOG_DIR" ]]; then
    if [[ -f "$SCRIPT" ]]; then
        sha256sum "$SCRIPT" > "${LOG_DIR}/script_sha256.txt"
    fi
    sha256sum "$R_SCRIPT" >> "${LOG_DIR}/script_sha256.txt"
fi

echo "[DONE] DADA2 completed: $(date)"
echo "[OUTPUT] $OUT_DIR"
echo "[OUTPUT] ${LOG_DIR}/read_tracking.tsv"

#!/usr/bin/env bash
#SBATCH --job-name=5_lulu
#SBATCH --output=5_lulu_%j.out
#SBATCH --error=5_lulu_%j.err
#SBATCH --partition=Node4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=30-00:00:00

# ==============================================================================
# PURPOSE
#   Curate Swarm MOTUs using VSEARCH sequence matches and LULU.
#   Reconstruct curated MOTU sequences, sample abundances and parent mappings.
#
# INPUT (example: d=2)
#   /storage2/jihoonkim/Arielina_eDNA/analyses/1_20261103/4_swarm_motu/ALL_d2/ALL_d2_MOTU.fasta
#   /storage2/jihoonkim/Arielina_eDNA/analyses/1_20261103/4_swarm_motu/ALL_d2/ALL_d2_MOTU_table.tsv
#
# OUTPUT (example: d=2)
#   5_lulu/ALL_d2/ALL_d2_vsearch_matchlist.tsv   VSEARCH candidate pairs
#   5_lulu/ALL_d2/LULU_curated_MOTU_table.tsv    Curated sample abundances
#   5_lulu/ALL_d2/LULU_curated_MOTU.fasta        Retained representative sequences
#   5_lulu/ALL_d2/LULU_MOTU_map.tsv             MOTU-to-parent mapping
#   5_lulu/ALL_d2/LULU_sample_summary.tsv       Per-sample retention
#   5_lulu/ALL_d2/LULU_zero_samples.tsv         Zero-read samples
#   5_lulu/ALL_d2/LULU_summary.tsv              Global curation summary
#   5_lulu/ALL_d2/LULU_result.rds               Original LULU result object
#   5_lulu/ALL_d2/logs/                         Parameters and script checksum
#
# SETTINGS
#   Swarm d=2 by default; override with first sbatch script argument.
#   VSEARCH global identity >=97%, iddef=1, query coverage >=90%.
#   LULU minimum_match=97%, minimum_ratio_type=min, minimum_ratio=1,
#   minimum_relative_cooccurence=0.95.
#   Preserve sample IDs and read totals; no read removal.
#   Existing output directory is protected against overwriting.
# ==============================================================================
set -euo pipefail

WORK_DIR="/storage2/jihoonkim/Arielina_eDNA/analyses/1_20261103"
D_VAL="${1:-2}"
THREADS="${SLURM_CPUS_PER_TASK:-4}"
MINIMUM_MATCH=97
MINIMUM_RATIO=1
MINIMUM_COOCCURRENCE=0.95
SCRIPT="${WORK_DIR}/5_lulu.sh"

POOL="ALL_d${D_VAL}"
SWARM_DIR="${WORK_DIR}/4_swarm_motu/${POOL}"
FASTA="${SWARM_DIR}/${POOL}_MOTU.fasta"
TABLE="${SWARM_DIR}/${POOL}_MOTU_table.tsv"
OUT_DIR="${WORK_DIR}/5_lulu/${POOL}"
LOG_DIR="${OUT_DIR}/logs"
MATCHLIST="${OUT_DIR}/${POOL}_vsearch_matchlist.tsv"

# ==============================================================================
# PREFLIGHT
# ==============================================================================

[[ "$D_VAL" =~ ^[1-9][0-9]*$ ]] || {
    echo "[ERROR] Swarm d must be a positive integer: $D_VAL" >&2
    exit 1
}
[[ -d "$WORK_DIR" && -d "$SWARM_DIR" ]] || {
    echo "[ERROR] Missing work or Swarm directory" >&2
    exit 1
}
for file in "$FASTA" "$TABLE"; do
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
for exe in vsearch Rscript; do
    command -v "$exe" >/dev/null 2>&1 || {
        echo "[ERROR] Missing executable: $exe" >&2
        exit 1
    }
done
Rscript --vanilla -e '
  if (!requireNamespace("lulu", quietly=TRUE)) stop("R package lulu is missing")
  cat("[INFO] LULU version: ", as.character(packageVersion("lulu")), "\n", sep="")
'

mkdir -p "$LOG_DIR"
cd "$OUT_DIR"
{
    printf 'parameter\tvalue\n'
    printf 'date\t%s\n' "$(date -Is)"
    printf 'work_dir\t%s\n' "$WORK_DIR"
    printf 'swarm_d\t%s\n' "$D_VAL"
    printf 'input_fasta\t%s\n' "$FASTA"
    printf 'input_table\t%s\n' "$TABLE"
    printf 'vsearch_version\t%s\n' "$(vsearch --version 2>&1 | head -n 1)"
    printf 'lulu_version\t%s\n' "$(Rscript --vanilla -e 'cat(as.character(packageVersion("lulu")))')"
    printf 'vsearch_identity_threshold\t0.97\n'
    printf 'vsearch_iddef\t1\n'
    printf 'vsearch_query_coverage\t0.90\n'
    printf 'vsearch_maxaccepts\t0\n'
    printf 'vsearch_maxrejects\t0\n'
    printf 'lulu_minimum_match\t%s\n' "$MINIMUM_MATCH"
    printf 'lulu_minimum_ratio_type\tmin\n'
    printf 'lulu_minimum_ratio\t%s\n' "$MINIMUM_RATIO"
    printf 'lulu_minimum_relative_cooccurence\t%s\n' "$MINIMUM_COOCCURRENCE"
    printf 'threads\t%s\n' "$THREADS"
} > "${LOG_DIR}/parameters.tsv"

echo "[INFO] Start: $(date)"
echo "[INFO] Work directory: $WORK_DIR"
echo "[INFO] Input Swarm d: $D_VAL"
echo "[INFO] VSEARCH: $(vsearch --version 2>&1 | head -n 1)"
echo "[INFO] VSEARCH minimum identity: ${MINIMUM_MATCH}%"
echo "[INFO] LULU minimum ratio: $MINIMUM_RATIO (min)"
echo "[INFO] LULU minimum relative co-occurrence: $MINIMUM_COOCCURRENCE"

# ==============================================================================
# STEP 1. FIND HIGH-IDENTITY NONSELF MOTU SEQUENCE MATCHES
# ==============================================================================

echo "[STEP 1] VSEARCH global alignment of MOTU representatives"

vsearch \
    --usearch_global "$FASTA" \
    --db "$FASTA" \
    --strand plus \
    --self \
    --id 0.97 \
    --iddef 1 \
    --query_cov 0.90 \
    --maxaccepts 0 \
    --maxrejects 0 \
    --threads "$THREADS" \
    --userout "$MATCHLIST" \
    --userfields query+target+id

[[ -f "$MATCHLIST" ]] || {
    echo "[ERROR] Missing VSEARCH match list: $MATCHLIST" >&2
    exit 1
}

echo "[QC] VSEARCH match records: $(wc -l < "$MATCHLIST")"

# ==============================================================================
# STEP 2. LULU CURATION, PARENT MAPPING AND READ CONSERVATION
# ==============================================================================

echo "[STEP 2] LULU curation and abundance validation"

Rscript --vanilla - \
    "$TABLE" "$FASTA" "$MATCHLIST" "$OUT_DIR" "$D_VAL" \
    "$MINIMUM_MATCH" "$MINIMUM_RATIO" "$MINIMUM_COOCCURRENCE" <<'RSCRIPT'

suppressPackageStartupMessages(library(lulu))

args <- commandArgs(trailingOnly=TRUE)
input_table <- args[1]
input_fasta <- args[2]
match_file <- args[3]
out_dir <- args[4]
d <- as.integer(args[5])
minimum_match <- as.numeric(args[6])
minimum_ratio <- as.numeric(args[7])
minimum_cooccur <- as.numeric(args[8])

read_fasta <- function(path) {
  txt <- readLines(path, warn=FALSE)
  heads <- which(startsWith(txt, ">"))
  if (!length(heads)) stop("No FASTA headers: ", path)
  ids <- sub("^>([^[:space:]]+).*$", "\\1", txt[heads])
  if (anyDuplicated(ids) || any(!nzchar(ids))) stop("Duplicate/empty FASTA IDs")
  ends <- c(heads[-1L]-1L, length(txt))
  seqs <- vapply(seq_along(heads), function(i) {
    ix <- seq.int(heads[i]+1L, ends[i])
    if (!length(ix) || heads[i]+1L > ends[i]) stop("Empty FASTA sequence: ", ids[i])
    paste0(txt[ix], collapse="")
  }, character(1))
  seqs <- toupper(seqs)
  if (any(!grepl("^[ACGT]+$", seqs))) stop("Invalid nucleotide sequences in FASTA")
  setNames(seqs, ids)
}

seqs <- read_fasta(input_fasta)
raw <- read.delim(input_table, check.names=FALSE, stringsAsFactors=FALSE)
if (ncol(raw) < 2L || nrow(raw) < 1L) stop("Invalid MOTU abundance table")
ids <- as.character(raw[[1]])
if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) stop("Invalid MOTU IDs")
if (!setequal(ids, names(seqs))) stop("FASTA and MOTU abundance table have different IDs")
sample_ids <- names(raw)[-1L]
if (anyDuplicated(sample_ids) || any(!nzchar(sample_ids))) stop("Duplicate/empty sample names")

otu_full <- as.data.frame(raw[-1L], check.names=FALSE, stringsAsFactors=FALSE)
rownames(otu_full) <- ids
for (sample in sample_ids) {
  vals <- suppressWarnings(as.numeric(otu_full[[sample]]))
  if (anyNA(vals) || any(!is.finite(vals)) || any(vals < 0) ||
      any(vals != floor(vals))) stop("Invalid read counts in sample: ", sample)
  otu_full[[sample]] <- vals
}
if (any(rowSums(otu_full) <= 0)) stop("Input MOTUs with zero reads")
input_totals <- colSums(otu_full)
input_reads <- sum(input_totals)
input_motus <- nrow(otu_full)

zero_samples <- names(input_totals)[input_totals == 0]
write.table(data.frame(sample=zero_samples),
            file.path(out_dir, "LULU_zero_samples.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
work_samples <- setdiff(sample_ids, zero_samples)
if (!length(work_samples)) stop("All input samples contain zero reads")
otu <- otu_full[, work_samples, drop=FALSE]

matches <- if (file.size(match_file) == 0) {
  data.frame(OTUid=character(), hit=character(), match=numeric())
} else {
  read.delim(match_file, header=FALSE, sep="\t", stringsAsFactors=FALSE,
             col.names=c("OTUid", "hit", "match"))
}
if (anyNA(matches$match) || any(!is.finite(matches$match)) ||
    any(matches$match < 0) || any(matches$match > 100)) {
  stop("Invalid VSEARCH percentage identity")
}
if (!all(matches$OTUid %in% ids) || !all(matches$hit %in% ids)) {
  stop("VSEARCH match list contains unknown MOTU IDs")
}
if (any(matches$OTUid == matches$hit)) stop("Self-hits found in VSEARCH match list")

if (nrow(matches)) {
  result <- lulu::lulu(
    otu, matches,
    minimum_ratio_type="min",
    minimum_ratio=minimum_ratio,
    minimum_match=minimum_match,
    minimum_relative_cooccurence=minimum_cooccur
  )
} else {
  # No match candidates: retain all input MOTUs unchanged.
  result <- list(
    curated_table=otu,
    curated_count=nrow(otu),
    curated_otus=rownames(otu),
    discarded_count=0L,
    discarded_otus=character(),
    otu_map=data.frame(
      total=rowSums(otu), spread=rowSums(otu > 0),
      parent_id=rownames(otu), curated=rep("parent", nrow(otu)),
      row.names=rownames(otu)
    ),
    original_table=otu
  )
}

curated <- as.data.frame(result$curated_table, check.names=FALSE)
map <- as.data.frame(result$otu_map, check.names=FALSE)
retained <- rownames(curated)
if (!length(retained) || anyDuplicated(retained) ||
    !setequal(retained, as.character(result$curated_otus)) ||
    nrow(curated) != result$curated_count ||
    nrow(curated) + result$discarded_count != input_motus) {
  stop("LULU curated MOTU accounting failed")
}
if (!setequal(rownames(map), ids) || !"parent_id" %in% names(map)) {
  stop("LULU parent mapping is missing or inconsistent")
}
if (!all(map$parent_id %in% retained)) {
  stop("LULU map includes a parent absent from the curated table")
}
if (!setequal(colnames(curated), work_samples)) {
  stop("LULU curated sample names disagree")
}
curated <- curated[, work_samples, drop=FALSE]

# Restore samples with zero reads so output shape matches input.
for (s in zero_samples) curated[[s]] <- rep(0, nrow(curated))
curated <- curated[, sample_ids, drop=FALSE]

# Independently reconstruct the curated abundance table from LULU parent mapping.
expected <- rowsum(as.matrix(otu_full),
                   group=map[ids, "parent_id"], reorder=FALSE)
if (!setequal(rownames(expected), rownames(curated))) {
  stop("Curated MOTU IDs do not agree with the parent map")
}
expected <- expected[rownames(curated), sample_ids, drop=FALSE]
if (!identical(dim(expected), dim(as.matrix(curated))) ||
    any(as.matrix(curated) != expected)) {
  stop("Curated read counts are not explained by LULU parent mapping")
}
if (any(colSums(curated) != input_totals)) {
  stop("Per-sample read conservation failed")
}
if (sum(as.matrix(curated)) != input_reads) {
  stop("Global read conservation failed")
}

# Curated output: MOTUs as rows, samples as columns (LULU-compatible).
write.table(
  data.frame(MOTU_ID=rownames(curated), curated, check.names=FALSE),
  file.path(out_dir, "LULU_curated_MOTU_table.tsv"),
  sep="\t", quote=FALSE, row.names=FALSE
)
write.table(
  data.frame(MOTU_ID=rownames(map), map, check.names=FALSE),
  file.path(out_dir, "LULU_MOTU_map.tsv"),
  sep="\t", quote=FALSE, row.names=FALSE
)

with_fasta <- file(file.path(out_dir, "LULU_curated_MOTU.fasta"), "w")
for (id in rownames(curated)) {
  writeLines(c(paste0(">", id), seqs[[id]]), with_fasta)
}
close(with_fasta)
saveRDS(result, file.path(out_dir, "LULU_result.rds"))

sample_qc <- data.frame(
  sample=sample_ids,
  input_reads=as.numeric(input_totals),
  curated_reads=as.numeric(colSums(curated)),
  input_detected_MOTUs=as.numeric(colSums(otu_full > 0)),
  curated_detected_MOTUs=as.numeric(colSums(curated > 0))
)
write.table(sample_qc, file.path(out_dir, "LULU_sample_summary.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

summary <- data.frame(
  swarm_d=d,
  minimum_match=minimum_match,
  minimum_ratio_type="min",
  minimum_ratio=minimum_ratio,
  minimum_relative_cooccurence=minimum_cooccur,
  vsearch_match_rows=nrow(matches),
  input_MOTUs=input_motus,
  retained_MOTUs=nrow(curated),
  merged_MOTUs=input_motus-nrow(curated),
  input_samples=length(sample_ids),
  zero_read_samples=length(zero_samples),
  input_reads=input_reads,
  curated_reads=sum(as.matrix(curated)),
  read_conservation="PASS",
  stringsAsFactors=FALSE
)
write.table(summary, file.path(out_dir, "LULU_summary.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

cat("\n", paste(rep("=", 65), collapse=""), "\n", sep="")
cat("LULU CURATION SUMMARY\n")
cat(paste(rep("=", 65), collapse=""), "\n", sep="")
for (name in names(summary)) {
  cat(sprintf("%-36s %s\n", name, as.character(summary[[name]][1])))
}
cat("\nSAMPLE-LEVEL READ CONSERVATION\n")
for (i in seq_len(nrow(sample_qc))) {
  cat(sprintf("  %s: %s -> %s reads; %s -> %s MOTUs\n",
              sample_qc$sample[i],
              format(sample_qc$input_reads[i], big.mark=",", scientific=FALSE),
              format(sample_qc$curated_reads[i], big.mark=",", scientific=FALSE),
              sample_qc$input_detected_MOTUs[i],
              sample_qc$curated_detected_MOTUs[i]))
}
cat("\n[QC PASS] All MOTUs mapped; curated counts and sample reads conserved\n")
RSCRIPT

# ==============================================================================
# STEP 3. VALIDATE OUTPUT FILES AND RECORD SCRIPT CHECKSUM
# ==============================================================================

for file in \
    "${OUT_DIR}/LULU_curated_MOTU_table.tsv" \
    "${OUT_DIR}/LULU_curated_MOTU.fasta" \
    "${OUT_DIR}/LULU_MOTU_map.tsv" \
    "${OUT_DIR}/LULU_sample_summary.tsv" \
    "${OUT_DIR}/LULU_zero_samples.tsv" \
    "${OUT_DIR}/LULU_summary.tsv" \
    "${OUT_DIR}/LULU_result.rds"; do
    [[ -s "$file" ]] || {
        echo "[ERROR] Empty or missing final output: $file" >&2
        exit 1
    }
done

if [[ -f "$SCRIPT" ]]; then
    sha256sum "$SCRIPT" > "${LOG_DIR}/script_sha256.txt"
fi

echo "[QC PASS] LULU output files validated"
echo "[DONE] LULU curation completed: $(date)"
echo "[OUTPUT] ${OUT_DIR}/LULU_curated_MOTU_table.tsv"
echo "[OUTPUT] ${OUT_DIR}/LULU_curated_MOTU.fasta"
echo "[OUTPUT] ${OUT_DIR}/LULU_MOTU_map.tsv"
echo "[OUTPUT] ${OUT_DIR}/LULU_summary.tsv"

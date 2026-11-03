#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dada2))

# ==============================================================================
# PURPOSE
#   Infer COI amplicon sequence variants from paired trimmed reads.
#   Filter reads, estimate errors with standard LOESS, denoise with pseudo-pooling,
#   merge pairs, remove chimeras, and retain 307-319 bp sequences.
#
# INPUT
#   /storage2/jihoonkim/Arielina_eDNA/analyses/1_20261103/1_trimming/3_primer_trim/*_R[12]_001.fastq
#
# OUTPUT
#   2_dada2/seqtab.prechim.rds      ASVs before chimera removal
#   2_dada2/seqtab.nochim.rds       ASVs after chimera removal
#   2_dada2/seqtab.final.rds        307-319 bp nonchimeric ASVs
#   2_dada2/ASVs.nochim.fasta       Sequences in seqtab.final.rds
#   2_dada2/ASVs.length_pass.fasta  Candidate ASV FASTA
#   2_dada2/ASV_ID_to_sequence.tsv  ASV IDs and sequences
#   2_dada2/sample_by_ASV.tsv      Sample-by-ASV abundance table
#   2_dada2/ASV_by_sample.tsv      ASV-by-sample abundance table
#   2_dada2/filtered/              Quality-filtered paired reads
#   2_dada2/plots/                 Quality and error-model plots
#   2_dada2/logs/                  Read tracking and length QC
#
# SETTINGS
#   WORK_DIR: /storage2/jihoonkim/Arielina_eDNA/analyses/1_20261103
#   DADA2: truncLen 240/200, maxEE 2/2, pseudo-pooling, LOESS errors.
#   Merging: minOverlap 100 bp, maxMismatch 0, trimOverhang FALSE.
#   Chimeras: consensus removal before final 307-319 bp length selection.
#   ASV IDs: descending total abundance, then lexicographic sequence order.
#   COI_WORKDIR is passed by the SLURM wrapper; otherwise use script directory.
# ==============================================================================

args <- commandArgs(trailingOnly=FALSE)
script_arg <- grep("^--file=", args, value=TRUE)
WORKDIR <- Sys.getenv("COI_WORKDIR", unset="")

if (!nzchar(WORKDIR)) {
  if (length(script_arg) != 1L) {
    stop("Set COI_WORKDIR or run via Rscript")
  }
  WORKDIR <- dirname(normalizePath(sub("^--file=", "", script_arg)))
}

WORKDIR <- normalizePath(WORKDIR, mustWork=TRUE)
IN_DIR <- file.path(WORKDIR, "1_trimming", "3_primer_trim")
OUT_DIR <- file.path(WORKDIR, "2_dada2")
PLOTS <- file.path(OUT_DIR, "plots")
LOGS <- file.path(OUT_DIR, "logs")
FILT <- file.path(OUT_DIR, "filtered")

THREADS <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "12"))
if (is.na(THREADS) || THREADS < 1L) THREADS <- 12L

# DADA2 filtering, error learning, merging and final sequence selection.
TRUNC_LEN <- c(240L, 200L)
MAX_EE <- c(2, 2)
TRUNC_Q <- 2L
MAX_N <- 0L
RM_PHIX <- TRUE
POOL_MODE <- "pseudo"
MIN_OVERLAP <- 100L
MAX_MISMATCH <- 0L
NBases <- 2e8
LEN_MIN <- 307L
LEN_MAX <- 319L

write_tsv <- function(df, path) {
  write.table(
    df, file=path, sep="\t",
    row.names=FALSE, quote=FALSE
  )
}

getN <- function(x) sum(getUniques(x))
mergedN <- function(x) {
  if (nrow(x) > 0L) sum(x$abundance) else 0
}

# ==============================================================================
# PREFLIGHT
# ==============================================================================

if (!dir.exists(IN_DIR)) {
  stop("Trimmed FASTQ directory missing: ", IN_DIR)
}

fnFs <- sort(list.files(
  IN_DIR, pattern="_R1_001\\.fastq$", full.names=TRUE
))
fnRs <- sort(list.files(
  IN_DIR, pattern="_R2_001\\.fastq$", full.names=TRUE
))

if (!length(fnFs) || length(fnFs) != length(fnRs)) {
  stop("No paired R1/R2 FASTQs, or unequal number of files: ", IN_DIR)
}

samples <- sub("_R1_001\\.fastq$", "", basename(fnFs))
expectedRs <- file.path(
  IN_DIR, paste0(samples, "_R2_001.fastq")
)

if (!identical(fnRs, sort(expectedRs))) {
  stop("R1/R2 sample pairing mismatch")
}

if (anyDuplicated(samples) ||
    any(!grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", samples))) {
  stop("Invalid or duplicate sample identifiers")
}

if (anyNA(file.size(c(fnFs, fnRs))) ||
    any(file.size(c(fnFs, fnRs)) == 0)) {
  stop("Empty FASTQ input")
}

if (file.exists(OUT_DIR)) {
  stop("Existing output is protected: ", OUT_DIR)
}

set.seed(1L)

dir.create(FILT, recursive=TRUE, showWarnings=FALSE)
dir.create(PLOTS, recursive=TRUE, showWarnings=FALSE)
dir.create(LOGS, recursive=TRUE, showWarnings=FALSE)

writeLines(c(
  paste("Start:", Sys.time()),
  paste("Working directory:", WORKDIR),
  paste("R:", R.version.string),
  paste("DADA2:", as.character(packageVersion("dada2"))),
  paste("Samples:", paste(samples, collapse=",")),
  paste("Threads:", THREADS),
  "Error model: loessErrfun (standard LOESS)",
  paste("nbases:", NBases),
  "randomize: FALSE",
  "MAX_CONSIST: 10",
  "truncLen: 240,200",
  "maxEE: 2,2",
  "truncQ: 2",
  "maxN: 0",
  "rm.phix: TRUE",
  "pool: pseudo",
  "minOverlap: 100",
  "maxMismatch: 0",
  "trimOverhang: FALSE",
  "Chimera: consensus (all lengths)",
  "Final length range after chimera: 307-319 inclusive",
  "Seed: 1"
), file.path(LOGS, "parameters.txt"))

cat("[INFO] Input samples:", length(samples), "\n")

# ==============================================================================
# STEP 1. GENERATE PAIRED-READ QUALITY PROFILES
# ==============================================================================

cat("[STEP 1] Generate paired-read quality profiles\n")

nplot <- min(6L, length(fnFs))

pdf(
  file.path(PLOTS, "01_qualityProfile_R1.pdf"),
  width=10, height=7
)
print(plotQualityProfile(fnFs[seq_len(nplot)]))
invisible(dev.off())

pdf(
  file.path(PLOTS, "02_qualityProfile_R2.pdf"),
  width=10, height=7
)
print(plotQualityProfile(fnRs[seq_len(nplot)]))
invisible(dev.off())

# ==============================================================================
# STEP 2. FILTER PAIRED READS BY LENGTH, QUALITY, Ns AND PHIX
# ==============================================================================

filtFs <- file.path(
  FILT, paste0(samples, "_F_filt.fastq.gz")
)
filtRs <- file.path(
  FILT, paste0(samples, "_R_filt.fastq.gz")
)

cat("[STEP 2] Filter trimmed read pairs (filterAndTrim)\n")

filter_result <- filterAndTrim(
  fnFs, filtFs, fnRs, filtRs,
  truncLen=TRUNC_LEN,
  maxN=MAX_N,
  maxEE=MAX_EE,
  truncQ=TRUNC_Q,
  rm.phix=RM_PHIX,
  compress=TRUE,
  multithread=THREADS,
  verbose=TRUE
)

if (nrow(filter_result) != length(samples) ||
    any(filter_result[, "reads.out"] <= 0)) {
  stop("Filtering failed or zero filtered reads")
}

write.csv(
  filter_result,
  file.path(LOGS, "filterAndTrim_readCounts.csv"),
  quote=FALSE
)

# ==============================================================================
# STEP 3. LEARN INDEPENDENT R1 AND R2 ERROR MODELS
# ==============================================================================

cat("[STEP 3] Learn R1/R2 errors using standard LOESS\n")

errF <- learnErrors(
  filtFs, nbases=NBases,
  errorEstimationFunction=loessErrfun,
  randomize=FALSE,
  MAX_CONSIST=10,
  multithread=THREADS,
  verbose=TRUE
)

errR <- learnErrors(
  filtRs, nbases=NBases,
  errorEstimationFunction=loessErrfun,
  randomize=FALSE,
  MAX_CONSIST=10,
  multithread=THREADS,
  verbose=TRUE
)

saveRDS(errF, file.path(OUT_DIR, "errF.rds"))
saveRDS(errR, file.path(OUT_DIR, "errR.rds"))

for (direction in c("R1", "R2")) {
  err <- if (direction == "R1") errF else errR

  pdf(
    file.path(PLOTS, paste0("03_errorModel_", direction, ".pdf")),
    width=10, height=7
  )
  print(plotErrors(err, nominalQ=TRUE))
  invisible(dev.off())
}

# ==============================================================================
# STEP 4. INFER ASVs USING PSEUDO-POOLING AND MERGE READ PAIRS
# ==============================================================================

cat("[STEP 4] Pseudo-pooled denoising\n")

dadaFs <- dada(
  filtFs, err=errF,
  pool=POOL_MODE, multithread=THREADS
)
cat("\n")

dadaRs <- dada(
  filtRs, err=errR,
  pool=POOL_MODE, multithread=THREADS
)
cat("\n")

names(dadaFs) <- samples
names(dadaRs) <- samples

saveRDS(dadaFs, file.path(OUT_DIR, "dadaFs.rds"))
saveRDS(dadaRs, file.path(OUT_DIR, "dadaRs.rds"))

cat("[STEP 4] mergePairs\n")

mergers <- mergePairs(
  dadaFs, filtFs,
  dadaRs, filtRs,
  minOverlap=MIN_OVERLAP,
  maxMismatch=MAX_MISMATCH,
  trimOverhang=FALSE,
  verbose=TRUE
)

names(mergers) <- samples
saveRDS(mergers, file.path(OUT_DIR, "mergers.rds"))

# ==============================================================================
# STEP 5. CONSTRUCT ASV TABLE AND REMOVE CHIMERAS
# ==============================================================================

cat("[STEP 5] Pre-chimera sequence table\n")

seqtab <- makeSequenceTable(mergers)
rownames(seqtab) <- samples

if (ncol(seqtab) == 0L) {
  stop("No merged ASVs")
}

saveRDS(
  seqtab,
  file.path(OUT_DIR, "seqtab.prechim.rds")
)

cat("[STEP 5] Consensus chimera removal (all lengths)\n")

seqtab.nochim <- removeBimeraDenovo(
  seqtab,
  method="consensus",
  multithread=THREADS,
  verbose=TRUE
)

if (ncol(seqtab.nochim) == 0L) {
  stop("No nonchimeric ASVs")
}

saveRDS(
  seqtab.nochim,
  file.path(OUT_DIR, "seqtab.nochim.rds")
)

# ==============================================================================
# STEP 6. APPLY ASV LENGTH FILTER AND EXPORT SEQUENCES/TABLES
# ==============================================================================

cat("[STEP 6] Length filter: 307-319 bp\n")

keep <- nchar(colnames(seqtab.nochim)) >= LEN_MIN &
        nchar(colnames(seqtab.nochim)) <= LEN_MAX

seqtab.final <- seqtab.nochim[, keep, drop=FALSE]

if (ncol(seqtab.final) == 0L) {
  stop("No ASVs in 307-319 bp range")
}

# Assign ASV identifiers in descending abundance and sequence order.
ord <- order(
  -colSums(seqtab.final), colnames(seqtab.final)
)
seqtab.final <- seqtab.final[, ord, drop=FALSE]

saveRDS(
  seqtab.final,
  file.path(OUT_DIR, "seqtab.final.rds")
)

seqs <- colnames(seqtab.final)
ids <- paste0("ASV_", seq_along(seqs))
totals <- colSums(seqtab.final)

writeLines(
  as.vector(rbind(paste0(">", ids), seqs)),
  file.path(OUT_DIR, "ASVs.nochim.fasta")
)

stopifnot(file.copy(
  file.path(OUT_DIR, "ASVs.nochim.fasta"),
  file.path(OUT_DIR, "ASVs.length_pass.fasta"),
  overwrite=FALSE
))

write_tsv(
  data.frame(
    ASV_ID=ids,
    sequence=seqs,
    length_bp=nchar(seqs),
    total_reads=as.numeric(totals)
  ),
  file.path(OUT_DIR, "ASV_ID_to_sequence.tsv")
)

sample_table <- data.frame(
  Sample=samples, seqtab.final,
  check.names=FALSE
)
colnames(sample_table) <- c("Sample", ids)

write_tsv(
  sample_table,
  file.path(OUT_DIR, "sample_by_ASV.tsv")
)

asv_table <- data.frame(
  ASV_ID=ids,
  t(seqtab.final),
  Total_reads=as.numeric(totals),
  check.names=FALSE
)
colnames(asv_table) <- c("ASV_ID", samples, "Total_reads")

write_tsv(
  asv_table,
  file.path(OUT_DIR, "ASV_by_sample.tsv")
)

# ==============================================================================
# STEP 7. SUMMARIZE READ RETENTION AND ASV LENGTH DISTRIBUTIONS
# ==============================================================================

cat("[STEP 7] Read tracking and length summaries\n")

track <- data.frame(
  sample=samples,
  input=as.numeric(filter_result[, "reads.in"]),
  filtered=as.numeric(filter_result[, "reads.out"]),
  denoisedF=as.numeric(vapply(dadaFs, getN, numeric(1))),
  denoisedR=as.numeric(vapply(dadaRs, getN, numeric(1))),
  merged=as.numeric(vapply(mergers, mergedN, numeric(1))),
  nonchim=as.numeric(rowSums(seqtab.nochim)),
  lenfilt=as.numeric(rowSums(seqtab.final))
)

write_tsv(
  track,
  file.path(LOGS, "read_tracking.tsv")
)

write.csv(
  track,
  file.path(LOGS, "track_reads.csv"),
  row.names=FALSE, quote=FALSE
)

stages <- c(
  "input", "filtered", "merged", "nonchim", "lenfilt"
)
reads <- vapply(
  stages,
  function(s) sum(track[[s]]),
  numeric(1)
)

write.csv(
  data.frame(
    stage=stages,
    reads=reads,
    percent_of_input=100*reads/reads[1]
  ),
  file.path(LOGS, "global_read_retention.csv"),
  row.names=FALSE, quote=FALSE
)

length_rows <- list()

for (stage in c("prechim", "nochim", "final")) {
  tab <- switch(
    stage,
    prechim=seqtab,
    nochim=seqtab.nochim,
    final=seqtab.final
  )

  L <- nchar(colnames(tab))
  abund <- colSums(tab)

  for (n in sort(unique(L))) {
    idx <- L == n

    length_rows[[length(length_rows)+1L]] <-
      data.frame(
        stage=stage,
        length_bp=n,
        ASVs=sum(idx),
        reads=sum(abund[idx])
      )
  }
}

write_tsv(
  do.call(rbind, length_rows),
  file.path(LOGS, "ASV_length_distribution_all.tsv")
)

write_tsv(
  data.frame(
    sample=samples,
    n_ASVs=rowSums(seqtab.final > 0),
    reads_final=rowSums(seqtab.final)
  ),
  file.path(LOGS, "asv_counts_per_sample.tsv")
)

# Validate read counts, sequence identifiers, and filtered lengths.
stopifnot(
  sum(seqtab.final) == sum(totals),
  ncol(seqtab.final) == length(ids),
  !anyDuplicated(seqs),
  all(grepl("^[ACGT]+$", seqs)),
  all(nchar(seqs) >= LEN_MIN & nchar(seqs) <= LEN_MAX),
  all(track$input >= track$filtered),
  all(track$filtered >= track$merged),
  all(track$merged >= track$nonchim),
  all(track$nonchim >= track$lenfilt)
)

capture.output(
  sessionInfo(),
  file=file.path(LOGS, "sessionInfo.txt")
)

writeLines(
  c(
    paste("End:", Sys.time()),
    paste("Samples:", length(samples)),
    paste("ASVs prechim:", ncol(seqtab)),
    paste("ASVs postchim:", ncol(seqtab.nochim)),
    paste("ASVs final:", ncol(seqtab.final)),
    paste("Reads final:", sum(seqtab.final))
  ),
  file.path(LOGS, "run_summary.txt")
)

cat("[QC PASS] Final ASV tables and read tracking validated\n")
cat(
  "[SUMMARY] ASVs prechim:", ncol(seqtab),
  "postchim:", ncol(seqtab.nochim),
  "length pass:", ncol(seqtab.final), "\n"
)
cat("[SUMMARY] Final reads:", sum(seqtab.final), "\n")

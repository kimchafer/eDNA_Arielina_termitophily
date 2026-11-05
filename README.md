[README.md](https://github.com/user-attachments/files/34424405/README.md)
# Body-surface eDNA metabarcoding pipeline for *Arielina peringueyi*

This repository provides the analysis scripts and a step-by-step tutorial for processing COI amplicons recovered from material adhering to the body surface of *Arielina peringueyi*. The workflow was developed to screen for potentially associated organisms, including candidate termite hosts. Sequence-based taxonomic assignments identify candidates for further study; they do **not**, on their own, demonstrate host identity, ecological association, or direct physical contact.

The tutorial follows the original four-library Illumina MiSeq paired-end dataset, from demultiplexed raw reads through ASV inference, MOTU curation, and taxonomic assignment. **The commands shown under “Core analysis” illustrate the relevant parts of the published scripts; run the complete scripts rather than running these excerpts separately.**

## Workflow

```text
NCBI SRA / demultiplexed paired-end FASTQ
                 |
                 v
1. Cutadapt          Adapter and primer trimming
                 |
                 v
2. DADA2             Quality filtering, error learning, denoising,
                     read merging, chimera removal, ASV selection
                 |
                 v
3. Subsequence       Exact ASV containment and conservative resolution
                 |
                 v
4. Swarm             ASV clustering into MOTUs (default d = 2)
                 |
                 v
5. VSEARCH + LULU    Sequence matching and abundance-based MOTU curation
                 |
                 v
6. MMseqs2           Taxonomic assignment at 80% and 97% minimum identity
                 |
                 v
7. Python            MOTU-level taxonomy comparison (TSV + XLSX)
```

## Contents

- [Repository structure](#repository-structure)
- [Software and computing requirements](#software-and-computing-requirements)
- [Raw data acquisition](#raw-data-acquisition)
- [Preparing the analysis workspace](#preparing-the-analysis-workspace)
- [Step 1. Adapter and primer trimming](#step-1-adapter-and-primer-trimming)
- [Step 2. ASV inference with DADA2](#step-2-asv-inference-with-dada2)
- [Step 3. ASV subsequence resolution](#step-3-asv-subsequence-resolution)
- [Step 4. MOTU clustering with Swarm](#step-4-motu-clustering-with-swarm)
- [Step 5. MOTU curation with VSEARCH and LULU](#step-5-motu-curation-with-vsearch-and-lulu)
- [Step 6. Taxonomic assignment with MMseqs2](#step-6-taxonomic-assignment-with-mmseqs2)
- [Step 7. Comparing 80% and 97% assignments](#step-7-comparing-80-and-97-assignments)
- [Example results](#example-results)
- [Interpretation and reproducibility notes](#interpretation-and-reproducibility-notes)

## Repository structure

```text
eDNA_Arielina_termitophily/
├── README.md
├── 1_sequence_processing/
│   ├── 1_trimming.sh
│   ├── 2_dada2_run.sh
│   ├── 2_dada2.r
│   └── 3_subsequence.sh
├── 2_motu_curation/
│   ├── 4_swarm.sh
│   └── 5_lulu.sh
└── 3_taxonomic_assignment/
    ├── 6_mmseqs_lca4.sh
    └── 7_taxonomy_comparison.py
```

The repository contains **source code**, not the full sequencing dataset or a prebuilt MMseqs2 reference database. Intermediate and final analysis files are generated in the configured working directory.

## Software and computing requirements

The supplied shell scripts were developed for a **Linux HPC cluster running Slurm**. The principal software dependencies are:

| Step | Main software | Role |
| --- | --- | --- |
| Raw data | [NCBI SRA Toolkit](https://github.com/ncbi/sra-tools) | Retrieve and convert archived reads |
| 1 | [Cutadapt](https://cutadapt.readthedocs.io/) | Trim adapters and degenerate COI primers |
| 2 | [DADA2](https://benjjneb.github.io/dada2/) in R | Infer ASVs and remove chimeras |
| 3 | Python 3 standard library | Identify and resolve exact subsequence containment |
| 4 | [Swarm](https://github.com/torognes/swarm), Python 3 | Cluster ASVs into MOTUs |
| 5 | [VSEARCH](https://github.com/torognes/vsearch), R [LULU](https://github.com/tobiasgf/lulu) | Identify similar MOTUs and perform curation |
| 6 | [MMseqs2](https://github.com/soedinglab/MMseqs2) | Assign taxonomy using a nucleotide taxonomy database |
| 7 | Python 3 standard library | Join and export taxonomy comparison tables |

The scripts assume that their required environments have already been installed. For a different cluster or computer, update the hard-coded Conda activation paths, `WORK_DIR`, `RAW_DIR`, MMseqs2 reference database, and `#SBATCH` partition/resource directives. A nucleotide MMseqs2 database with associated taxonomy information is required for Step 6; it is **not** downloaded by this pipeline. The dataset's specific software settings and available version information are recorded in each stage's logs.

**Important:** Steps 1–6 are submitted to Slurm; Step 7 is a short Python post-processing task and can be run directly. The published scripts also protect existing output directories from being overwritten. Complete each step and inspect its QC output before submitting the next one.

## Raw data acquisition

This tutorial starts from paired-end FASTQ files corresponding to the original dataset. Users can retrieve the sequencing reads from NCBI SRA under BioProject **[PRJNA1554401](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1554401)**. Open the BioProject or its [SRA Run Selector](https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA1554401), identify the appropriate **SRR run accessions**, and record which accession corresponds to each sample (`KPC0541`–`KPC0544`).

Install the [SRA Toolkit](https://github.com/ncbi/sra-tools) and, **for each SRR accession**, run:

```bash
prefetch <SRR_accession>
fasterq-dump <SRR_accession> --split-files -e 8
gzip <SRR_accession>_1.fastq <SRR_accession>_2.fastq
```

Replace `<SRR_accession>` with an actual run accession (without angle brackets). `prefetch` downloads the SRA run; `fasterq-dump --split-files` converts the paired reads to separate FASTQ files; `-e 8` uses eight processing threads. Since `fasterq-dump` produces uncompressed FASTQ, `gzip` is run separately. The SRA extraction step can require substantial temporary disk space.

The resulting files are ordinarily named `SRRxxxxxx_1.fastq.gz` and `SRRxxxxxx_2.fastq.gz`. **Rename them according to the verified SRR-to-sample mapping** and place the files in the configured `raw_data/` directory before running Step 1. For example, *only when the selected SRR is confirmed to represent KPC0541*:

```bash
mv SRRxxxxxx_1.fastq.gz raw_data/KPC0541_E_Lib_S99_L001_R1_001.fastq.gz
mv SRRxxxxxx_2.fastq.gz raw_data/KPC0541_E_Lib_S99_L001_R2_001.fastq.gz
```

The original input naming pattern was:

```text
raw_data/
├── KPC0541_E_Lib_S99_L001_R1_001.fastq.gz
├── KPC0541_E_Lib_S99_L001_R2_001.fastq.gz
├── KPC0542_E_Lib_S100_L001_R1_001.fastq.gz
├── KPC0542_E_Lib_S100_L001_R2_001.fastq.gz
├── KPC0543_E_Lib_S101_L001_R1_001.fastq.gz
├── KPC0543_E_Lib_S101_L001_R2_001.fastq.gz
├── KPC0544_E_Lib_S102_L001_R1_001.fastq.gz
└── KPC0544_E_Lib_S102_L001_R2_001.fastq.gz
```

The trimming script automatically discovers files ending in `_R1_001.fastq.gz` and `_R2_001.fastq.gz` and recognizes the original `_E_Lib_S..._L...` pattern to derive sample IDs such as `KPC0541`. **Do not assign SRR files to these sample names without checking the SRA metadata.** Users applying the tutorial to different datasets may use different sample prefixes, but should review the assumptions about sample naming and the four fixed sample columns in Step 7.

## Preparing the analysis workspace

Clone the repository:

```bash
git clone https://github.com/kimchafer/eDNA_Arielina_termitophily.git
cd eDNA_Arielina_termitophily
```

The workflow uses a common working directory for scripts and generated outputs. As an example, copy the source scripts into a working directory:

```bash
REPO_DIR="$PWD"
WORK_DIR=/path/to/your/analysis
RAW_DIR=/path/to/your/raw_data
mkdir -p "$WORK_DIR" "$RAW_DIR"

cp "$REPO_DIR/1_sequence_processing/"{1_trimming.sh,2_dada2_run.sh,2_dada2.r,3_subsequence.sh} "$WORK_DIR/"
cp "$REPO_DIR/2_motu_curation/"{4_swarm.sh,5_lulu.sh} "$WORK_DIR/"
cp "$REPO_DIR/3_taxonomic_assignment/"{6_mmseqs_lca4.sh,7_taxonomy_comparison.py} "$WORK_DIR/"
```

**Edit the paths inside the copied scripts before executing them.** Defining `WORK_DIR` in the terminal does **not** automatically override hard-coded `WORK_DIR` values inside each shell script. Step 2's Slurm wrapper expects `2_dada2.r` in the working directory. The original MMseqs2 reference path refers to an institution-specific database and must be replaced with a suitable local database.

The examples below assume these paths have been configured and the raw data are present. Commands shown in **Core analysis** blocks are **excerpts for explanation** and are not a substitute for the executable scripts.

---

## Step 1. Adapter and primer trimming

**Source:** [`1_trimming.sh`](1_sequence_processing/1_trimming.sh)  
**Input:** demultiplexed paired `*.fastq.gz` reads  
**Output:** primer-trimmed paired FASTQ files and Cutadapt/QC reports

The target is the approximately 313-bp COI amplicon amplified with the degenerate primers:

| Primer | Sequence (5′–3′) |
| --- | --- |
| mlCOIintF | `GGWACWGGWTGAACWGTWTAYCCYCC` |
| jgHCO2198 | `TAIACYTCIGGRTGICCRAARAAYCA` |

The reverse primer contains inosine (`I`). The script retains the original reverse-primer specification for 5′ matching and replaces `I` with `N` when building the reverse-complement 3′ search sequence.

### Core analysis

Cutadapt is run in **three successive passes**, rather than in a single combined command. The excerpts below retain the principal options from the published script:

**Pass 1 — sequencing adapters**

```bash
cutadapt -j "$THREADS" \
    -a "$ADAPT_R1" -A "$ADAPT_R2" \
    -O 15 -e 0.20 --no-indels \
    --trim-n --minimum-length 50 --pair-filter=any \
    -o "$adapter_r1" -p "$adapter_r2" \
    "$raw_r1" "$raw_r2"
```

**Pass 2 — 5′ COI primers**

```bash
cutadapt -j "$THREADS" \
    -g "X${FWD};e=0.20" -G "X${REV_I};e=0.20" \
    -O 17 --no-indels --match-read-wildcards \
    --discard-untrimmed --pair-filter=both --trim-n \
    -o "$fivep_r1" -p "$fivep_r2" \
    "$adapter_r1" "$adapter_r2"
```

**Pass 3 — 3′ COI primers and final filters**

```bash
cutadapt -j "$THREADS" \
    -a "$REV_RC_N" -A "$RC_FWD" \
    -O 15 -e 0.20 --no-indels --match-read-wildcards \
    --trim-n --max-n 0 --minimum-length 50 --pair-filter=any \
    -o "$final_r1" -p "$final_r2" \
    "$fivep_r1" "$fivep_r2"
```

### What the principal parameters do

| Option | Purpose and rationale |
| --- | --- |
| `-a` / `-A`, `-g` / `-G` | Specify adapter or primer sequences for R1 and R2, respectively. The separate passes make the trimming logic explicit. |
| `-O 15` or `-O 17` | Require at least 15 aligned bases for adapter/3′ primer matches, or 17 for 5′ primers. A minimum overlap helps avoid short, accidental matches while allowing incomplete primer sequences. |
| `-e 0.20` | Permit a maximum error rate of 20% within an adapter/primer match, accommodating mismatches in divergent or imperfect reads. |
| `--no-indels` | Prevent insertions/deletions within adapter/primer matching. This is a *matching* rule, not a claim that COI amplicons contain no indels. |
| `X` prefix and `--match-read-wildcards` | Use non-internal 5′ matching (including partial matches at the read end) and allow ambiguous nucleotide codes in matching. |
| `--discard-untrimmed --pair-filter=both` | During the **5′ primer pass only**, discard a read pair if neither read has a recognized 5′ primer. One matched primer is sufficient to retain the pair. |
| `--pair-filter=any` | In passes with length/ambiguous-base filters, discard a pair if **either** member fails the specified filter. |
| `--minimum-length 50`, `--max-n 0` | Require both final reads to be at least 50 nt long and to contain no ambiguous `N` bases. |
| `--trim-n` | Remove terminal `N` bases before subsequent filtering. |

The 3′ primer match is **optional**: absence of a 3′ primer does not, by itself, cause rejection. The script also validates FASTQ integrity, paired-read synchronization, and length distributions. These diagnostic checks do not introduce an additional, separate read-removal step.

### Run the stage

```bash
cd "$WORK_DIR"
sbatch 1_trimming.sh
```

**Main files to inspect:**

```text
1_trimming/3_primer_trim/<sample>_R1_001.fastq
1_trimming/3_primer_trim/<sample>_R2_001.fastq
1_trimming/logs/read_counts_summary.tsv
1_trimming/logs/<sample>_adapter.txt
1_trimming/logs/<sample>_5prime.txt
1_trimming/logs/<sample>_3prime.txt
```

```bash
head -n 6 1_trimming/logs/read_counts_summary.tsv
```

**QC checkpoint:** Verify sample pairing and the fraction of read pairs retained after each trimming pass. The final uncompressed paired FASTQ files are used by DADA2.

---

## Step 2. ASV inference with DADA2

**Sources:** [`2_dada2_run.sh`](1_sequence_processing/2_dada2_run.sh) and [`2_dada2.r`](1_sequence_processing/2_dada2.r)  
**Input:** paired FASTQ from Step 1  
**Output:** denoised, merged, nonchimeric, length-selected ASVs and abundance tables

The Slurm wrapper launches the R script. DADA2 models sequence-specific errors rather than treating all unique read sequences as distinct biological variants. Its major operations are quality filtering, error estimation, denoising, read merging, chimera removal, and final length selection.

### Core analysis

**Quality filtering:**

```r
TRUNC_LEN <- c(240L, 200L)
MAX_EE <- c(2, 2)
TRUNC_Q <- 2L
MAX_N <- 0L
RM_PHIX <- TRUE

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
```

**Error learning and denoising:**

```r
errF <- learnErrors(
  filtFs, nbases=2e8,
  errorEstimationFunction=loessErrfun,
  randomize=FALSE, MAX_CONSIST=10,
  multithread=THREADS, verbose=TRUE
)

errR <- learnErrors(
  filtRs, nbases=2e8,
  errorEstimationFunction=loessErrfun,
  randomize=FALSE, MAX_CONSIST=10,
  multithread=THREADS, verbose=TRUE
)

dadaFs <- dada(filtFs, err=errF, pool="pseudo", multithread=THREADS)
dadaRs <- dada(filtRs, err=errR, pool="pseudo", multithread=THREADS)
```

**Merge, remove chimeras, and select amplicon lengths:**

```r
mergers <- mergePairs(
  dadaFs, filtFs, dadaRs, filtRs,
  minOverlap=100L, maxMismatch=0L,
  trimOverhang=FALSE, verbose=TRUE
)

seqtab <- makeSequenceTable(mergers)
seqtab.nochim <- removeBimeraDenovo(
  seqtab, method="consensus", multithread=THREADS, verbose=TRUE
)

keep <- nchar(colnames(seqtab.nochim)) >= 307L &
        nchar(colnames(seqtab.nochim)) <= 319L
seqtab.final <- seqtab.nochim[, keep, drop=FALSE]
```

### What the principal parameters do

| Setting | Purpose and rationale |
| --- | --- |
| `truncLen=c(240,200)` | Truncate R1 and R2 to 240 and 200 nt. This preserves enough overlapping sequence for a ~313-bp fragment while limiting lower-quality read tails. |
| `maxEE=c(2,2)` | Require no more than two expected errors in each retained read, a quality-aware criterion based on per-base Phred scores. |
| `truncQ=2` | Truncate at the first base whose quality is at or below the specified threshold; also subject to length filtering. |
| `maxN=0` | Reject ambiguous bases in filtered reads. |
| `rm.phix=TRUE` | Remove reads recognized as PhiX control sequence. |
| `learnErrors(..., nbases=2e8, loessErrfun)` | Learn separate substitution-error models for forward and reverse reads using up to the requested number of bases and the standard LOESS fit. |
| `pool="pseudo"` | Share evidence of variants across samples approximately, improving sensitivity to low-abundance variants supported elsewhere without full pooling. |
| `minOverlap=100`, `maxMismatch=0` | Require a substantial, mismatch-free overlap when merging paired reads. |
| `method="consensus"` | Remove likely chimeras based on agreement across samples rather than relying on a single pooled sample. |
| `307–319 bp` | Keep amplicons within 313 ± 6 bp after chimera removal. This allows limited length variation, including potential short COI indels, without accepting arbitrary amplicon lengths. |

The pipeline does **not** apply an additional, explicit singleton-read filter after DADA2. ASV inference itself depends on the learned error model and the selected pseudo-pooling strategy.

### Run the stage

```bash
cd "$WORK_DIR"
sbatch 2_dada2_run.sh
```

**Main files to inspect:**

```text
2_dada2/seqtab.final.rds
2_dada2/ASVs.nochim.fasta
2_dada2/ASV_ID_to_sequence.tsv
2_dada2/sample_by_ASV.tsv
2_dada2/ASV_by_sample.tsv
2_dada2/plots/01_qualityProfile_R1.pdf
2_dada2/plots/02_qualityProfile_R2.pdf
2_dada2/logs/read_tracking.tsv
```

```bash
head -n 6 2_dada2/logs/read_tracking.tsv
head -n 6 2_dada2/ASV_ID_to_sequence.tsv
```

**QC checkpoint:** Examine read retention at filtering, merging, chimera removal, and length selection. Inspect the quality and error-model plots. Despite its filename, `ASVs.nochim.fasta` contains the **final length-selected** nonchimeric ASVs.

---

## Step 3. ASV subsequence resolution

**Source:** [`3_subsequence.sh`](1_sequence_processing/3_subsequence.sh)  
**Input:** ASV FASTA, ASV sequence mapping, and sample-by-ASV table from Step 2  
**Output:** resolved ASVs and explicit original-to-final ASV mapping

Some ASVs can be shorter fragments fully contained within longer ASVs. This custom Python stage tests **all shorter–longer pairs** for 100% exact nucleotide containment and applies conservative, transparent criteria before transferring counts.

### Core decision logic

The following excerpt expresses the decision implemented by the script:

```python
if not parents:
    status, parent = 'retained_no_exact_parent', ''
elif (len(seq[child]) - 313) % 3 == 0:
    status, parent = 'manual_review_313_frame_class', ''
elif len(parents) > 1:
    status, parent = 'manual_review_multiple_parents', ''
elif n_internal or not n_terminal:
    status, parent = 'manual_review_internal_match', ''
else:
    status, parent = 'merged_unique_terminal', parents[0][0]
    merge_to[child] = parent
```

| Rule | Why it is applied |
| --- | --- |
| 100% identity, no gaps | Avoid merging sequences that differ at nucleotide positions. |
| Strictly longer parent | Only resolve a shorter, contained ASV into a longer one; equal-length variants are kept distinct. |
| Unique parent | Prevent arbitrary selection among multiple possible longer ASVs. |
| Terminal, not internal, containment | Restrict automatic merging to cases consistent with a missing sequence end rather than a wholly internal fragment. |
| `(child_length − 313) % 3 != 0` | Automatically resolve only the designated non-codon-length class; matched sequences in the 313-frame class are flagged for review rather than merged. This is a conservative **workflow rule**, not a demonstration that a given sequence is an error. |
| Sample-specific transfer | Preserve each sample's original read abundance and the identity of the contributing ASV. |

If an ASV is automatically merged, its per-sample counts are transferred to the ultimately retained parent. Ambiguous matches remain available in the manual-review report.

### Run the stage

```bash
cd "$WORK_DIR"
sbatch 3_subsequence.sh
```

**Main files to inspect:**

```text
3_subsequence/asvs.subsequence_matches.tsv
3_subsequence/asvs.subsequence_resolution.tsv
3_subsequence/asvs.subsequence_manual_review.tsv
3_subsequence/asvs.subsequence_resolved.fasta
3_subsequence/asv_id_to_sequence_resolved.tsv
3_subsequence/sample_by_asv_resolved.tsv
3_subsequence/asvs.subsequence_summary.tsv
```

```bash
cat 3_subsequence/asvs.subsequence_summary.tsv
```

**QC checkpoint:** Check proposed automatic parent–child matches and any manual-review entries. The number of ASVs may decrease, but **global and per-sample read totals must not change**.

---

## Step 4. MOTU clustering with Swarm

**Source:** [`4_swarm.sh`](2_motu_curation/4_swarm.sh)  
**Input:** resolved ASV FASTA and abundances from Step 3  
**Output:** MOTU clusters, representative sequences, abundance tables, and ASV membership

Swarm groups related ASVs using a local sequence-difference threshold instead of assigning every sequence pair to an independent fixed-identity bin. The script creates abundance-annotated input labels of the form `>ASV_ID;size=N;` and then calls Swarm.

### Core analysis

```bash
# d=2: the default analysis
swarm -d "$D_VAL" -t "$THREADS" -z \
    -o "$CLUSTERS" -s "$STATS" -w "$SEEDS" "$SWARM_FASTA"

# d=1: the script additionally supplies -f (fastidious mode)
```

| Option | Purpose and rationale |
| --- | --- |
| `-d 2` | Set the local clustering distance to two nucleotide differences; this is the reference analysis setting. It is **not** a global 98% identity cutoff. |
| `-f` when `d=1` | Use Swarm's fastidious clustering mode for the one-difference analysis. It is not used for the default `d=2` run. |
| `-z` | Read sequence abundance from the `;size=N;` FASTA annotations. |
| `-t` | Use the number of threads allocated by Slurm. |
| `-o`, `-s`, `-w` | Write cluster membership, cluster statistics, and representative/seed sequences, respectively. |

The post-processing code verifies that each input ASV belongs to exactly one MOTU, builds MOTU-by-sample counts, and records each MOTU's seed ASV. Swarm clustering **does not perform taxonomic assignment** and does not itself determine biological species boundaries.

### Run the stage

```bash
cd "$WORK_DIR"
sbatch 4_swarm.sh 2
```

The `2` is the clustering distance. Separate analyses can be run with `sbatch 4_swarm.sh 1` or `sbatch 4_swarm.sh 3`; the script writes to separate `ALL_d1/`, `ALL_d2/`, or `ALL_d3/` directories.

**Main files to inspect (`d=2`):**

```text
4_swarm_motu/ALL_d2/ALL_d2.swarm
4_swarm_motu/ALL_d2/ALL_d2.stats
4_swarm_motu/ALL_d2/ALL_d2_MOTU.fasta
4_swarm_motu/ALL_d2/ALL_d2_MOTU_table.tsv
4_swarm_motu/ALL_d2/ALL_d2_members.tsv
4_swarm_motu/ALL_d2/ALL_d2_cluster_summary.tsv
4_swarm_motu/ALL_d2/ALL_d2_summary.tsv
```

```bash
cat 4_swarm_motu/ALL_d2/ALL_d2_summary.tsv
```

**QC checkpoint:** Confirm cluster accounting, seed identity, and total/sample read conservation. The `singleton_ASV_MOTUs` statistic counts MOTUs containing **one ASV**, not necessarily MOTUs supported by only one read.

---

## Step 5. MOTU curation with VSEARCH and LULU

**Source:** [`5_lulu.sh`](2_motu_curation/5_lulu.sh)  
**Input:** Swarm MOTU representative FASTA and abundance table  
**Output:** curated MOTUs, parent mapping, and abundance/QC summaries

LULU identifies lower-abundance MOTUs that may represent erroneous variants of more abundant MOTUs. The script first identifies sufficiently similar sequence pairs with VSEARCH, then uses their relative abundance and sample co-occurrence patterns to decide whether a child MOTU should be merged into a parent.

### Core analysis

**VSEARCH candidate matching:**

```bash
vsearch \
    --usearch_global "$FASTA" --db "$FASTA" \
    --strand plus --self \
    --id 0.97 --iddef 1 --query_cov 0.90 \
    --maxaccepts 0 --maxrejects 0 \
    --threads "$THREADS" \
    --userout "$MATCHLIST" \
    --userfields query+target+id
```

**LULU curation in R:**

```r
result <- lulu::lulu(
  otu, matches,
  minimum_ratio_type="min",
  minimum_ratio=1,
  minimum_match=97,
  minimum_relative_cooccurence=0.95
)
```

| Setting | Purpose and rationale |
| --- | --- |
| `--usearch_global` | Globally align MOTU representatives to identify candidate related sequences. |
| `--strand plus` | Compare the forward orientation; the COI amplicons are expected in a consistent orientation after primer processing. |
| `--self` | Suppress self-comparisons when the query and reference are the same FASTA. |
| `--id 0.97`, `--iddef 1` | Request candidate matches with at least 97% identity under VSEARCH identity definition 1. Identity is computed over the alignment, not simply by dividing mismatches by 313. |
| `--query_cov 0.90` | Require at least 90% of the query sequence to be represented in the alignment. |
| `--maxaccepts 0`, `--maxrejects 0` | Avoid prematurely limiting the number of accepted/rejected candidate hits. |
| `minimum_match=97` | Apply the LULU similarity condition to the VSEARCH match list. |
| `minimum_ratio_type="min"`, `minimum_ratio=1` | Require the candidate parent to be at least as abundant as the child under LULU's minimum-ratio rule. |
| `minimum_relative_cooccurence=0.95` | Require strong relative co-occurrence of the child with its proposed parent across the available samples. |

High similarity and co-occurrence support **curation decisions**, not proof that two MOTUs are the same biological species. In particular, observing both candidates in only one sample can satisfy a high relative co-occurrence threshold without strong replication.

### Run the stage

```bash
cd "$WORK_DIR"
sbatch 5_lulu.sh 2
```

**Main files to inspect (`d=2`):**

```text
5_lulu/ALL_d2/ALL_d2_vsearch_matchlist.tsv
5_lulu/ALL_d2/LULU_curated_MOTU.fasta
5_lulu/ALL_d2/LULU_curated_MOTU_table.tsv
5_lulu/ALL_d2/LULU_MOTU_map.tsv
5_lulu/ALL_d2/LULU_sample_summary.tsv
5_lulu/ALL_d2/LULU_summary.tsv
```

```bash
cat 5_lulu/ALL_d2/LULU_summary.tsv
```

**QC checkpoint:** Review the parent/daughter mapping and per-sample read conservation. LULU reassigns counts to retained MOTUs; it does not remove those reads. Subsequent taxonomy uses the **LULU-curated representative FASTA**, not the pre-LULU Swarm FASTA.

---

## Step 6. Taxonomic assignment with MMseqs2

**Source:** [`6_mmseqs_lca4.sh`](3_taxonomic_assignment/6_mmseqs_lca4.sh)  
**Input:** LULU-curated MOTU representatives and a nucleotide taxonomy reference database  
**Output:** two MOTU-level taxonomic assignment TSVs and search/QC logs

This stage compares two minimum sequence identity thresholds, **80%** and **97%**, while keeping the remaining settings identical. The lower threshold provides a broader taxonomic placement for divergent or poorly represented sequences; the more stringent threshold assesses which assignments persist under more restrictive matching. These are **two separate analyses**, not a sequential filtering of 80% hits to 97%.

### Core analysis

The full script constructs an MMseqs2 query database and runs the following taxonomy command **twice**, with `identity=0.80` and `identity=0.97`:

```bash
mmseqs createdb "$QUERY_FASTA" "$QUERY_DB"

mmseqs taxonomy \
    "$QUERY_DB" "$REF_DB" "$result" "$tmp" \
    --tax-lineage 1 --threads "$THREADS" \
    --search-type 3 --orf-filter 0 \
    --lca-mode 4 --alignment-mode 3 --strand 2 \
    -s 5 -e 1e-5 \
    --min-seq-id "$identity" --min-aln-len 300 \
    -c 0.80 --cov-mode 2 \
    --max-accept 20 --max-rejected 50

mmseqs createtsv "$QUERY_DB" "$result" "${result}.tsv" \
    --threads "$THREADS"
```

Here, `$tmp` is shorthand for the run-specific temporary directory used by the actual script. The executable uses distinct result and temporary paths for the two identity thresholds.

| Option | Purpose and rationale |
| --- | --- |
| `--search-type 3` | Perform a nucleotide search against the nucleotide reference database. |
| `--orf-filter 0` | Disable the ORF-filter setting for these short nucleotide amplicons. |
| `--lca-mode 4` | Apply the selected MMseqs2 LCA classification mode consistently in both comparisons. |
| `--alignment-mode 3` | Use the specified alignment calculation mode for the search. |
| `--strand 2` | Search both sequence orientations. |
| `-s 5` | Specify search sensitivity; higher settings may increase sensitivity and computational cost. |
| `-e 1e-5` | Apply the E-value threshold to candidate alignments. |
| `--min-seq-id 0.80` / `0.97` | Compare permissive and stringent minimum identity conditions. These values are **acceptance thresholds, not the observed identity of each MOTU**. |
| `--min-aln-len 300` | Require an alignment spanning at least 300 nt of the approximately 313-nt COI target. |
| `-c 0.80 --cov-mode 2` | Apply a minimum 80% query-coverage criterion. |
| `--max-accept 20`, `--max-rejected 50` | Limit the number of accepted and rejected search candidates considered by the search procedure. |
| `--tax-lineage 1` | Include taxonomic lineage information in the output. |

An assignment that becomes `unclassified` at 97% should **not** be reinterpreted as evidence that the sequence is absent or belongs outside the originally assigned higher taxon. Taxonomy is reference-dependent; the two conditions assess the stability of the reported assignments under the configured criteria.

### Run the stage

```bash
cd "$WORK_DIR"
sbatch 6_mmseqs_lca4.sh
```

**Main files to inspect:**

```text
6_mmseqs_lca4/results/taxonomy_lca4_a20_80.tsv
6_mmseqs_lca4/results/taxonomy_lca4_a20_97.tsv
6_mmseqs_lca4/logs/input_metrics.tsv
6_mmseqs_lca4/logs/taxonomy_qc.tsv
```

```bash
head -n 5 6_mmseqs_lca4/results/taxonomy_lca4_a20_80.tsv
head -n 5 6_mmseqs_lca4/results/taxonomy_lca4_a20_97.tsv
cat 6_mmseqs_lca4/logs/taxonomy_qc.tsv
```

**QC checkpoint:** Check that both assignments report each of the input curated MOTU IDs exactly once, that output column counts are correct, and that the MMseqs2 reference database and parameter settings have been recorded. The taxonomy TSV does **not** by itself contain per-hit alignment identity, coverage, or a full list of reference hits.

---

## Step 7. Comparing 80% and 97% assignments

**Source:** [`7_taxonomy_comparison.py`](3_taxonomic_assignment/7_taxonomy_comparison.py)  
**Input:** both MMseqs2 taxonomy TSVs and the LULU-curated MOTU abundance table  
**Output:** joined 21-column TSV and plain XLSX worksheet, plus QC logs

The final Python script joins the two taxonomy conditions **by MOTU ID** (not by row position) and appends the original sample-by-MOTU read counts. The merged table retains taxonomic names, ranks, NCBI TaxIDs, order/family information derived from lineage, and whether each assignment was retained at the more stringent threshold.

### Core comparison logic

```python
def compare(t80, t97):
    a, b = t80[0], t97[0]  # NCBI TaxIDs
    if a == 0 and b == 0:
        return "Unclassified both"
    if a != 0 and b == 0:
        return "Lost at 97%"
    if a == 0 and b != 0:
        return "Gained at 97%"
    return "Same TaxID" if a == b else "Changed TaxID"
```

| Comparison value | Meaning |
| --- | --- |
| `Same TaxID` | Both conditions returned the same NCBI taxon. |
| `Lost at 97%` | The 80% analysis returned a taxon but the 97% analysis returned no assignment. |
| `Gained at 97%` | The 97% analysis returned a taxon but the 80% analysis did not. |
| `Changed TaxID` | Both analyses assigned taxa, but to different NCBI TaxIDs. |
| `Unclassified both` | Neither condition returned a taxonomic assignment. |

This stage does **not** reclassify sequences or change their read counts. The Excel file is deliberately a plain worksheet without custom colors, borders, filters, or frozen panes.

### Run the stage

```bash
cd "$WORK_DIR"
python3 7_taxonomy_comparison.py
```

**Main files to inspect:**

```text
7_taxonomy_comparison/MOTU_comparison.tsv
7_taxonomy_comparison/MOTU_comparison.xlsx
7_taxonomy_comparison/logs/comparison_summary.tsv
```

```bash
cat 7_taxonomy_comparison/logs/comparison_summary.tsv
```

**QC checkpoint:** All MOTU IDs in both taxonomy files must match those in the LULU abundance table; totals must remain unchanged. The script expects the example study's four sample columns `KPC0541`–`KPC0544`; adapt the header configuration if analyzing a different sample set.

---

## Example results

The following are **observed results for the four-library dataset**, not required values for another study:

| Stage | Retained units | Total reads |
| --- | ---: | ---: |
| DADA2 | 293 ASVs | 579,596 |
| Subsequence resolution | 292 ASVs | 579,596 |
| Swarm (`d=2`) | 268 MOTUs | 579,596 |
| LULU | 267 MOTUs | 579,596 |

At 80% minimum sequence identity, 40 of the 267 curated MOTUs were assigned within Insecta. At 97%, nine MOTUs retained an Insecta assignment. The comparison stage reported 48 unchanged TaxIDs, 163 assignments lost at 97%, and 56 unclassified MOTUs under both settings.

These numbers are presented as **tutorial checkpoints** only. The code does not require another dataset to reproduce them.

## Interpretation and reproducibility notes

- **Read conservation:** Steps 3–5 redistribute reads among ASVs or MOTUs rather than removing reads. A difference between input and output read totals indicates an accounting problem; it is not evidence of biological loss.
- **ASVs versus MOTUs:** ASVs are denoised sequence variants; MOTUs are sequence clusters curated by defined similarity and abundance/co-occurrence rules. Neither automatically represents a biological species.
- **Low-abundance taxonomic signals:** A few reads assigned to a particular taxon warrant careful examination of sample handling, negative controls, possible contamination or index misassignment, reference coverage, and independent ecological verification.
- **Reference database:** MMseqs2 classifications are dependent on the composition, taxonomy mapping, and version of the nucleotide database. Record the database provenance if you apply this workflow to a different reference.
- **Threshold interpretation:** `--min-seq-id` sets a minimum search identity; it is not a per-MOTU identity estimate. A lack of assignment at 97% is not proof of non-detection.
- **Coding-sequence QC:** No MACSE or other translation/frameshift-filtering stage is part of this published workflow.
- **Hardware and paths:** The supplied Slurm resources and absolute paths describe the original HPC environment and should be adapted locally. Do not create stage output directories in advance; the scripts have overwrite-protection safeguards.

The linked source scripts are the authoritative record of executable settings. The code excerpts in this tutorial are intentionally shorter and provided to explain how the principal methods and parameters are used.

# Body-surface eDNA metabarcoding pipeline for *Arielina peringueyi*
<img width="1286" height="807" alt="cover" src="https://github.com/user-attachments/assets/1d8ab01b-4cad-47e3-b36a-b99e4746e5b6" />

This repository provides the analysis scripts and a step-by-step tutorial for processing COI amplicons recovered from material adhering to the body surface of *Arielina peringueyi* beetles (Insecta: Coleoptera: Scarabaeidae: Cetoniinae: Cremastocheilini). As this species belongs to Cremastocheilini, a tribe known for social parasitism, the workflow was developed to screen for associated organisms, with particular emphasis on identifying potential termite hosts.

The tutorial follows the original four-library Illumina MiSeq paired-end dataset, from demultiplexed raw reads through ASV inference, MOTU curation, and taxonomic assignment. **The commands shown under “Core analysis” illustrate the relevant parts of the published scripts; run the complete scripts rather than running these excerpts separately.** ChatGPT (OpenAI) assisted with script development and troubleshooting, while all analytical decisions, parameter settings, and validation were performed by the authors.

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

The repository contains **source code**, not the full sequencing dataset or a prebuilt MMseqs2 reference database. Intermediate and final analysis files are generated in the working directory specified by the user when running the pipeline.

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

The trimming script expects paired FASTQ files ending in `_R1_001.fastq.gz` and `_R2_001.fastq.gz` and recognizes the original `_E_Lib_S..._L...` naming pattern to identify samples. Users working with differently named files should either modify the filename-parsing logic in `1_trimming.sh` or rename their files to match the expected pattern, ensuring that sample identities are correctly assigned. Users applying the pipeline to other datasets should also review sample-naming assumptions in downstream scripts, including the fixed sample columns in Step 7.

## Preparing the analysis workspace

Before running the pipeline, prepare separate directories for raw sequencing data and analysis outputs. Download the eight script files covering Steps 1–7 from this repository and place them in a common working directory, following the structure below.

```text
project/
├── raw_data/
│   ├── KPC0541_E_Lib_S99_L001_R1_001.fastq.gz
│   ├── KPC0541_E_Lib_S99_L001_R2_001.fastq.gz
│   └── ... (remaining paired-end FASTQ files)
│
└── analyses/
    └── run_01/                         # WORK_DIR
        ├── 1_trimming.sh
        ├── 2_dada2_run.sh
        ├── 2_dada2.r
        ├── 3_subsequence.sh
        ├── 4_swarm.sh
        ├── 5_lulu.sh
        ├── 6_mmseqs_lca4.sh
        └── 7_taxonomy_comparison.py
```

For example, create the directories using:

```bash
PROJECT_DIR="/path/to/your/project"
mkdir -p "$PROJECT_DIR/raw_data" "$PROJECT_DIR/analyses/run_01"
```

The scripts are available in the `1_sequence_processing/`, `2_motu_curation/`, and `3_taxonomic_assignment/` folders of this repository.

**Before execution, update the directory paths defined within the scripts** to match your local environment. In particular, verify the working directory (`WORK_DIR`), raw-data paths, and the MMseqs2 reference database path. Setting `WORK_DIR` in the terminal does not automatically override paths hard-coded inside individual scripts. The DADA2 wrapper (`2_dada2_run.sh`) also expects `2_dada2.r` in the working directory.

Intermediate and final results will be generated in step-specific subdirectories under the configured working directory.

The following tutorial assumes that all input files and paths have been configured. Commands shown under **Core analysis** are excerpts for explanation; execute the complete scripts to reproduce the analyses.

---

## Step 1. Adapter and primer trimming

**Script:** [`1_trimming.sh`](1_sequence_processing/1_trimming.sh)  
**Input:** demultiplexed paired FASTQ files (`*.fastq.gz`)

Cutadapt removes sequencing adapters and COI primers in three passes: (1) adapters, (2) 5′ primers, and (3) 3′ primers with final read filtering. The two target primers are:

| Primer | Sequence (5′–3′) |
| --- | --- |
| mlCOIintF | `GGWACWGGWTGAACWGTWTAYCCYCC` |
| jgHCO2198 | `TAIACYTCIGGRTGICCRAARAAYCA` |

**Key command — 5′ primer trimming (excerpt):**

```bash
cutadapt -j "$THREADS" \
    -g "X${FWD};e=0.20" -G "X${REV_I};e=0.20" \
    -O 17 --no-indels --match-read-wildcards \
    --discard-untrimmed --pair-filter=both --trim-n \
    -o "$fivep_r1" -p "$fivep_r2" \
    "$adapter_r1" "$adapter_r2"
```

`-O 17` requires at least 17 matched primer bases, while `e=0.20` permits up to 20% errors in a match. `--no-indels` disallows indels *during primer matching*. With `--pair-filter=both`, a pair is retained if **at least one** read has a recognized 5′ primer. In the final pass, `--minimum-length 50` and `--max-n 0` remove pairs with a short or ambiguous read; 3′ primer detection itself is optional. See the full script for adapter and 3′ primer commands.

**Trimming parameter considerations.** The minimum overlap setting (`-O`) can substantially affect primer detection and trimming efficiency. Users are encouraged to compare different overlap settings based on their primer sequences and sequencing characteristics. Particular attention should be paid to the proportion of reads longer than the expected amplicon by approximately the length of a primer, differences in trimming performance between R1 and R2, and reads containing only partial primer sequences, which may remain undetected under overly stringent overlap settings.

**Run:**

```bash
cd "$WORK_DIR"
sbatch 1_trimming.sh
```

**Main outputs:**

```text
1_trimming/3_primer_trim/<sample>_R1_001.fastq
1_trimming/3_primer_trim/<sample>_R2_001.fastq
1_trimming/logs/read_counts_summary.tsv
```

**Example Slurm output (`1_trimming_<jobid>.out`, excerpt):**

```text
[QC PASS] 4 unique paired-end samples
[INFO] Cutadapt 5.1; threads=12
[QC PASS] KPC0541: 115,638 pairs; below DADA2 truncLen: 54
[QC PASS] KPC0542: 149,321 pairs; below DADA2 truncLen: 107
[QC PASS] KPC0543: 129,552 pairs; below DADA2 truncLen: 79
[QC PASS] KPC0544: 224,354 pairs; below DADA2 truncLen: 140
[QC PASS] Total final pairs: 618,865
```

The “below DADA2 truncLen” values are **diagnostic counts**, not additional trimming exclusions.

---

## Step 2. ASV inference with DADA2

**Scripts:** [`2_dada2_run.sh`](1_sequence_processing/2_dada2_run.sh) and [`2_dada2.r`](1_sequence_processing/2_dada2.r)  
**Input:** primer-trimmed paired FASTQ files from Step 1

DADA2 filters reads, learns error models, infers ASVs, merges paired reads, removes chimeras, and applies a final length filter. As demonstrated in the [official DADA2 tutorial](https://benjjneb.github.io/dada2/tutorial.html), users should inspect read quality profiles before selecting truncation lengths and evaluate the learned error models afterward. Our script integrates these steps into a single workflow, and the following QC section illustrates the recommended checks using the example dataset.

<details>
<summary><strong>Recommended QC: Read quality profiles and error models</strong></summary>

> ### 1. Inspect read quality profiles
>
> Read quality should be examined separately for R1 and R2 before selecting truncation lengths. The following plots show the quality profiles of primer-trimmed reads from the four samples, prior to DADA2 quality filtering and truncation.
>
> **R1 (forward reads)**
>
> <img width="2000" height="1400" alt="01_qualityProfile_R1" src="https://github.com/user-attachments/assets/dece047a-9eca-47a7-8fd8-adf6f107cc9e" />
>
> **R2 (reverse reads)**
>
> <img width="2000" height="1400" alt="02_qualityProfile_R2" src="https://github.com/user-attachments/assets/066cfcaa-bf81-4325-9c93-2b46dff6b5cb" />
>
> R1 generally maintained high sequencing quality, whereas R2 showed a more pronounced decline toward the 3′ end. Based on these profiles, truncation lengths of 240 bp (R1) and 200 bp (R2) were used. These settings retain a theoretical overlap of 121–133 bp for the selected 307–319 bp amplicons, exceeding the minimum 100 bp required for paired-end merging.
>
> Truncation lengths should be reassessed for other datasets based on read quality, read-length distributions, and the overlap required for successful merging.
>
> ### 2. Evaluate the learned error models
>
> After error learning, the estimated error rates should be compared with the observed error frequencies using `plotErrors()`.
>
> **R1 (forward error model)**
>
> <img width="2000" height="1400" alt="03_errorModel_R1" src="https://github.com/user-attachments/assets/c13a5bbb-6d79-4978-b839-16d6353afb4a" />
>
> **R2 (reverse error model)**
>
> <img width="2000" height="1400" alt="03_errorModel_R2" src="https://github.com/user-attachments/assets/aba71298-9f42-4ef6-98ba-1c1819fe4d63" />
>
> Gray points represent observed error frequencies, black curves represent the fitted error models, and red lines indicate nominal quality-score-based error expectations.
>
> The fitted models generally captured the quality-dependent trends in observed error frequencies, although deviations were evident for some substitution types and quality-score ranges. Both models converged after four rounds, but convergence alone does not establish model adequacy. Users should therefore examine error-model plots alongside read-filtering, merging, and retention statistics before interpreting the final ASV results.
>
> **Note:** The integrated script generates these QC plots during execution. For new datasets, preliminary read quality inspection should be performed before running the full pipeline, followed by evaluation of the learned error models and processing statistics.

</details>


**Key R settings and functions (excerpts):**

```r
TRUNC_LEN <- c(240L, 200L)
MAX_EE <- c(2, 2)
TRUNC_Q <- 2L
MAX_N <- 0L

# Filtering: filterAndTrim(..., truncLen=TRUNC_LEN, maxEE=MAX_EE, ...)
# Error learning: learnErrors(..., nbases=2e8, ...)
dadaFs <- dada(filtFs, err=errF, pool="pseudo", multithread=THREADS)
dadaRs <- dada(filtRs, err=errR, pool="pseudo", multithread=THREADS)

# Pair merging: mergePairs(..., minOverlap=100L, maxMismatch=0L)
# Chimera removal: removeBimeraDenovo(..., method="consensus")
# Final length range: 307–319 bp (313 ± 6 bp)
```

`maxEE=2` limits expected errors per read, while `pool="pseudo"` shares evidence of sequence variants across samples. Paired-end merging requires at least 100 bp of overlap with no mismatches. The final length range of 313 ± 6 bp accommodates limited amplicon-length variation, including possible short COI indels. No separate singleton-removal filter was applied (no singleton ASVs were observed in this study).

**Run:**

```bash
cd "$WORK_DIR"
sbatch 2_dada2_run.sh
```

**Main outputs:**

```text
2_dada2/ASVs.nochim.fasta
2_dada2/ASV_ID_to_sequence.tsv
2_dada2/sample_by_ASV.tsv
2_dada2/logs/read_tracking.tsv
```

**Example Slurm output (`2_dada2_<jobid>.out`, excerpt):**

```text
[INFO] DADA2: 1.38.0
[STEP 3] Learn R1/R2 errors using standard LOESS
Convergence after  4  rounds.
[STEP 4] Pseudo-pooled denoising
[SUMMARY] ASVs prechim: 375 postchim: 343 length pass: 293
[SUMMARY] Final reads: 579596
[QC PASS] Final ASV tables and read tracking validated
```

**QC assessment:**

Of the 618,865 read pairs retained after primer trimming, 589,476 (95.3%) passed DADA2 quality filtering, and 579,596 (93.7%) remained after merging, chimera removal, and length selection. Read retention was consistently high across the four samples (93.1–94.6%). Although some deviations between observed and fitted error rates were evident, both error models converged after four rounds, and no substantial read loss occurred during subsequent processing. The resulting 293 ASVs, with no singletons, were therefore carried forward to the subsequence resolution step (Step 3). The final sequences and per-sample abundances are provided in `ASVs.nochim.fasta` and `sample_by_ASV.tsv`, respectively.

---

## Step 3. ASV subsequence resolution

**Script:** [`3_subsequence.sh`](1_sequence_processing/3_subsequence.sh)  
**Input:** DADA2 ASV FASTA, sequence mapping, and sample-by-ASV abundance table

In Step 2, an amplicon length range of 313 ± 6 bp was allowed to accommodate potential biological length variation, including short indels in insect COI sequences. However, this range may also retain artificially shortened ASVs that do not represent genuine biological variation. To address this before MOTU clustering, shorter ASVs were examined for exact sequence containment within longer ASVs. Automatic merging was restricted to ASVs showing 100% nucleotide identity to a unique longer parent at a terminal position and a length difference from 313 bp that was not a multiple of three. Any matches not meeting these criteria would be flagged for manual review rather than automatically merged.

**Key decision (excerpt):**

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

Only the ASV label and abundance allocation change; reads are transferred **within their original samples**. The rule is a conservative curation criterion, not proof that the shorter ASV is erroneous.

**Run:**

```bash
cd "$WORK_DIR"
sbatch 3_subsequence.sh
```

**Main outputs:**

```text
3_subsequence/asvs.subsequence_resolved.fasta
3_subsequence/asvs.subsequence_resolution.tsv
3_subsequence/sample_by_asv_resolved.tsv
3_subsequence/asvs.subsequence_summary.tsv
```

**Example Slurm output (`3_subsequence_<jobid>.out`, excerpt):**

```text
input_asvs                                       293
exact_match_child_asvs                             1
auto_merged_child_asvs                             1
auto_merged_child_reads                           15
output_asvs                                      292
output_reads                                 579,596
  ASV_99 (312 bp, 15 reads) -> ASV_89 (313 bp), final=ASV_89
[QC PASS] All ASVs accounted for; all sample read totals conserved
```

In this dataset, one shorter ASV (ASV_99, 312 bp, 15 reads) was merged into an exact-matching 313-bp ASV (ASV_89), reducing the total from 293 to 292 ASVs without changing the total read count. The resolved ASVs were then used for Swarm clustering in Step 4.

---

## Step 4. MOTU clustering with Swarm

**Script:** [`4_swarm.sh`](2_motu_curation/4_swarm.sh)  
**Input:** resolved ASV FASTA and abundance table from Step 3

Swarm clusters ASVs into MOTUs using a local nucleotide-difference threshold. ASV input headers include abundance annotations (`;size=N;`).

**Key command (excerpt):**

```bash
swarm -d "$D_VAL" -t "$THREADS" -z \
    -o "$CLUSTERS" -s "$STATS" -w "$SEEDS" "$SWARM_FASTA"
```

The default is `d=2` (up to two differences per local link), **not** a global percentage-identity cutoff. `-z` reads the abundance annotations, and `-o`, `-s`, and `-w` generate membership, statistics, and representative sequences. For `d=1`, the script also enables fastidious mode (`-f`).

**Clustering parameter considerations.** Users are encouraged to test multiple `d` values (e.g., 1–10) and compare changes in MOTU counts and cluster composition. An elbow-like pattern in MOTU counts or a marked increase in cluster merging may help identify a suitable range, although neither provides a definitive criterion. The final `d` value should be selected considering the characteristics of the target marker and dataset. This analysis used `d=2`.

**Run:**

```bash
cd "$WORK_DIR"
sbatch 4_swarm.sh 2
```

**Main outputs (`d=2`):**

```text
4_swarm_motu/ALL_d2/ALL_d2_MOTU.fasta
4_swarm_motu/ALL_d2/ALL_d2_MOTU_table.tsv
4_swarm_motu/ALL_d2/ALL_d2_members.tsv
4_swarm_motu/ALL_d2/ALL_d2_summary.tsv
```

**Example Slurm output (`4_swarm_<jobid>.out`, excerpt):**

```text
[INFO] Swarm: Swarm 3.1.4
[INFO] Distance d: 2
input_ASVs                                  292
MOTUs                                       268
singleton_ASV_MOTUs                         259
multi_ASV_MOTUs                               9
output_reads                            579,596
[QC PASS] All ASVs assigned once; every sample and all reads conserved
```

Here, `singleton_ASV_MOTUs` means **MOTUs containing one ASV**, not MOTUs with one read.

---

## Step 5. MOTU curation with VSEARCH and LULU

**Script:** [`5_lulu.sh`](2_motu_curation/5_lulu.sh)  
**Input:** Swarm representative FASTA and MOTU abundance table

VSEARCH identifies similar MOTU representative sequences; LULU then uses sequence similarity, relative abundance, and sample co-occurrence to merge candidate daughter MOTUs into parents.

**Key commands (excerpts):**

```bash
vsearch --usearch_global "$FASTA" --db "$FASTA" \
    --strand plus --self --id 0.97 --iddef 1 --query_cov 0.90 \
    --maxaccepts 0 --maxrejects 0 \
    --threads "$THREADS" --userout "$MATCHLIST" \
    --userfields query+target+id
```

```r
result <- lulu::lulu(
  otu, matches,
  minimum_ratio_type="min", minimum_ratio=1,
  minimum_match=97, minimum_relative_cooccurence=0.95
)
```

VSEARCH uses a **97% identity** threshold and at least **90% query coverage** to identify candidates. LULU requires matching similarity, a parent-to-child abundance ratio of at least 1 under its selected rule, and ≥95% relative co-occurrence.

**Run:**

```bash
cd "$WORK_DIR"
sbatch 5_lulu.sh 2
```

**Main outputs (`d=2`):**

```text
5_lulu/ALL_d2/LULU_curated_MOTU.fasta
5_lulu/ALL_d2/LULU_curated_MOTU_table.tsv
5_lulu/ALL_d2/LULU_MOTU_map.tsv
5_lulu/ALL_d2/LULU_summary.tsv
```

**Example Slurm output (`5_lulu_<jobid>.out`, excerpt; progress percentages omitted):**

```text
[INFO] LULU version: 0.1.0
[QC] VSEARCH match records: 16
input_MOTUs                          268
retained_MOTUs                       267
merged_MOTUs                         1
input_reads                          579596
curated_reads                        579596
read_conservation                    PASS
[QC PASS] All MOTUs mapped; curated counts and sample reads conserved
```

---

## Step 6. Taxonomic assignment with MMseqs2

**Script:** [`6_mmseqs_lca4.sh`](3_taxonomic_assignment/6_mmseqs_lca4.sh)  
**Input:** LULU-curated MOTU representatives and a taxonomy-enabled nucleotide reference database

Because this study required species- and genus-level taxonomic information to investigate potential host associations, assignments obtained using the 97% identity threshold were prioritized for interpretation. However, many MOTUs remained unclassified under this condition. An additional analysis at 80% identity was therefore performed to explore the broader taxonomic affiliations of these MOTUs.

**Key command (excerpt; `$identity` is 0.80 or 0.97):**

```bash
mmseqs taxonomy "$QUERY_DB" "$REF_DB" "$result" "$tmp" \
    --tax-lineage 1 --threads "$THREADS" \
    --search-type 3 --orf-filter 0 --lca-mode 4 \
    --alignment-mode 3 --strand 2 \
    -s 5 -e 1e-5 --min-seq-id "$identity" \
    --min-aln-len 300 -c 0.80 --cov-mode 2 \
    --max-accept 20 --max-rejected 50
```

`--lca-mode 4` assigns taxonomy based on the lowest common ancestor of equally scoring top hits and was applied to both identity conditions. `--min-aln-len 300` requires a long match to the ~313 bp target; `-c 0.80 --cov-mode 2` sets query coverage. `-s 5` controls search sensitivity, and `--max-accept 20` limits the number of accepted search hits.

**Run:**

```bash
cd "$WORK_DIR"
sbatch 6_mmseqs_lca4.sh
```

**Main outputs:**

```text
6_mmseqs_lca4/results/taxonomy_lca4_a20_80.tsv
6_mmseqs_lca4/results/taxonomy_lca4_a20_97.tsv
6_mmseqs_lca4/logs/taxonomy_qc.tsv
```

**Example Slurm output (`6_mmseqs_lca4_<jobid>.out`, excerpt; search progress omitted):**

```text
[INFO] MMseqs2: 13.45111
[QC PASS] Curated input: 267 MOTUs, 579,596 reads
[STEP 2] Create query database (267 MOTUs)
[QC PASS] Query database created
[QC PASS] 80: 267 unique taxonomy rows; all MOTUs represented
[QC PASS] 97: 267 unique taxonomy rows; all MOTUs represented
[QC PASS] Both taxonomy assignments validated for 267 MOTUs
```

---

## Step 7. Comparing 80% and 97% assignments

**Script:** [`7_taxonomy_comparison.py`](3_taxonomic_assignment/7_taxonomy_comparison.py)  
**Input:** both MMseqs2 taxonomy TSVs and the LULU-curated abundance table

This Python stage joins the two sets of taxonomic assignments **by MOTU ID** (not row position), preserves sample-level abundances, and labels changes between the conditions.

**Key comparison logic (excerpt):**

```python
def compare(t80, t97):
    a, b = t80[0], t97[0]
    if a == 0 and b == 0:
        return "Unclassified both"
    if a != 0 and b == 0:
        return "Lost at 97%"
    if a == 0 and b != 0:
        return "Gained at 97%"
    return "Same TaxID" if a == b else "Changed TaxID"
```

A `Lost at 97%` assignment means that the stricter search did not return a TaxID; it does **not** mean the MOTU or its reads were removed. The output XLSX is a plain worksheet with no custom styling or filters.

**Run (directly in the terminal; no Slurm submission):**

```bash
cd "$WORK_DIR"
python3 7_taxonomy_comparison.py
```

**Main outputs:**

```text
7_taxonomy_comparison/MOTU_comparison.tsv
7_taxonomy_comparison/MOTU_comparison.xlsx
7_taxonomy_comparison/logs/comparison_summary.tsv
```

**Example terminal output (excerpt):**

```text
[QC PASS] Input IDs agree: 267 MOTUs, 579,596 reads
[QC PASS] Total MOTUs: 267; reads: 579,596
[QC] Insecta: 80% = 40; 97% = 9
[QC] Same TaxID: 48
[QC] Lost at 97%: 163
[QC] Gained at 97%: 0
[QC] Changed TaxID: 0
[QC] Unclassified both: 56
```

---

## Example results

| Stage | Retained units | Reads |
| --- | ---: | ---: |
| DADA2 | 293 ASVs | 579,596 |
| Subsequence resolution | 292 ASVs | 579,596 |
| Swarm (`d=2`) | 268 MOTUs | 579,596 |
| LULU | 267 MOTUs | 579,596 |

The displayed Slurm and terminal excerpts are taken **from this dataset's run logs**, with progress messages and insignificant whitespace omitted. For troubleshooting or detailed QC, inspect the complete log and output tables for the relevant stage.

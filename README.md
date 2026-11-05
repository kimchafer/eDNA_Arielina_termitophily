[README.md](https://github.com/user-attachments/files/34422758/README.md)
# eDNA_Arielina_termitophily
Analysis pipeline for surface-soil eDNA metabarcoding to infer putative termite hosts of Arielina beetles.
# Body-surface eDNA metabarcoding pipeline for *Arielina peringueyi*

This repository documents a step-by-step workflow for processing short COI amplicons recovered from material adhering to the body surface of *Arielina peringueyi*. The workflow was developed to screen for potentially associated organisms, including candidate termite hosts. Taxonomic assignments provide candidates for subsequent validation; they do not, by themselves, establish host identity or a biological association.

The analysis starts with demultiplexed Illumina paired-end FASTQ files and ends with a comparison of taxonomic assignments at two minimum sequence-identity thresholds. Scripts are organized by analytical stage and include input validation, output checks, and QC reports.

**Workflow:**

```text
Demultiplexed paired-end FASTQ
    |
    v
1. Cutadapt: adapter and COI primer trimming
    |
    v
2. DADA2: quality filtering, denoising, paired-end merging,
   chimera removal, and ASV inference
    |
    v
3. ASV subsequence resolution
    |
    v
4. Swarm: MOTU clustering
    |
    v
5. VSEARCH + LULU: MOTU curation
    |
    v
6. MMseqs2: taxonomic assignment (80% and 97% identity)
    |
    v
7. Taxonomic comparison: combined TSV and Excel tables
```

## Contents

- [Repository organization](#repository-organization)
- [Requirements](#requirements)
- [Input data and setup](#input-data-and-setup)
- [Step 1. Adapter and primer trimming](#step-1-adapter-and-primer-trimming)
- [Step 2. ASV inference with DADA2](#step-2-asv-inference-with-dada2)
- [Step 3. ASV subsequence resolution](#step-3-asv-subsequence-resolution)
- [Step 4. Swarm clustering](#step-4-swarm-clustering)
- [Step 5. LULU curation](#step-5-lulu-curation)
- [Step 6. Taxonomic assignment with MMseqs2](#step-6-taxonomic-assignment-with-mmseqs2)
- [Step 7. Compare the two assignment thresholds](#step-7-compare-the-two-assignment-thresholds)
- [Example analysis summary](#example-analysis-summary)
- [Interpretation and reproducibility](#interpretation-and-reproducibility)

## Repository organization

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
├── 3_taxonomic_assignment/
│   ├── 6_mmseqs_lca4.sh
│   └── 7_taxonomy_comparison.py
├── Fig2A/
├── Fig2B/
└── Fig2C/
```

The numbered folders contain the source scripts. The intermediate and final analysis outputs are generated in a separate working directory and are **not** included in the repository.

## Requirements

The workflow was written for a Linux HPC environment using the [Slurm](https://slurm.schedmd.com/documentation.html) job scheduler. It uses:

| Stage | Software |
| --- | --- |
| 1 | Cutadapt, Python 3 |
| 2 | R, DADA2 |
| 3 | Python 3 |
| 4 | Swarm, Python 3 |
| 5 | VSEARCH, R, LULU |
| 6 | MMseqs2 with a nucleotide taxonomy database |
| 7 | Python 3 standard library |

Python **3.9 or later** is recommended for compatibility with the supplied scripts. The shell scripts activate pre-existing Conda environments at paths/names used on the original HPC system; these are not environment-definition files. Install the required software and change the environment activation commands for your system before running the scripts. Each stage records relevant software versions or parameters in its output logs.

**Before running:** the scripts contain absolute paths to the original project's working directory, raw FASTQ directory, Conda installation, and MMseqs2 reference database. These paths, as well as the Slurm partition and resource requests, **must be adapted for a different computing environment**. Several stages protect existing output directories from being overwritten.

## Input data and setup

The original dataset consisted of four demultiplexed Illumina paired-end libraries (`KPC0541`–`KPC0544`). The target was an approximately 313-bp COI fragment amplified using:

- **mlCOIintF**: `5'-GGWACWGGWTGAACWGTWTAYCCYCC-3'`
- **jgHCO2198**: `5'-TAIACYTCIGGRTGICCRAARAAYCA-3'`

The raw paired FASTQ files should be placed in the configured `RAW_DIR`, with matching `_R1_001.fastq.gz` and `_R2_001.fastq.gz` names. For example:

```text
raw_data/
├── KPC0541_E_Lib_S99_L001_R1_001.fastq.gz
├── KPC0541_E_Lib_S99_L001_R2_001.fastq.gz
├── KPC0542_E_Lib_S100_L001_R1_001.fastq.gz
├── KPC0542_E_Lib_S100_L001_R2_001.fastq.gz
└── ...
```

Download the source code:

```bash
git clone https://github.com/kimchafer/eDNA_Arielina_termitophily.git
cd eDNA_Arielina_termitophily
```

The scripts expect to find each other's outputs under a shared `WORK_DIR`. For the original HPC setup, the relevant paths were:

```bash
WORK_DIR=/storage2/jihoonkim/Arielina_eDNA/analyses/1_20261103
RAW_DIR=/storage2/jihoonkim/Arielina_eDNA/raw_data
```

For a new analysis, create your own working and raw-data directories and copy the eight scripts into your working directory. For example, from the repository root:

```bash
REPO_DIR="$PWD"
WORK_DIR=/path/to/your/analysis
mkdir -p "$WORK_DIR"
cp "$REPO_DIR/1_sequence_processing/"{1_trimming.sh,2_dada2_run.sh,2_dada2.r,3_subsequence.sh} "$WORK_DIR/"
cp "$REPO_DIR/2_motu_curation/"{4_swarm.sh,5_lulu.sh} "$WORK_DIR/"
cp "$REPO_DIR/3_taxonomic_assignment/"{6_mmseqs_lca4.sh,7_taxonomy_comparison.py} "$WORK_DIR/"
```

Then edit the path configuration in each copied script. In particular, `2_dada2_run.sh` expects `2_dada2.r` directly in `WORK_DIR`, rather than in a subdirectory. When retaining the published wrapper names, also check the optional `SCRIPT` checksum path in Steps 1 and 2, which refers to an earlier filename. **Do not create stage output directories in advance**, because the scripts use overwrite-protection checks.

The following commands assume that the scripts have already been copied to `WORK_DIR`, that the raw data and reference database are available, and that the absolute paths and software environments have been configured. Run each stage **only after the preceding stage has finished successfully**.

## Step 1. Adapter and primer trimming

**Script:** [`1_trimming.sh`](1_sequence_processing/1_trimming.sh)

Trim Illumina adapters, then the 5′ and 3′ COI primers with Cutadapt. The workflow requires at least one 5′ primer match per read pair and applies the final minimum-length (50 bp) and ambiguous-base (`N`) filters to both reads. Paired-read synchronization and length distributions are checked.

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

Inspect the trimming report before proceeding:

```bash
head -n 6 1_trimming/logs/read_counts_summary.tsv
```

The final paired FASTQ files serve as the input to Step 2.

## Step 2. ASV inference with DADA2

**Scripts:** [`2_dada2_run.sh`](1_sequence_processing/2_dada2_run.sh) and [`2_dada2.r`](1_sequence_processing/2_dada2.r)

The Slurm wrapper launches the R pipeline. Reads are quality-filtered (`truncLen=240,200`; `maxEE=2,2`), denoised with a learned error model and pseudo-pooling, merged (`minOverlap=100`, `maxMismatch=0`), and checked for chimeras. The final candidate sequences are selected within the **307–319 bp** length range.

```bash
cd "$WORK_DIR"
sbatch 2_dada2_run.sh
```

**Main outputs:**

```text
2_dada2/ASVs.nochim.fasta
2_dada2/ASV_ID_to_sequence.tsv
2_dada2/sample_by_ASV.tsv
2_dada2/seqtab.final.rds
2_dada2/logs/read_tracking.tsv
```

`ASVs.nochim.fasta` contains the final length-selected, nonchimeric ASVs despite its shorter filename. Inspect read tracking and the final ASV abundance table:

```bash
head -n 6 2_dada2/logs/read_tracking.tsv
head -n 6 2_dada2/ASV_ID_to_sequence.tsv
```

## Step 3. ASV subsequence resolution

**Script:** [`3_subsequence.sh`](1_sequence_processing/3_subsequence.sh)

This step searches for 100% exact containment of shorter ASVs within longer ASVs. It conservatively transfers reads when a shorter sequence matches a single longer parent exclusively at a sequence end and meets the script's length-class rule. Ambiguous relationships are retained for manual review rather than automatically merged.

```bash
cd "$WORK_DIR"
sbatch 3_subsequence.sh
```

**Main outputs:**

```text
3_subsequence/asvs.subsequence_resolved.fasta
3_subsequence/sample_by_asv_resolved.tsv
3_subsequence/asvs.subsequence_resolution.tsv
3_subsequence/asvs.subsequence_manual_review.tsv
3_subsequence/asvs.subsequence_summary.tsv
```

```bash
cat 3_subsequence/asvs.subsequence_summary.tsv
```

This is an abundance **redistribution** step. Global and per-sample read totals should remain unchanged.

## Step 4. Swarm clustering

**Script:** [`4_swarm.sh`](2_motu_curation/4_swarm.sh)

Cluster the resolved ASVs into molecular operational taxonomic units (MOTUs) using [Swarm](https://github.com/torognes/swarm). The default clustering distance is **`d=2`**; `d=1` uses fastidious mode. The input ASVs and their abundances remain traceable through an ASV-to-MOTU membership table.

```bash
cd "$WORK_DIR"
sbatch 4_swarm.sh 2
```

**Main outputs for `d=2`:**

```text
4_swarm_motu/ALL_d2/ALL_d2_MOTU.fasta
4_swarm_motu/ALL_d2/ALL_d2_MOTU_table.tsv
4_swarm_motu/ALL_d2/ALL_d2_members.tsv
4_swarm_motu/ALL_d2/ALL_d2_summary.tsv
```

```bash
cat 4_swarm_motu/ALL_d2/ALL_d2_summary.tsv
```

Other distances can be tested in separate output directories (e.g., `sbatch 4_swarm.sh 1` or `sbatch 4_swarm.sh 3`).

## Step 5. LULU curation

**Script:** [`5_lulu.sh`](2_motu_curation/5_lulu.sh)

Generate pairwise MOTU candidate matches with VSEARCH and curate the Swarm clusters using LULU. The reference settings use a **97% sequence-identity candidate threshold**, **90% query coverage**, `minimum_ratio_type="min"`, `minimum_ratio=1`, and `minimum_relative_cooccurence=0.95`.

```bash
cd "$WORK_DIR"
sbatch 5_lulu.sh 2
```

**Main outputs:**

```text
5_lulu/ALL_d2/LULU_curated_MOTU.fasta
5_lulu/ALL_d2/LULU_curated_MOTU_table.tsv
5_lulu/ALL_d2/LULU_MOTU_map.tsv
5_lulu/ALL_d2/LULU_summary.tsv
```

```bash
cat 5_lulu/ALL_d2/LULU_summary.tsv
```

LULU redistributes daughter-MOTU abundances into retained parent MOTUs without removing reads. The parent mapping records these changes.

## Step 6. Taxonomic assignment with MMseqs2

**Script:** [`6_mmseqs_lca4.sh`](3_taxonomic_assignment/6_mmseqs_lca4.sh)

Classify the **LULU-curated MOTU representative sequences** against the configured nucleotide taxonomy database using MMseqs2 in **LCA mode 4**. Two independent analyses differ only in the minimum sequence-identity threshold: **0.80** and **0.97**. Other important parameters include `--min-aln-len 300`, `-c 0.80`, `--cov-mode 2`, and `--max-accept 20`.

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

```bash
head -n 5 6_mmseqs_lca4/results/taxonomy_lca4_a20_80.tsv
head -n 5 6_mmseqs_lca4/results/taxonomy_lca4_a20_97.tsv
```

A minimum-identity threshold is a **search criterion**, not a reported per-query alignment identity. An unclassified result under the 97% setting does not establish that the underlying sequence is non-insect or absent from the sample.

## Step 7. Compare the two assignment thresholds

**Script:** [`7_taxonomy_comparison.py`](3_taxonomic_assignment/7_taxonomy_comparison.py)

Join the two MMseqs2 taxonomy tables by MOTU ID and add the LULU-curated sample abundances. The script produces a **21-column** comparison, including taxonomic assignment, rank, NCBI TaxID, lineage, insect classification, and whether the assignment is retained or lost at the more stringent identity setting.

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

The Excel workbook intentionally contains a **plain, unformatted worksheet**, without custom coloring, borders, table styles, filters, or frozen panes.

```bash
cat 7_taxonomy_comparison/logs/comparison_summary.tsv
```

## Example analysis summary

For the four-library analysis used to develop this workflow, the main sequence-count checkpoints were:

| Processing stage | Retained sequence units | Total reads |
| --- | ---: | ---: |
| DADA2 | 293 ASVs | 579,596 |
| Subsequence resolution | 292 ASVs | 579,596 |
| Swarm (`d=2`) | 268 MOTUs | 579,596 |
| LULU | 267 MOTUs | 579,596 |

At minimum sequence identity 80%, MMseqs2 assigned 40 of the 267 curated MOTUs within Insecta; at 97%, nine were assigned within Insecta. These counts describe the **specific analysis shown here**, not values that should be hard-coded for other datasets.

## Interpretation and reproducibility

- **Low-abundance detections require caution.** Taxonomic assignments from a small number of reads are candidate signals and should be interpreted with extraction/PCR controls, sample handling, geographic context, and independent validation in mind.
- **MOTUs are not automatically biological species.** Swarm clustering and LULU curation follow sequence-similarity and abundance/co-occurrence rules; they do not establish species boundaries.
- **Read conservation is a computational check.** Steps 3–5 redistribute abundance among sequence units. An unexpected change in total or per-sample reads indicates a problem with accounting or data processing.
- **No MACSE step is part of this published workflow.** The final analysis consists of the seven stages above.
- **Raw sequencing reads and the nucleotide taxonomy database are not distributed in this code repository.** The MMseqs2 reference database must be obtained or built separately, and its identity/version should be documented for reproducibility.
- **These scripts are tailored to a particular HPC environment.** Adjust absolute paths, environments, Slurm resources, sample names where relevant, and output directories before reusing them.

Each script contains further details under `PURPOSE`, `INPUT`, `OUTPUT`, and `SETTINGS`. Consult those source files as the authoritative record of the implemented analysis parameters.

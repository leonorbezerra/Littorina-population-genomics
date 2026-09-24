# Thesis-2026: Population Genomics and Pool-Seq Analysis in *Littorina*

> ⚠️ **IMPORTANT NOTICE:** These scripts are **NOT** "plug-and-play" pipelines. They are heavily customized codes and templates. To use them, you **MUST manually open the scripts and edit** paths, local file names, chromosome naming conventions, and sample-specific biological parameters (e.g., pool sizes and coverage depth).

This repository contains tools and pipelines for short-read variant processing, sliding-window population genomics statistics, and downstream visualization for marine snails of the genus *Littorina*. 

The workflow for individual analysis is adapted from *A standard pipeline for processing short-read sequencing data from Littorina snails V.3* (James Reeve et al.).
The analysis in pools use a sliding-window genomic paradigm to analyze allele frequencies (AF), nucleotide diversity (π), Tajima's D, and differentiation (\(F_{ST}\)).

---

## Repository Contents

*   `run_gatk+filters_popgen_pipeline.sh` - Main Bash script for alignment processing, GATK HaplotypeCaller variant calling, and `bcftools` quality filtering. Commands a priori for adapter removal and mapping followed this pipeline
*   `poolseq-popgen-analysis.sh` - Script to execute multi-population sliding-window calculations (AF, Pi, Tajima's D, and \(F_{ST}\)) for pool data.
*   `plot_poolseq_popgen_stats.R` - R script for genome-wide and chromosome-specific Manhattan plots and statistical tests.
*   `individual_genomic_analyses.txt` - Core comprehensive code workflow for downstream analyses (LD Heatmaps, Regional/Window PCAs, GWAS via GEMMA, pi and \(d_{XY}\) calculations, and Phylogenetic Trees), including codes used both in R and in Bash (Linux)
*   `calc_window_heterozygosity.R` - R script using `SeqArray` to compute individual observed heterozygosity in sliding windows.

*   `LICENSE` - Open-source MIT License.

---

## Requirements

### Command Line Tools
The scripts require the following command-line tools to be installed:

*   **Trimmomatic** (PE read trimming)
*   **BWA** (Genome indexing and `mem` mapping)
*   **Samtools** & **Bcftools** (Alignment sorting/indexing, pileup generation, and VCF filtering)
*   **Picard Tools** (`MarkDuplicates` and `AddOrReplaceReadGroups`)
*   **PoPoolation1** (`Variance-sliding.pl` for window-based \(\pi\) and Tajima's D)
*   **PoPoolation2** (`mpileup2sync.jar`, `fst-sliding.pl`, and `snp-frequency-diff.pl`)
*   **GATK4** (`HaplotypeCaller`, `CombineGVCFs`, `GenotypeGVCFs`, and `GatherVcfs`)
*   **PLINK 1.9** (Linkage Disequilibrium matrix and dataset pruning)
*   **GEMMA** (Genome-wide Association Studies via LMM)
*   **IQ-TREE 3** (Phylogenetic tree inference: ML and BIONJ)

### R Packages (CRAN & Bioconductor)
The scripts require the following R packages to be installed:

```r
# From CRAN
install.packages(c("ggplot2", "data.table", "dplyr", "tidyr", 
                   "readr", "readxl", "stringr", "patchwork", 
                   "plotly", "qqman", "ggtext"))

# From Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("SeqArray", "SNPRelate", "vcfR", "adegenet"))
```

## Usage Guide


### 1. `poolseq-popgen-analysis.sh` (Pool-Seq Architecture)
**What it does:**
Processes raw FASTQ data from pools. It trims adapters, maps reads, strips PCR duplicates, and estimates average coverage. It then feeds the resulting alignments into the PoPoolation environment to calculate windowed Nucleotide Diversity (\(\pi\)), Tajima's D, sliding-window \(F_{ST}\), and allele frequency differences.

**What MUST be adapted inside the script:**
*   **Executable Paths and file names:** All `path/to/trimmomatic.jar`, `path/to/picard.jar`, the paths to PoPoolation `.pl` scripts or `.jar` tools, and all the output and inpute file names.
*   **Biological & Sequencing Specs:**
    *   `--pool-size` (PoPoolation1 and PoPoolation2) must be adjusted to the exact number of haploid chromosomes inside your pools.
    *   `--max-coverage` must be recalculated as \(2.5\times\) your pool coverage average.


### 2. `plot_poolseq_popgen_stats.R` (Metric Plotting & Local Statistics)
**What it does:**
Loads external files for \(F_{ST}\), AF, and Tajima's D, combines chromosome boundaries to map a single continuous genome track in Megabases, and runs localized Wilcoxon tests between inversion zones and background chromosomes.

**What MUST be adapted inside the script:**
*   **Input Files:** file names like `"input.fst"`, `"input.af_pwc"`, and `"input.D"` must be changed accordingly
*   **Coordinates:** Chromosome boundaries are customized specifically for the target inversion blocks on Chromosome 11 on the Littorina pools used during this work

### 3. `run_gatk+filters_popgen_pipeline.sh` (Individual Variant Calling)
**What it does:**
Iterates through individual duplicate-removed BAM files to add Read Groups via `Picard`, calls genomic variants (`-ERC GVCF`) per chromosome using GATK `HaplotypeCaller`, merges individual blocks (`CombineGVCFs` / `GenotypeGVCFs`), and passes the raw VCF through a strict multi-stage `bcftools` filtering cascade.

**What MUST be adapted inside the script:**
*   **Paths:** The `REF`, `BAMS`, `FINAL_VCF_DIR`, and `LOCAL_TMP` variables.
*   **Region Parameters:** The `$GATK_REGION` and `$CHR` variables used to restrict `HaplotypeCaller` targets.
*   **Chromosome Identifiers:** Specific ENA accession numbers (`'ENA|OZ249416|OZ249416.1'`) mutch match the vcf header.
*   **Coverage Filter:** The depth expression `INFO/DP>2.5xAverage Depth` is a placeholder text line. It must be replaced with the real value.


### 4. `individual_genomic_analyses.txt` (Advanced Analytical Notebook)
**What it does:**
Combines distinct bash and R blocks to evaluate Linkage Disequilibrium heatmaps (`PLINK`), regional & window-based PCAs via `adegenet` and `SNPRelate` packages, nucleoide diversity and divergence (\(d_{XY}\) by winows, Mixed-Model Association Maps (`GEMMA`), and alignments for maximum-likelihood tree structures (`IQ-TREE`).

**What MUST be adapted inside the script:**
*   **Local Spreadsheets:** Commands like `read_excel("samples.xlsx")` and `read.csv("excel.csv")` require local metadata sheets containing sample IDs matched with color phenotypes.
*   **Variables:** Variable definitions like `putative_homo_1`, sample series, and target matrix offsets must be manually rewritten to map your sample directory.

  

### 5. `calc_window_heterozygosity.R` (Observed Heterozygosity Arrays)
**What it does:**
Opens a processed genomic GDS file using the `SeqArray` engine, and evaluates heterozygosity values in 250kb windows.

**What MUST be adapted inside the script:**
*   **GDS Target:** The `seqOpen("chr11")` target file name.
*   **Window Bounds:** The `start_pos`, `end_pos`, and `window_size` variables that define the loop array.
---

## 📄 License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

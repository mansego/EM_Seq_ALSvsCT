# Multi-Omics Analysis of cfDNA in ALS
**Universitat Autònoma de Barcelona (UAB)** 
**Bioinformatics Master's Thesis – 2025**

## 📌 Overview 
This repository contains bioinformatics pipelines and analysis scripts for the integrative study of **cell-free DNA (cfDNA)** from **Amyotrophic Lateral Sclerosis (ALS)** patients and healthy controls, using **Enzymatic Methyl-seq (EM-seq)** technology. The project aims to identify **ALS-associated epigenetic and fragmentation biomarkers**, investigate the **tissue origin of cfDNA**, and explore **epigenetic differences between ALS clinical phenotypes**, specifically **bulbar-onset vs spinal-onset ALS**. Analyses include genome-wide methylation profiling, fragmentomics, deconvolution of cfDNA origin, and clinical data integration.

Key analytical components include:

- **Genome-wide profiling of differentially methylated loci and regions (DMLs/DMRs) associated with ALS and ALS clinical phenotypes**.

- **Profilling cfDNA Fragmentation Patterns associated with ALS**.

- **Tissue-of-origin deconvolution using supervised and unsupervised reference-free methods**.

- **Clinical and phenotypic data integration**.


## 📂 Directory Structure  

```
EM_seq_ALSvsCT/
├── PrePro/               # Preprocessing with nf-core/methylseq (FastQC, Bismark, Trim Galore, MultiQC)
├── fragment_analysis/    # Fragment size distributions and GC bias analysis
├── DMA/                  # Differential Methylation Analysis (BSseq, DSS)
├── Deconvolution/        # Tissue-of-origin estimation (UXM, DeconvR, HiBED)
├── pData/                # Clinical and phenotypic metadata + exploratory statistics
├── references/           # Genome reference files (GRCh38, SNPs, annotations)
├── data/                 # Raw input data (FASTQ)
├── Suplementary_data.xlsx  # Supplementary tables summarizing key results (DMLs, DMRs, sample metrics, cell-type composition)
└── README.md             # Main project documentation
```

## 🛠️ Installation & Dependencies  
1. **Environment Setup**: 
Use `mamba` (recommended over conda for performance):
  ```bash
conda config --add channels conda-forge
conda config --set channel_priority strict
conda install mamba

mamba create -n tfm_als -c bioconda -c conda-forge -c defaults \
  samtools picard  r-ggplot2 r-tidyverse python nextflow nf-core singularity
conda activate tfm_als
   ```
2. **R/Python Libraries**: 
  ```r
  # R packages
   install.packages(c("ggplot2", "tidyverse", "ggpubr","UpSetR", "pheatmap"))  # R  
   pip install seaborn pandas numpy     
   ```
```bash
# Python packages
pip install pandas numpy seaborn 
   ```

## 🚀 Quick Start

## 📁 Preprocessing (nf-core/methylseq with Bismark for EM-seq)

```bash
cd PrePro/
nextflow run nf-core/methylseq -profile singularity \
  -params-file params.yaml \
  -resume
```
EM-seq mode is activated with `--em_seq true`, using GRCh38 as reference genome.
Output includes: `bismark/`, `fastqc/`, `multiqc/`, `trimgalore/`

## 📊 Fragment Size & GC Bias Analysis
 ```bash
# Run fragment length extraction  
bash scripts/fragmentomics/01_extract_lengths.sh  
# Generate plots  
Rscript scripts/fragmentomics/02_plot_distributions.R  
```
## 📈 Differential Methylation Analysis (DMA)
Genome-wide differential methylation analysis was conducted on EM-seq cfDNA data using the DSS package. Two comparisons were explored:

- ALS vs Control samples

- Bulbar-onset vs Spinal-onset ALS phenotypes

After rigorous quality control (SNP exclusion, coverage filtering, batch correction), the following results were obtained:
|Comparison	|Significant DMLs	|Significant DMRs|
|-----------|-----------------|----------------|
|ALS vs Control|	69,749	|182|
|Bulbar vs Spinal ALS	|233	|1|
Functional enrichment (GO, KEGG, DisGeNET) highlighted neuronal and ALS-related pathways. Scripts and outputs are organized in the [`DMA/`](DMA) folder and include annotated DMLs/DMRs, volcano plots, PCA, and enrichment summaries.

## 🧬 Deconvolution of Tissue-of-Origin
Methods used:

- **UXM** (Unsupervised reference-free)

- **DeconvR** (Reference-free)

- **HiBED** (Supervised for cfDNA tissues)
```bash
Rscript Deconvolution/Deconvolution.R
```

Output:

- `Deconvolution/tables/`: Cell-type estimates, statistics

- `Deconvolution/figures/`: Heatmaps, barplots by phenotype (e.g., bulbar vs spinal)

## 🧾 Clinical & Phenotypic Metadata
Contains demographic and clinical annotations, and basic statistical comparisons between groups (ALS vs Control, Bulbar vs Spinal):

Files:

- ALSvsCT_clinical.html/pdf

- pData_ALSvsCT.csv, pData_BvsS.csv

## 📄 Supplementary Data

The file Suplementary_data.xlsx contains seven worksheets that summarize key quantitative results supporting the main analyses:
| Table | Title | Description |
|-------|-------|-------------|
| **Supplementary Table 1** | Sample Information and EM-Seq Summary | Contains metadata and summary metrics from the Enzymatic Methyl-seq protocol: input DNA amount, library yield, conversion rate, etc. |
| **Supplementary Table 2** | Alignment and Methylation Metrics | Includes sequencing and methylation statistics per sample: total reads, mapping efficiency, duplication rate, CpG coverage, and global methylation percentage. |
| **Supplementary Table 3** | ALS vs Control DMLs | Differentially methylated loci (DMLs) between ALS and control samples (FDR < 0.05 and Δβ > 0.1). |
| **Supplementary Table 4** | ALS vs Control DMRs | Differentially methylated regions (DMRs) between ALS and controls (P < 1e-5 and Δβ > 0.1). |
| **Supplementary Table 5** | Bulbar vs Spinal DMLs | DMLs between bulbar-onset and spinal-onset ALS phenotypes (FDR < 0.05 and Δβ > 0.1). |
| **Supplementary Table 6** | Bulbar vs Spinal DMRs | DMRs between ALS phenotypes (FDR < 0.05 and Δβ > 0.1). |
| **Supplementary Table 8** | Cell-Type Composition Estimates | Comparative analysis of estimated cfDNA cell-type proportions between ALS and control groups, and within ALS phenotypes (bulbar vs spinal), using deconvolution methods. |

All statistical results are derived from genome-wide methylation analysis and cfDNA deconvolution approaches. Session information and software versions used to generate these results are included in supplementary text files (`session_info_*.txt`) to ensure reproducibility.

All results were generated from EM-Seq data of cfDNA, processed using Bismark, BSseq, DSS and deconvolution tools (UXM, HiBED, DeconvR).

For reproducibility, corresponding session information is included in:

- `fragment_analysis/session_info_FA.txt`

- `fragment_analysis/session_info_GC.txt`

- `Deconvolution/session_info.txt`
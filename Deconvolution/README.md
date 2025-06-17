# Cell Type Deconvolution from cfDNA Methylation Data

This project focuses on **cell type deconvolution** of cfDNA methylation profiles from EM-seq data. The objective is to estimate the proportions of cell types contributing to plasma cfDNA in ALS patients and controls using three complementary tools:

- **Deconv**: R-based reference deconvolution [BIMSBbioinfo/deconvR](https://github.com/BIMSBbioinfo/deconvR)
- **UXM**: Unsupervised deconvolution of methylation patterns using fragment-level WGBS data [nloyfer/UXM_deconv](https://github.com/nloyfer/UXM_deconv)
- **HiBED**: Hierarchical brain cell-type deconvolution [SalasLab/HiBED](https://github.com/SalasLab/HiBED)


## UXM Workflow

### 1. Environment Setup

```bash
# Create and activate Conda environment
conda create -n UXM python=3.10 pandas=1.5 numpy scipy samtools htslib bedtools -c conda-forge -c bioconda -y
conda activate UXM
```
### 2. Install wgbs_tools

```bash
cd /mnt/mydisk/EM_SEQ_ALSvsCT/Deconvolution/UXM_deconv/
git clone https://github.com/nloyfer/wgbs_tools.git
cd wgbs_tools
python setup.py install
export PATH="${PATH}:$PWD"
```
### 3. Genome Initialization
```bash
wgbs_tools.py init_genome hg38
```
### 4. Generate `.pat.gz` and `beta` files from Deduplicated BAM files

To convert deduplicated and sorted BAM files to `.pat.gz` and beta format using `wgbstools`, run the following loop:

```bash
mkdir -p logs/bam2pat

for sample in BACW_42 BACW_44 BACW_45 BACW_47 BACW_48 BACW_50 BACW_52 BACW_53 \
              BACW_55 BACW_56 BACW_57 BACW_58 BACW_59 BACW_61 BACW_64 BACW_65; do
  echo "Processing $sample..."
  wgbstools bam2pat ../../PrePro/bismark/deduplicated/${sample}.deduplicated.sorted.bam \
    > logs/bam2pat/${sample}.log 2>&1

  if [[ $? -ne 0 ]]; then
    echo "❌ Error processing $sample. Check logs/bam2pat/${sample}.log"
  else
    echo "✅ Done with $sample"
  fi
done

mkdir ref_dir
mv *beta ref_dir/
 
    

```

### 5. Install and Configure UXM

```bash
git clone https://github.com/nloyfer/UXM_deconv.git
cd ../UXM_deconv/
export PATH="${PATH}:$PWD"
```
### 6. Run UXM Deconvolution
```bash
uxm deconv --atlas supplemental/Atlas.U25.l4.hg38.full.tsv ../wgbs_tools/*pat.gz -o output.csv
```

### 7. Plot Results
```bash
uxm plot output.csv -o output.pdf

```
⚠️ If `uxm plot` fails, use the following conda environment:

```bash
conda create -n UXM2 python pandas=1.5.3 matplotlib scipy samtools tabix bedtools -y
conda activate UXM2
```
## R Unified Deconvolution Workflow: DeconvR, HiBED, and UXM

The [`Deconvolution.R`](Deconvolution.R)  script performs a **comprehensive pipeline** for cell-type deconvolution using EM-seq data. It integrates:

### 📋 Main Steps:

1. **Mapping** of cfDNA EM-seq signal to EPIC probe space using `bsseq`.
2. **Deconvolution**:
   - `deconvR`: Uses the `HumanCellTypeMethAtlas` reference.
   - `HiBED`: Brain-specific hierarchical model with `HiBED_Libraries`.
   - `UXM`: Uses previously generated `output.csv` from UXM (unsupervised).
3. **Statistical comparisons**:
   - ALS vs Control (`Condition`)
   - Bulbar vs Spinal within ALS (`ALS_Phenotype`)
4. **Visualization**:
   - Boxplots for significantly different cell types.
5. **Export**:
   - Tables in Excel `.xlsx` format.
   - PNG plots per tool and contrast.
   - Session info log.

### 📂 Input Requirements

| File | Description |
|------|-------------|
| `pData_ALSvsCT.csv` | Metadata with `Sample`, `Condition`, `ALS_Phenotype` columns |
| `bsseq_chr*.rds` | Chromosome-wise BSseq objects (chr1–22) |
| `UXM_deconv/output.csv` | Output file from UXM deconvolution |

### 📁 Output Structure

All outputs are saved in the `Deconvolution/` directory:

```swift
Deconvolution/
├── Deconvolution.R # Main R script
├── session_info.txt # R environment info
├── data/
│ └── meth_matrix.rda # Probe-mapped methylation matrix
├── figures/
│ ├── *_boxplot.png # Boxplots per tool/contrast
├── tables/
│ ├── *_combined_results.xlsx # Summary statistics for each tool
│ ├── merged_deconv_results.tsv # Merged table across methods
├── UXM_deconv/
│ └── output.csv # Input for UXM analysis
```

### 📜 Usage

Place the script in the project folder and run:

```r
source("Deconvolution.R")
```

📎 Note: The UXM output must be available as UXM_deconv/output.csv before running this script.
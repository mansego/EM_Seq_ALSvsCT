# cfDNA Fragmentomics Analysis Pipeline

**Objective**: Quantify cfDNA fragment size distributions and identify ALS-specific patterns.

**Input**:
- Deduplicated BAM files (`*.deduplicated.sorted.bam`).
- Sample metadata (ALS vs. Controls).

**Tools**:
- `samtools`, `Picard`, `R` (tidyverse, ggplot2).

**Steps**:

 ## 0.Set Up Environment

 ```bash
 # Directory structure
mkdir -p fragment_analysis/raw_lengths fragment_analysis/plots fragment_analysis/stats

# List of samples (ALS vs. Controls)
ALS_SAMPLES=("BACW_42" "BACW_45" "BACW_47" "BACW_50" "BACW_53" "BACW_56" "BACW_58" "BACW_65")
CTRL_SAMPLES=("BACW_44" "BACW_48" "BACW_52" "BACW_55" "BACW_57" "BACW_59" "BACW_61" "BACW_64")
```
##  1.Extract Fragment Lengths from BAMs

`samtools` is used to extract insert sizes (ignoring unmapped/mate-unmapped reads):

 ```bash
for sample in "${ALS_SAMPLES[@]}" "${CTRL_SAMPLES[@]}"; do
  samtools view -@ 8 -F 0x404 -q 20 "/mnt/mydisk/EM_Seq_ALSvsCT/PrePro/bismark/deduplicated/${sample}.deduplicated.sorted.bam" | 
    awk '$9 > 0 {print $9}' > "fragment_analysis/raw_lengths/${sample}.lengths.txt"
done

for sample in "${ALS_SAMPLES[@]}" "${CTRL_SAMPLES[@]}"; do
samtools view -@ 8 -F 0x404 -q 20 "/mnt/mydisk/EM_Seq_ALSvsCT/PrePro/bismark/deduplicated/${sample}.deduplicated.sorted.bam" | 
awk '!($3 ~ /^chrUn_|_alt$|_random$/) && $9 > 0' |
bioawk -c fastx '{if(gsub(/[Nn]/, "", $seq)/length($seq) < 0.1) print $9}' > "fragment_analysis/raw_lengths/${sample}.lengths.txt"
done
 ```

## 2.Calculate Fragment Statistics
`Picard Tools` is used to generate fragment length metrics
```bash
for sample in "${ALS_SAMPLES[@]}" "${CTRL_SAMPLES[@]}"; do
  java -jar picard.jar CollectInsertSizeMetrics \
    I="/mnt/mydisk/EM_Seq_ALSvsCT/PrePro/bismark/deduplicated/${sample}.deduplicated.sorted.bam" \
    O="fragment_analysis/stats/${sample}_insert_metrics.txt" \
    H="fragment_analysis/plots/${sample}_insert_histogram.pdf" \
    M=0.5
done
```
**Outputs**:

- `*_insert_size_metrics.txt`: Mean/median fragment size, % short/long fragments.

- `*_insert_size_histogram.pdf`: Visual distribution.


## 3.Aggregate and Analyze Data (R)

R Script (`fragment_analysis/fragment_analysis.R`):

### 3.1.Merge Fragment Lengths into a Single DataFrame


```r
library(tidyverse)

# Read all length files
als_samples <- c("BACW_42", "BACW_45", "BACW_47", "BACW_50", "BACW_53", "BACW_56", "BACW_58", "BACW_65")
ctrl_samples <- c("BACW_44", "BACW_48", "BACW_52", "BACW_55", "BACW_57", "BACW_59", "BACW_61", "BACW_64")

base_dir <- "fragment_analysis/raw_lengths"

als_files <- file.path(base_dir, paste0(als_samples, ".lengths.txt"))
ctrl_files <- file.path(base_dir, paste0(ctrl_samples, ".lengths.txt"))


read_lengths <- function(files, group) {
  map_dfr(files, ~ read_tsv(.x, col_names = "length") %>% 
           mutate(sample = basename(.x)), .id = "file_id") %>%
    mutate(group = group)
}

als_data <- read_lengths(als_files, "ALS")
ctrl_data <- read_lengths(ctrl_files, "Control")
combined_data <- bind_rows(als_data, ctrl_data)

# Save merged data
write_tsv(combined_data, "fragment_analysis/stats/combined_fragment_lengths.tsv")
```
### 3.2. Statistical Comparisons
```r
# Summary statistics
summary_stats <- combined_data %>%
  group_by(group) %>%
  summarise(
    mean_length = mean(length),
    median_length = median(length),
    sd_length = sd(length),
    prop_short = sum(length < 150) / n(),
    prop_long = sum(length > 250) / n()
  )

# Kolmogorov-Smirnov test
ks_test <- ks.test(
  x = filter(combined_data, group == "ALS")$length,
  y = filter(combined_data, group == "Control")$length
)

# Save results
write_tsv(summary_stats, "fragment_analysis/stats/summary_stats.tsv")
sink("fragment_analysis/stats/ks_test.txt")
print(ks_test)
sink()
```
### 3.3. Visualization

**Density Plot (Group Comparison)**:
```r
ggplot(combined_data, aes(x = length, fill = group)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = c(166, 320), linetype = "dashed", color = "gray40") +
  labs(title = "cfDNA Fragment Length Distribution",
       x = "Fragment Length (bp)",
       y = "Density",
       fill = "Group") +
  theme_minimal()
ggsave("fragment_analysis/plots/fragment_length_density.png", width = 8, height = 6)
```
**Boxplot (Short vs. Long Fragments)**:
```r
combined_data %>%
  mutate(category = case_when(
    length < 150 ~ "Short (<150 bp)",
    length > 250 ~ "Long (>250 bp)",
    TRUE ~ "Medium"
  )) %>%
  filter(category != "Medium") %>%
  ggplot(aes(x = group, y = length, fill = category)) +
  geom_boxplot() +
  labs(title = "Short vs. Long cfDNA Fragments by Group")
ggsave("fragment_analysis/plots/short_vs_long_boxplot.png")
#Short/Long Fragment Ratios
ratio_plot <- summary_stats %>%
  pivot_longer(cols = c(prop_short, prop_long), names_to = "category") %>%
  ggplot(aes(x = group, y = value, fill = category)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Proportion of Short/Long cfDNA Fragments", 
       y = "Proportion", 
       x = "Group") +
  theme_minimal()
ggsave("fragment_analysis/plots/short_long_ratio.png", ratio_plot, width = 8, height = 6)
```
## Integration with Other Modules
### 1. GC Bias Correction:

- Run *after* fragment length extraction to normalize coverage biases.
```bash
java -jar picard.jar CollectGcBiasMetrics \
   I="/mnt/mydisk/EM_Seq_ALSvsCT/PrePro/bismark/deduplicated/${sample}.deduplicated.sorted.bam" \
   O="fragment_analysis/stats/${sample}_gc_bias_metrics.txt" \
   CHART="fragment_analysis/plots/${sample}_gc_bias_chart.pdf" \
   S="fragment_analysis/stats/${sample}_gc_summary.txt" \
   R="../references/Homo_sapiens/NCBI/GRCh38/Sequence/BismarkIndex/genome.fa"
```
### 2. Nucleosome Positioning:

- Use fragment length data (combined_data) to detect periodicity:
```r
spectrum <- spectrum(combined_data$length, plot = FALSE)
nucleosome_period <- 1 / spectrum$freq[which.max(spectrum$spec)]  # Expected: ~166 bp

```
### 3. Downstream Analyses:

- Feed fragmentomics results into DMA (e.g., filter DMRs by fragment size) or deconvolution (e.g., cell-type-specific fragmentation).

## Expected ALS-Specific Signals (TODO)
- [x] Increased short fragments (<150 bp) in ALS due to apoptotic debris.
- [x] Altered nucleosome peaks (e.g., 166 bp nucleosomal vs. 320 bp subnucleosomal) suggesting chromatin remodeling.

## Final Outputs
- fragment_analysis/stats/:
  - summary_stats.tsv (mean/median/short-long ratios).
  - ks_test.txt (group differences).
- fragment_analysis/plots/:
  - Density plots, boxplots, histograms.
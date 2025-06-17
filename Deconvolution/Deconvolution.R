# -------------------------------------------------------------
# Unified Deconvolution Analysis: deconvR, UXM and HiBED
# Author: ML Mansegp 
# Description: Cell-type deconvolution using reference-based methods
# -------------------------------------------------------------

# Load libraries
library(bsseq)
library(GenomicRanges)
library(deconvR)
library(HiBED)
library(tidyverse)
library(writexl)

# Load phenotype data
pheno <- read_csv("/mnt/mydisk/EM_Seq_ALSvsCT/pData/pData_ALSvsCT.csv") %>%
  select(Sample, Condition, ALS_Phenotype)

# Directories
input_dir <- "/mnt/mydisk/EM_Seq_ALSvsCT/DMA/RDS/NoFilter/"
output_dir <- "/mnt/mydisk/EM_Seq_ALSvsCT/Deconvolution/"
dir.create(output_dir, showWarnings = FALSE)

# -------------------------------------------------------------
# Function
# -------------------------------------------------------------
run_deconv_stats <- function(props_df, pheno_df, contrast_col, group1, group2, prefix, outdir) {
  df_long <- props_df %>%
    pivot_longer(cols = -c(Sample, Condition, ALS_Phenotype),
                 names_to = "CellType", values_to = "Proportion") %>%
    filter(!is.na(Proportion))
  
  res <- df_long %>%
    group_by(CellType) %>%
    summarise(
      n_group1 = sum(!!sym(contrast_col) == group1),
      n_group2 = sum(!!sym(contrast_col) == group2),
      median_group1 = median(Proportion[!!sym(contrast_col) == group1], na.rm = TRUE),
      median_group2 = median(Proportion[!!sym(contrast_col) == group2], na.rm = TRUE),
      p_value = ifelse(n_group1 >= 2 & n_group2 >= 2,
                       wilcox.test(Proportion[!!sym(contrast_col) == group1],
                                   Proportion[!!sym(contrast_col) == group2],
                                   exact = FALSE)$p.value,
                       NA_real_)
    ) %>%
    mutate(FDR = p.adjust(p_value, method = "fdr")) %>%
    select(CellType,
           !!paste0("median_", group1) := median_group1,
           !!paste0("median_", group2) := median_group2,
           p_value, FDR)

  #write_csv(res, file.path(outdir, paste0(prefix, ".csv")))
  return(res)
}

plot_deconv_boxplots <- function(df_long, contrast_col, out_prefix, outdir, signif_only = FALSE, stats_table = NULL) {
  df_long <- df_long %>%
    filter(!is.na(Proportion), !is.na(!!sym(contrast_col))) %>%
    mutate(Grupo = !!sym(contrast_col))

  # Filter only cell types significativos 
  if (signif_only && !is.null(stats_table)) {
    signif_celltypes <- stats_table %>%
      filter(!is.na(p_value) & p_value < 0.1) %>%
      pull(CellType)
    df_long <- df_long %>% filter(CellType %in% signif_celltypes)
  }

  p <- ggplot(df_long, aes(x = Grupo, y = Proportion, fill = Grupo)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.15, size = 1.2, alpha = 0.6) +
    facet_wrap(~ CellType, scales = "free_y") +
    theme_minimal(base_size = 11) +
    labs(title = paste("Cell-type proportions by", contrast_col),
         x = NULL, y = "Proportion") +
    theme(strip.text = element_text(face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")

  ggsave(filename = file.path(outdir, paste0(out_prefix, "_boxplot.png")),
         plot = p, width = 12, height = 8, dpi = 300)
}

# -------------------------------------------------------------
# Step 1: Map EM-seq data to EPIC probes
# -------------------------------------------------------------
data("HumanCellTypeMethAtlas")
data("IlluminaMethEpicB5ProbeIDs")
data("HiBED_Libraries")

message("Mapping EM-Seq data to EPIC probes...")
all_mapped <- list()
for (chr in 1:22) {
  message(sprintf("Processing chromosome %d", chr))
  bsseq_chr <- readRDS(file.path(input_dir, paste0("bsseq_chr", chr, ".rds")))
  beta <- getMeth(bsseq_chr, type = "raw")
  beta[getCoverage(bsseq_chr) < 0] <- NA
  gr <- granges(bsseq_chr)
  mcols(gr) <- as.data.frame(beta)
  all_mapped[[chr]] <- BSmeth2Probe(probe_id_locations = IlluminaMethEpicB5ProbeIDs, WGBS_data = gr, multipleMapping = TRUE,cutoff = 10 ) #,cutoff = 10
}
meth_matrix <- do.call(rbind, all_mapped)
save(meth_matrix, file = file.path(output_dir, "data/meth_matrix.rda"))
#load(file.path(output_dir, "data/meth_matrix.rda"))

# -------------------------------------------------------------
# 1. deconvR Analysis
# -------------------------------------------------------------
message("Running deconvR...")
deconvR_props <- as.data.frame(deconvolute(reference = HumanCellTypeMethAtlas, bulk = meth_matrix)$proportions)
deconvR_props$Sample <- rownames(deconvR_props)
deconvR_props <- left_join(deconvR_props, pheno, by = "Sample")

stats_deconvR <- run_deconv_stats(deconvR_props, pheno, "Condition", "ALS", "Control", "deconvR_ALSvsCT", output_dir)
stats_deconvR_pheno <- run_deconv_stats(deconvR_props %>% filter(Condition == "ALS"), pheno, "ALS_Phenotype", "Bulbar", "Spinal", "deconvR_ALS_phenotype", output_dir)


write_xlsx((full_join(stats_deconvR, stats_deconvR_pheno, by = "CellType", suffix = c("_ALS_CT", "_Bulbar_Spinal"))), file.path(output_dir, "deconvR_combined_results.xlsx"))

deconvR_long <- deconvR_props %>%
  pivot_longer(cols = -c(Sample, Condition, ALS_Phenotype),
               names_to = "CellType", values_to = "Proportion")

plot_deconv_boxplots(deconvR_long, "Condition", "deconvR_ALS_CT", output_dir, signif_only = TRUE, stats_table = stats_deconvR)
plot_deconv_boxplots(filter(deconvR_long, Condition == "ALS"), "ALS_Phenotype", "deconvR_Bulbar_Spinal", output_dir, signif_only = TRUE, stats_table = stats_deconvR_pheno)

# -------------------------------------------------------------
# 2. HiBED Analysis
# -------------------------------------------------------------
message("Running HiBED...")
rownames(meth_matrix) <- meth_matrix[, 1]
meth_matrix <- meth_matrix[, -1]
HiBED_props <- as.data.frame(HiBED_deconvolution(meth_matrix, h = 2))
HiBED_props$Sample <- rownames(HiBED_props)
HiBED_props <- left_join(HiBED_props, pheno, by = "Sample")

stats_HiBED <- run_deconv_stats(HiBED_props, pheno, "Condition", "ALS", "Control", "HiBED_ALSvsCT", output_dir)
stats_HiBED_pheno <- run_deconv_stats(HiBED_props %>% filter(Condition == "ALS"), pheno, "ALS_Phenotype", "Bulbar", "Spinal", "HiBED_ALS_phenotype", output_dir)


write_xlsx((full_join(stats_HiBED, stats_HiBED_pheno, by = "CellType", suffix = c("_ALS_CT", "_Bulbar_Spinal"))), file.path(output_dir, "HiBED_combined_results.xlsx"))

HiBED_long <- HiBED_props %>%
  pivot_longer(cols = -c(Sample, Condition, ALS_Phenotype),
               names_to = "CellType", values_to = "Proportion")

plot_deconv_boxplots(HiBED_long, "Condition", "HiBED_ALS_CT", output_dir, signif_only = TRUE, stats_table = stats_HiBED)
plot_deconv_boxplots(filter(HiBED_long, Condition == "ALS"), "ALS_Phenotype", "HiBED_Bulbar_Spinal", output_dir, signif_only = TRUE, stats_table = stats_HiBED_pheno)

# -------------------------------------------------------------
# 3. UXM Analysis
# -------------------------------------------------------------
message("Processing UXM...")
uxm <- read_csv(file.path(output_dir, "UXM_deconv/output.csv"))
uxm_long <- uxm %>%
  pivot_longer(cols = -CellType, names_to = "Sample", values_to = "Proportion") %>%
  left_join(pheno, by = "Sample")

uxm_props <- uxm_long %>%
  pivot_wider(names_from = CellType, values_from = Proportion) %>%
  distinct(Sample, .keep_all = TRUE)

stats_uxm <- run_deconv_stats(uxm_props, pheno, "Condition", "ALS", "Control", "UXM_ALSvsCT", output_dir)
stats_uxm_pheno <- run_deconv_stats(uxm_props %>% filter(Condition == "ALS"), pheno, "ALS_Phenotype", "Bulbar", "Spinal", "UXM_ALS_phenotype", output_dir)

write_xlsx((full_join(stats_uxm, stats_uxm_pheno, by = "CellType", suffix = c("_ALS_CT", "_Bulbar_Spinal"))), file.path(output_dir, "UXM_combined_results.xlsx"))

plot_deconv_boxplots(uxm_long, "Condition", "UXM_ALS_CT", output_dir, signif_only = TRUE, stats_table = stats_uxm)
plot_deconv_boxplots(filter(uxm_long, Condition == "ALS"), "ALS_Phenotype", "UXM_Bulbar_Spinal", output_dir, signif_only = TRUE, stats_table = stats_uxm_pheno)



# ======================
#  SESSION INFO
# ======================

writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"))

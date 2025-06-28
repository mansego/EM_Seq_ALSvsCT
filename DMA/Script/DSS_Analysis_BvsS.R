#!/usr/bin/env Rscript

# ======================
# 1. INITIAL SETUP
# ======================

# Load required packages with installation if needed
required_packages <- c(
  "DSS", "bsseq", "annotatr", "pheatmap", "GenomicRanges",
  "data.table", "parallel", "magrittr", "dplyr", "ggplot2", 
  "stringr", "qvalue", "patchwork","ggpubr","limma","ggrepel"
)

invisible({
  suppressPackageStartupMessages({
    for (pkg in required_packages) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(pkg, repos = "https://cloud.r-project.org")
      }
      library(pkg, character.only = TRUE)
    }
  })
})

# ======================
# 2. CONFIGURATION
# ======================


# File paths
config <- list(
  work_dir = "/mnt/mydisk/EM_Seq_ALSvsCT/DMA/RDS/Filter_batch_SNP_cov10_BvsS_batch_onset",
  base_path = "/mnt/mydisk/EM_Seq_ALSvsCT/PrePro/bismark/methylation_call/methylation_coverage",
  snp_file = "/mnt/mydisk/EM_Seq_ALSvsCT/references/snp151Common.txt.gz",
  output_dir = "/mnt/mydisk/EM_Seq_ALSvsCT/DMA/Results/BulbarvsSpinal_ALS"
)

# Create directories if they don't exist
#dir.create(config$work_dir, showWarnings = FALSE, recursive = TRUE)
#dir.create(config$output_dir, showWarnings = FALSE, recursive = TRUE)

# Sample metadata
sample_meta <- read.csv("/mnt/mydisk/EM_Seq_ALSvsCT/pData/pData_ALSvsCT.csv", stringsAsFactors = FALSE) %>%
  filter(Condition == "ALS") %>%
  mutate(
    id = c(42, 45, 47, 50, 53, 56, 58, 65),
    name = c("B42", "B45", "S47", "B50", "S53","B56", "S58", "S65"),
    group = ifelse(grepl("^B", name), "Bulbar", "Spinal"),
    group = factor(group, levels = c("Spinal", "Bulbar")),
    #cov_file = file.path(config$base_path, paste0("BACW_", id, "_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz"))
  )

# Analysis parameters
params <- list(
  chromosomes = paste0("chr", 1:22),
  maf_cutoff = 0.10,
  min_coverage = 10,
  min_samples = 4,
  dmr_params = list(
    delta = 0.1,
    p.threshold = 1e-5,
    minlen = 50,
    minCG = 3,
    dis.merge = 100,
    pct.sig = 0.5
  )
)

# ======================
# 3. FUNCTIONS
# ======================

#' Process and filter SNPs by MAF
process_snps <- function(snp_file, maf_cutoff = 0.10) {
  message("Loading and filtering SNPs...")
  snp_df <- fread(
    snp_file, 
    header = FALSE, 
    select = c(2:5, 9, 23, 25),
    col.names = c("chr", "start", "end", "rsID","refalt", "alleles", "frequencies")
  )
  
  freqs <- tstrsplit(gsub(",$", "", snp_df$frequencies), ",", type.convert = TRUE)
  snp_df[, MAF := pmin(freqs[[1]], freqs[[2]], na.rm = TRUE)]
  snp_df <- snp_df[!is.na(MAF) & MAF >= maf_cutoff]
  snp_df$start1based <- snp_df$start + 1

  snp_gr <- GRanges(
    seqnames = snp_df$chr,
    ranges = IRanges(start = snp_df$start1based, width = 1)
  )
  seqlevelsStyle(snp_gr) <- "NCBI"
  
  return(snp_gr)
}



process_chromosome <- function(chr, sample_meta, snp_gr, min_coverage = 10, min_samples = 4) {
  message("\n===== Processing chromosome ", chr, " =====")
  
  # 1. Load Coverage Data
  message("\n[1/6] Loading coverage data...")
  current_chr <- chr
  cov_list <- lapply(sample_meta$cov_file, function(f) {
    message("  Reading file: ", basename(f))
    dat <- fread(f, header = FALSE,
                 col.names = c("chr", "start", "end", "meth_pct", "meth", "unmeth"))
    dat_chr <- dat[chr == current_chr, .(chr, pos = start, X = meth, N = meth + unmeth)]
    if (nrow(dat_chr) == 0) {
      message("    - No data for ", chr, " in this file")
      return(NULL)
    }
    message("    - ", nrow(dat_chr), " CpGs found")
    return(dat_chr)
  })
  
  # Load statistics
  cov_list <- cov_list[!sapply(cov_list, is.null)]
  message("\nLoad Summary:")
  message("- Files processed: ", length(sample_meta$cov_file))
  message("- Files with data: ", length(cov_list))
  
  if (length(cov_list) == 0) {
    message("\n[!] No data for ", chr, ", skipping...")
    return(NULL)
  }
  
  # 2. Create BSseq Object
  message("\n[2/6] Creating BSseq object...")
  bsseq_obj <- makeBSseqData(cov_list, sampleNames = sample_meta$name) %>%
    keepSeqlevels(value = chr, pruning.mode = "coarse")
  message("- Total CpGs: ", length(bsseq_obj))
  
  # 3. Filter SNPs
  message("\n[3/6] Filtering SNPs...")
  cpg_gr <- granges(bsseq_obj)
  seqlevelsStyle(cpg_gr) <- "NCBI"
  overlap <- findOverlaps(cpg_gr, snp_gr)
  message("- CpGs overlapping with SNPs: ", length(unique(queryHits(overlap))))
  bsseq_filtered <- bsseq_obj[-queryHits(overlap), ]
  message("- CpGs after SNP filtering: ", length(bsseq_filtered))
  
  # 4. Coverage Filtering
  message("\n[4/6] Applying coverage filters...")
  cov <- getCoverage(bsseq_filtered, type = "Cov")
  keep <- rowSums(cov[, sample_meta$group == "Bulbar"] >= min_coverage) >= min_samples &
          rowSums(cov[, sample_meta$group == "Spinal"] >= min_coverage) >= min_samples
  message("- CpGs removed due to low coverage: ", sum(!keep))
  bsseq_filtered <- bsseq_filtered[keep, ]
  message("- Remaining CpGs: ", length(bsseq_filtered))
  
  # 5. Batch Effect Correction
  message("\n[5/6] Applying batch effect correction...")
  Batch <- sample_meta$Batch
  Onset <- sample_meta$Onset_sympton
  meth <- getMeth(bsseq_filtered, type = "raw")
  
  # Calculate beta and M-values
  beta <- pmin(pmax(meth, 1e-5), 1 - 1e-5)
  Mvalues <- log2(beta / (1 - beta))
  
  # Identify rows with finite values
  finite_rows <- rowSums(is.finite(Mvalues)) == ncol(Mvalues)
  message("- CpGs with non-finite values: ", sum(!finite_rows))
  Mvalues_finite <- Mvalues[finite_rows, ]
  
  if (nrow(Mvalues_finite) == 0) {
    message("\n[!] No finite values after filtering for ", chr)
    return(NULL)
  }
  
  # Apply correction
  design <- model.matrix(~ 0 + group, data = sample_meta)
  Mvalue_corrected <- limma::removeBatchEffect(Mvalues_finite, batch = Batch, covariates=Onset,design = design)
  
  # Convert back to beta values
  beta_corrected <- 2^Mvalue_corrected / (1 + 2^Mvalue_corrected)
  beta_corrected[is.nan(beta_corrected)] <- 0
  message("- CpGs after correction: ", nrow(beta_corrected))
  
  # 6. Reconstruct BSseq Object
  message("\n[6/6] Reconstructing corrected BSseq object...")
  bsseq_finite <- bsseq_filtered[finite_rows, ]
  
  bsseq_corrected <- BSseq(
    chr = as.character(seqnames(bsseq_finite)),
    pos = start(bsseq_finite),
    M = beta_corrected * getCoverage(bsseq_finite),
    Cov = getCoverage(bsseq_finite),
    sampleNames = sampleNames(bsseq_finite)
  )
  message("- Reconstructed BSseq object with ", length(bsseq_corrected), " CpGs")
  
  # Perform DML Test
  message("\nPerforming DML test...")
  dml_test <- DMLtest(
    bsseq_filtered,
    group1 = which(sample_meta$group == "Bulbar"),
    group2 = which(sample_meta$group == "Spinal"),
    smoothing = TRUE
  )
  message("- DML test completed for ", nrow(dml_test), " CpGs")
  
  # Save Results
  save_files <- list(
    # Mvalue = file.path(config$work_dir, paste0("Mvalue_", chr, "_filtered.rds")),
    bsseq = file.path(config$work_dir, paste0("bsseq_", chr, "_filtered.rds")),
    dmlTest = file.path(config$work_dir, paste0("dmlTest_", chr, "_filtered.rds"))
  )
  
  message("\nSaving results:")
  # saveRDS(Mvalue_corrected, save_files$Mvalue)
  # message("- M-values saved to ", save_files$Mvalue)
  saveRDS(bsseq_corrected, save_files$bsseq)
  message("- BSseq object saved to ", save_files$bsseq)
  saveRDS(dml_test, save_files$dmlTest)
  message("- DML results saved to ", save_files$dmlTest)
  
  # Return summary
  list(
    chromosome = chr,
    CpGs_original = length(bsseq_obj),
    CpGs_removed_SNPs = length(unique(queryHits(overlap))),
    CpGs_removed_coverage = sum(!keep),
    CpGs_removed_nonfinite = sum(!finite_rows),
    CpGs_remaining = length(bsseq_corrected),
    files = save_files
  )
}


generate_summary_report <- function(results, output_dir) {
  # Convert results to data frame
  summary_df <- bind_rows(results)
  
  # Create report object
  report <- list(
    header = "METHYLATION PROCESSING REPORT",
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    parameters = list(
      min_coverage = params$min_coverage,
      min_samples = params$min_samples,
      chromosomes = params$chromosomes
    ),
    summary = summary_df,
    totals = list(
      total_CpGs_original = sum(summary_df$CpGs_original),
      total_CpGs_remaining = sum(summary_df$CpGs_remaining),
      retention_rate = round(sum(summary_df$CpGs_remaining) / sum(summary_df$CpGs_original) * 100, 2)
    )
  )
  
  # Save report
  report_file <- file.path(output_dir, "processing_report.json")
  writeLines(jsonlite::toJSON(report, pretty = TRUE), report_file)
  message("\nReport saved to: ", report_file)
  
  # Retention plot
  retention_plot <- ggplot(summary_df, aes(x = chromosome, y = CpGs_remaining/CpGs_original)) +
    geom_col(fill = "steelblue") +
    geom_text(aes(label = paste0(round(CpGs_remaining/CpGs_original*100, 1), "%")), 
              vjust = -0.5) +
    labs(title = "CpG Retention by Chromosome",
         y = "Retention Percentage", x = "Chromosome") +
    theme_minimal()
  
  plot_file <- file.path(output_dir, "retention_plot.png")
  ggsave(plot_file, retention_plot, width = 8, height = 6)
  message("Retention plot saved to: ", plot_file)
  
  return(report)
}


#' Simplify annotation types for visualization
simplify_annot_types <- function(annot_types) {
  annot_types %>%
    str_replace("hg38_genes_", "") %>%
    str_replace("hg38_cpg_", "CpG ") %>%
    recode(
      "3UTRs" = "3'UTRs",
      "5UTRs" = "5'UTRs",
      "exons" = "Exons",
      "introns" = "Introns",
      "promoters" = "Promoters",
      "intergenic" = "Intergenic",
      "CpG inter" = "CpG islands",
      "CpG shores" = "CpG shores",
      "CpG shelves" = "CpG shelves",
      "CpG islands" = "Inter-CpG-islands"
    )
}

#' Create annotation pie chart
plot_annotation_pie <- function(annot_df, title) {
  annot_df %>%
    count(annot.type) %>%
    mutate(
      annot_type = simplify_annot_types(annot.type),
      pct = round(n / sum(n) * 100, 1),
      label = sprintf("%s: %.1f%%", annot_type, pct)
    ) %>%
    ggplot(aes(x = "", y = n, fill = label)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar(theta = "y") +
    theme_void() +
    ggtitle(title) +
    theme(
      plot.margin = margin(10, 50, 10, 10),
      legend.text = element_text(size = 10)
    )
}

#' Plot fold enrichment
plot_fold_enrichment <- function(gr, annotations, title) {
  # Get observed annotations
  annot_df <- annotate_regions(regions = gr, annotations = annotations, ignore.strand = TRUE) %>% 
    as.data.frame()
  
  # Calculate observed percentages
  annot_table <- table(annot_df$annot.type)
  annot_pct <- annot_table / sum(annot_table)
  
  # Calculate expected percentages based on genome-wide feature lengths
  feature_lengths <- sapply(split(annotations, annotations$type), function(x) sum(width(x)))
  feature_pct <- feature_lengths / sum(feature_lengths)
  
  # Calculate fold enrichment
  fold_change <- annot_pct / feature_pct[names(annot_pct)]
  
  data.frame(
    annot_type = names(fold_change),
    fold_enrichment = as.numeric(fold_change),
    stringsAsFactors = FALSE
  ) %>%
    mutate(annot_type = simplify_annot_types(annot_type)) %>%
    arrange(desc(fold_enrichment)) %>%
    ggplot(aes(x = reorder(annot_type, fold_enrichment), y = fold_enrichment)) +
    geom_bar(stat = "identity", fill = "plum", width = 0.7) +
    coord_flip() +
    labs(
      title = title,
      x = "",
      y = "Fold enrichment (observed / expected)"
    ) +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 10))
}

#' Summarize DML/DMR results
summarize_regions <- function(df, diff_col, type = c("DML", "DMR")) {
  type <- match.arg(type)
  
  # Split into hyper and hypo
  hyper <- df[df[[diff_col]] > 0, ]
  hypo <- df[df[[diff_col]] < 0, ]
  
  # Common stats
  get_stats <- function(data) {
    genes <- unique(unlist(strsplit(na.omit(data$annot.symbol), ";")))
    genes <- genes[genes != "NA"]
    
    list(
      Count = nrow(data),
      Genes = length(genes)
    )
  }
  
  hyper_stats <- get_stats(hyper)
  hypo_stats <- get_stats(hypo)
  
  # Create summary data frame
  summary_df <- data.frame(
    Type = c(paste0("Hyper-", type, "s"), paste0("Hypo-", type, "s")),
    Count = c(hyper_stats$Count, hypo_stats$Count),
    Genes = c(hyper_stats$Genes, hypo_stats$Genes),
    stringsAsFactors = FALSE
  )
  
  # Add DMR-specific stats if needed
  if (type == "DMR") {
    summary_df$Length <- c(
      sprintf("%.0f (%.0f-%.0f)", median(hyper$length), min(hyper$length), max(hyper$length)),
      sprintf("%.0f (%.0f-%.0f)", median(hypo$length), min(hypo$length), max(hypo$length))
    )
    summary_df$CpGs <- c(
      sprintf("%.1f (%.1f-%.1f)", median(hyper$nCG), min(hyper$nCG), max(hyper$nCG)),
      sprintf("%.1f (%.1f-%.1f)", median(hypo$nCG), min(hypo$nCG), max(hypo$nCG))
    )
  }
  
  return(summary_df)
}
# ======================
# 4. MAIN ANALYSIS
# ======================

setwd(config$work_dir)

# 4.1 Process SNPs
snp_gr <- process_snps(config$snp_file, params$maf_cutoff)
saveRDS(snp_gr, file.path(config$work_dir, paste0("snp_gr_MAF", params$maf_cutoff, ".rds")))
#load(file.path(config$work_dir, paste0("snp_gr", params$maf_cutoff, ".rds")))

# 4.2 Process each chromosome

message("Starting processing of ", length(params$chromosomes), " chromosomes...")
results <- lapply(params$chromosomes, function(chr) {
  message("\n", strrep("=", 50))
  process_chromosome(chr, sample_meta, snp_gr, params$min_coverage, params$min_samples)
}) 

# Filter null results
results <- Filter(Negate(is.null), results)

# Generate summary report
if (length(results) > 0) {
  report <- generate_summary_report(results, config$output_dir)
  message("\nProcessing completed successfully!")
} else {
  message("\n[!] Warning: No chromosomes were processed successfully")
}

# 4.3 Quality Control (QC) Analyses
message("\n===== Running QC Analyses =====")

# Load all filtered BSseq objects and combine into a single BSseq object
bsseq_files <- list.files(pattern = "^bsseq_.*_filtered\\.rds$")
bsseq_all <- lapply(bsseq_files, readRDS)
bsseq_merged <- do.call(BiocGenerics::rbind, bsseq_all)
# saveRDS(bsseq_merged, file.path(config$work_dir, "bsseq_merged_filtered.rds"))

# 4.3.1. Global Methylation Distribution
meth <- getMeth(bsseq_merged, type = "raw")
meth_melted <- reshape2::melt(meth)
colnames(meth_melted) <- c("CpG", "Sample", "Beta")
meth_melted$Group <- sample_meta$group[match(meth_melted$Sample, sample_meta$name)]

p_density <- ggplot(meth_melted, aes(x = Beta, color = Group)) +
  geom_density(alpha = 0.5) +
  labs(title = "Global Methylation Distribution by Group") +
  theme_bw()

ggsave(file.path(config$output_dir, "QC_global_methylation_density.png"), p_density)

# 4.3.2. Principal component analysis (PCA)

epsilon <- 1e-5
Mvals <- log2(pmin(pmax(meth, epsilon), 1 - epsilon) / (1 - pmin(pmax(meth, epsilon), 1 - epsilon)))

# Filter rows and columns with valid variance
row_ok <- apply(Mvals, 1, function(x) all(is.finite(x)) && sd(x) > 0)
col_ok <- apply(Mvals[row_ok, ], 2, function(x) sd(x) > 0)
Mvals_clean <- Mvals[row_ok, col_ok]

# PCA
pca <- prcomp(t(Mvals_clean), scale. = TRUE)

pca_df <- data.frame(
  Sample = colnames(Mvals_clean),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  Group = sample_meta$group
)

p_pca<-ggplot(pca_df, aes(PC1, PC2, color = Group)) +
  stat_ellipse(aes(group = Group), type = "norm", linetype = 2, linewidth = 1, level = 0.95) +
  geom_point(size = 4, alpha = 0.7) +
  geom_text_repel(aes(label = Sample), size = 3, max.overlaps = 20, show.legend = FALSE) +
  scale_color_manual(values = c("Bulbar" = "#7CAE00", "Spinal" = "#C77CFF"))+
  labs(
    title = "PCA of Methylation Profiles (ALS: Bulbar vs Spinal) ",
    x = sprintf("PC1 (%.1f%% variance)", 100 * summary(pca)$importance[2, 1]),
    y = sprintf("PC2 (%.1f%% variance)", 100 * summary(pca)$importance[2, 2])
  ) +
  theme_minimal(base_size = 15) +
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5))

ggsave(file.path(config$output_dir, "QC_PCA.png"),p_pca,width = 8,height = 6,dpi = 300)

# 4.3.3. Coverage Uniformity
cov_data <- getCoverage(bsseq_merged, type = "Cov")

# Convert to data frame with chromosome info
cov_df <- data.frame(
  Chr = as.character(seqnames(bsseq_merged)),
  Coverage = rowMeans(cov_data, na.rm = TRUE),
  stringsAsFactors = FALSE
)

# Calculate mean coverage per chromosome
cov_by_chr <- cov_df %>%
  group_by(Chr) %>%
  summarise(Mean_Coverage = mean(Coverage, na.rm = TRUE)) %>%
  filter(Chr %in% paste0("chr", 1:22))  # Filter to standard chromosomes

# Plot
p_cov <- ggplot(cov_by_chr, aes(x = Chr, y = Mean_Coverage)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(
    title = "Mean Coverage by Chromosome",
    x = "Chromosome",
    y = "Mean Coverage"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(config$output_dir, "QC_coverage_by_chr.png"), p_cov, width = 8, height = 6)

# 4.3.4. Clustering samples
library(pheatmap)
cor_matrix <- cor(meth, use = "pairwise.complete.obs")
dist_matrix <- as.dist(1 - cor_matrix)
hc <- hclust(dist_matrix, method = "ward.D2")
p_feature<-plot(hc, main = "Sample Clustering (Ward.D2 + Correlation)", xlab = "")
ggsave(file.path(config$output_dir, "QC_methylation_by_feature.png"), p_feature, width = 8, height = 6)


# 4.3.4. Global Methylation by Condition

# Calculate average methylation per sample
meth_avg <- colMeans(meth, na.rm = TRUE)

global_meth <- data.frame(
  sample = colnames(meth),
  methylation = meth_avg
) %>%
  left_join(sample_meta, by = c("sample" = "name"))

# Plot
p_global <- ggplot(global_meth, aes(x = group, y = methylation, fill = group)) +
  #geom_violin(trim = FALSE, color = NA, alpha = 0.7) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.9) +
  scale_fill_manual(values = c("Control" = "#FFFACD", "ALS" = "#ADD8E6")) +
  stat_compare_means(method = "t.test",label = "p.format",label.x.npc = "center",label.y.npc = "top", size = 5) +
  theme_minimal() +
  labs(
    title = "Global methylation levels by group",
    y = "Average methylation (%)",
    x = ""
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme(legend.position = "none")

ggsave(file.path(config$output_dir, "Global_Methylation_Comparison.png"),
       p_global, width = 6, height = 5, dpi = 300)

# ==========================================
# 5. CALLING AND ANNOTATION OF DMLs AND DMRs
# ==========================================
# 5.1 Call DMLs
dml_files <- list.files(pattern = "dmlTest_.*_filtered\\.rds$")
dml_all <- lapply(dml_files, readRDS) %>% bind_rows()
dmls <- callDML(dml_all, delta=0.1 , p.threshold = 0.05)

# 5.2 Call DMRs
dmrs <- do.call(callDMR, c(list(dml_all), params$dmr_params))

# 5.3 Annotate DMLs and DMRs
annotations <- build_annotations("hg38", c("hg38_basicgenes", "hg38_genes_intergenic"))
annot_cpg <- build_annotations("hg38", annotations = c("hg38_cpgs"))

# DML annotation
dml_gr <- GRanges(seqnames = dmls$chr, ranges = IRanges(dmls$pos, dmls$pos + 1))
annot_dml <- annotate_regions(dml_gr, annotations, ignore.strand = TRUE) %>%
  as.data.frame() %>%
  dplyr::rename(chr = seqnames, pos = start) %>%
  group_by(chr, pos) %>%
  summarise(
    annot.symbol = paste(unique(annot.symbol), collapse = ";"),
    annot.type = paste(unique(annot.type), collapse = ";"),
    annot.gene_id = paste(unique(annot.gene_id), collapse = ";"),
    .groups = "drop"
  )

dmls_annot <- dmls %>%
  left_join(annot_dml, by = c("chr", "pos")) %>%
  dplyr::select(chr, pos, mu1, mu2, diff, pval, fdr, annot.symbol, annot.type, annot.gene_id) %>%
  arrange(fdr, desc(abs(diff)))

write.table(
  dmls_annot,
  file = file.path(config$output_dir, "DMLs_annotated.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# DMR annotation
dmr_gr <- GRanges(seqnames = dmrs$chr, ranges = IRanges(dmrs$start, dmrs$end))
annot_dmr <- annotate_regions(dmr_gr, annotations, ignore.strand = TRUE) %>%
  as.data.frame() %>%
  dplyr::rename(chr = seqnames, start = start, end = end) %>%
  group_by(chr, start, end) %>%
  summarise(
    annot.symbol = paste(unique(annot.symbol), collapse = ";"),
    annot.type = paste(unique(annot.type), collapse = ";"),
    .groups = "drop"
  )

dmrs_annot <- dmrs %>%
  left_join(annot_dmr, by = c("chr", "start", "end"))

write.table(
  dmrs_annot,
  file = file.path(config$output_dir, "DMRs_annotated.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# ======================
# 6. VISUALIZATION
# ======================

# 6.1 Create summary tables
dml_summary <- summarize_regions(dmls_annot, "diff", "DML")
dmr_summary <- summarize_regions(dmrs_annot, "diff.Methy", "DMR")

cat("Differentially Methylated Loci (DML) Summary:\n")
print(dml_summary, row.names = FALSE)

cat("\nDifferentially Methylated Regions (DMR) Summary:\n")
print(dmr_summary, row.names = FALSE)

# 6.2 DML visualization
dml_annot_df <- as.data.frame(annotate_regions(dml_gr, annotations, ignore.strand = TRUE))
dml_cpg_annot_df <- as.data.frame(annotate_regions(dml_gr, annot_cpg, ignore.strand = TRUE))

p1 <- plot_annotation_pie(dml_annot_df, "DMLs - Gene Annotation")
p2 <- plot_annotation_pie(dml_cpg_annot_df, "DMLs - CpG Annotation")
p3 <- plot_fold_enrichment(dml_gr, annotations, "DMLs - Enrichment (Genes)")
p4 <- plot_fold_enrichment(dml_gr, annot_cpg, "DMLs - Enrichment (CpG Features)")

p_dml <- (p1 | p2) / (p3 | p4)
ggsave(
  file.path(config$output_dir, "DML_annotation_summary.png"),
  p_dml, width = 12, height = 10, dpi = 300
)

# 6.3 DMR visualization
dmr_annot_df <- as.data.frame(annotate_regions(dmr_gr, annotations, ignore.strand = TRUE))
dmr_cpg_annot_df <- as.data.frame(annotate_regions(dmr_gr, annot_cpg, ignore.strand = TRUE))

p5 <- plot_annotation_pie(dmr_annot_df, "DMRs - Gene Annotation")
p6 <- plot_annotation_pie(dmr_cpg_annot_df, "DMRs - CpG Annotation")
p7 <- plot_fold_enrichment(dmr_gr, annotations, "DMRs - Enrichment (Genes)")
p8 <- plot_fold_enrichment(dmr_gr, annot_cpg, "DMRs - Enrichment (CpG Features)")

p_dmr <- (p5 | p6) / (p7 | p8)
ggsave(
  file.path(config$output_dir, "DMR_annotation_summary.png"),
  p_dmr, width = 12, height = 10, dpi = 300
)

# ======================
# 7. VOLCANO PLOT AND HEATMAP OF DMLs
# ======================
# vOLCANO PLOT

volcano_plot <- ggplot(dmls_annot, aes(x = diff, y = -log10(pval))) +
  geom_point(aes(color = fdr < 0.05), alpha = 0.5) +
  scale_color_manual(values = c("grey", "red"), name = "FDR < 0.05") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = "blue") +
  labs(
    title = "Volcano Plot of DMLs",
    x = "Difference in Methylation (ALS - Control)",
    y = "-log10(p-value)"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave(
  file.path(config$output_dir, "Volcano_DMLs.png"),
  volcano_plot, width = 8, height = 6, dpi = 300
)
# HEATMAP OF TOP DMLs
# Select top 100 DMLs by FDR
top_dmls <- dmls_annot %>%
  arrange(fdr, desc(abs(diff))) %>%
  head(50)

# 2. Create GRanges object from DML positions
top_gr <- GRanges(seqnames = top_dmls$chr,
                  ranges = IRanges(top_dmls$pos, width = 1))

# 3. Subset methylation matrix for these positions
bsseq_subset <- subsetByOverlaps(bsseq_merged, top_gr)
meth_top <- getMeth(bsseq_subset, type = "raw")

# 4. Match rows between methylation matrix and top DMLs
row_ids <- paste0(seqnames(granges(bsseq_subset)), ":", start(granges(bsseq_subset)))
target_ids <- paste0(top_dmls$chr, ":", top_dmls$pos)
meth_top <- meth_top[match(target_ids, row_ids, nomatch = NA), ]

# 5. Scale methylation values and reorder samples by group
sample_order <- sample_meta %>%
  arrange(group) %>%
  pull(name)

meth_scaled <- t(scale(t(meth_top), center = TRUE, scale = TRUE))
meth_scaled <- meth_scaled[complete.cases(meth_scaled), ]
meth_scaled <- meth_scaled[, sample_order]

# 6. Clean gene names or replace with coordinates if missing
gene_names <- top_dmls$annot.symbol

# If multiple names separated by ';', take the first valid one
gene_names <- sapply(strsplit(as.character(gene_names), ";"), function(x) {
  first_valid <- x[!is.na(x) & x != "NA"]
  if (length(first_valid) == 0) return(NA) else return(first_valid[1])
})

# Replace missing names with genomic coordinates
missing_genes <- is.na(gene_names) | gene_names == ""
gene_names[missing_genes] <- paste0(top_dmls$chr[missing_genes], ":", top_dmls$pos[missing_genes])

# Ensure unique row names
gene_names <- make.unique(gene_names)

# Set gene names as row names
rownames(meth_scaled) <- gene_names

# 7. Define column annotation (sample groups)
annotation_col <- data.frame(Group = sample_meta$group)
rownames(annotation_col) <- sample_meta$name
annotation_col <- annotation_col[sample_order, , drop = FALSE]

# 8. Generate heatmap
pheatmap(meth_scaled,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         annotation_col = annotation_col,
         main = "Top 50 DMLs (scaled methylation)",
         filename = file.path(config$output_dir, "Heatmap_DMLs.png"),
         width = 10, height = 10,
         show_rownames = TRUE,
         fontsize_row = 8)
      
# ======================
# 9. SESSION INFO
# ======================

writeLines(capture.output(sessionInfo()), file.path(config$output_dir, "session_info_BvsS.txt"))

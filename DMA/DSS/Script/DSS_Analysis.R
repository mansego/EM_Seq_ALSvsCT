#!/usr/bin/env Rscript

# ======================
# 1. INITIAL SETUP
# ======================

# Load required packages with installation if needed
required_packages <- c(
  "DSS", "bsseq", "annotatr", "pheatmap", "GenomicRanges",
  "data.table", "parallel", "magrittr", "dplyr", "ggplot2", 
  "stringr", "qvalue", "patchwork"
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
  work_dir = "/mnt/mydisk/EM_Seq_ALSvsCT/DMA/DSS/RDS/Filter_SNP_cov8/",
  base_path = "/mnt/mydisk/EM_Seq_ALSvsCT/PrePro/bismark/methylation_call/methylation_coverage/",
  snp_file = "/mnt/mydisk/EM_Seq_ALSvsCT/references/snp151Common.txt.gz",
  output_dir = "/mnt/mydisk/EM_Seq_ALSvsCT/DMA/DSS/Results"
)

# Create directories if they don't exist
dir.create(config$work_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(config$output_dir, showWarnings = FALSE, recursive = TRUE)

# Sample metadata
sample_meta <- data.frame(
  id = c(42, 44, 45, 47, 48, 50, 52, 53, 55, 56, 57, 58, 59, 61, 64, 65),
  name = c("ALS42", "CT44", "ALS45", "ALS47", "CT48", "ALS50", "CT52", "ALS53",
           "CT55", "ALS56", "CT57", "ALS58", "CT59", "CT61", "CT64", "ALS65"),
  stringsAsFactors = FALSE
) %>%
  mutate(
    group = ifelse(grepl("^ALS", name), "ALS", "Control"),
    group = factor(group, levels = c("Control", "ALS")),
    cov_file = file.path(config$base_path, paste0("BACW_", id, "_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz"))
  )

# Analysis parameters
params <- list(
  chromosomes = paste0("chr", 1:22),
  maf_cutoff = 0.10,
  min_coverage = 8,
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

#' Process methylation data for a single chromosome
process_chromosome <- function(chr, sample_meta, snp_gr, min_coverage = 8, min_samples = 4) {
  message("Processing ", chr)
  current_chr <- chr

  # Read and process coverage data
  cov_list <- lapply(sample_meta$cov_file, function(f) {
  dat <- fread(f, header = FALSE,
               col.names = c("chr", "start", "end", "meth_pct", "meth", "unmeth"))
  dat_chr <- dat[chr == current_chr, .(chr, pos = start, X = meth, N = meth + unmeth)]
  if (nrow(dat_chr) == 0) return(NULL)
  return(dat_chr)
})
  
  if (all(sapply(cov_list, nrow) == 0)) {
    message("No data for ", chr, ", skipping...")
    return(NULL)
  }
  
  # Create BSseq object and filter SNPs
  bsseq_obj <- makeBSseqData(cov_list, sampleNames = sample_meta$name) %>%
    keepSeqlevels(value = chr, pruning.mode = "coarse")
   
  # Filter SNPs
  cpg_gr <- granges(bsseq_obj)
  seqlevelsStyle(cpg_gr) <- "NCBI"
  
  # Find overlaps with SNPs
  overlap <- findOverlaps(cpg_gr, snp_gr)
  
  # Filter out overlapping CpGs
    bsseq_filtered <- bsseq_obj[-queryHits(overlap), ]

  # Filter by coverage
  cov <- getCoverage(bsseq_filtered, type = "Cov")
  keep <- rowSums(cov[, sample_meta$group == "ALS"] >= min_coverage) >= min_samples &
          rowSums(cov[, sample_meta$group == "Control"] >= min_coverage) >= min_samples
  bsseq_filtered <- bsseq_filtered[keep, ]
  
  # Perform DML test
  dml_test <- DMLtest(
    bsseq_filtered,
    group1 = which(sample_meta$group == "ALS"),
    group2 = which(sample_meta$group == "Control"),
    smoothing = TRUE
  )
  
  # Save results
  saveRDS(bsseq_filtered, file.path(config$work_dir, paste0("bsseq_", chr, "_filtered.rds")))
  saveRDS(dml_test, file.path(config$work_dir, paste0("dmlTest_", chr, "_filtered.rds")))
  
  # Return summary
  list(
    chromosome = chr,
    CpGs_original = length(bsseq_obj),
    CpGs_removed = length(unique(queryHits(overlap))),
    CpGs_remaining = length(bsseq_filtered)
  )
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

# 4.2 Process each chromosome
results <- lapply(params$chromosomes, function(chr) {
  process_chromosome(chr, sample_meta, snp_gr, params$min_coverage, params$min_samples)
}) %>% bind_rows()

write.csv(results, "CpG_filtered_by_SNP_cov8.csv", row.names = FALSE)
print(results)

# 4.3 Combine DML results
dml_files <- list.files(pattern = "dmlTest_.*_filtered\\.rds$")
dml_all <- lapply(dml_files, readRDS) %>% bind_rows()
dmls <- callDML(dml_all, delta=0.1 , p.threshold = 0.05)

# 4.4 Call DMRs
dmrs <- do.call(callDMR, c(list(dml_all), params$dmr_params))

# 4.5 Annotate DMLs and DMRs
annotations <- build_annotations("hg38", c("hg38_basicgenes", "hg38_genes_intergenic"))
annot_cpg <- build_annotations("hg38", annotations = c("hg38_cpgs"))

# DML annotation
dml_gr <- GRanges(seqnames = dmls$chr, ranges = IRanges(dmls$pos, dmls$pos + 1))
annot_dml <- annotate_regions(dml_gr, annotations, ignore.strand = TRUE) %>%
  as.data.frame() %>%
  rename(chr = seqnames, pos = start) %>%
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
  rename(chr = seqnames, start = start, end = end) %>%
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
# 5. VISUALIZATION
# ======================

# 5.1 Create summary tables
dml_summary <- summarize_regions(dmls_annot, "diff", "DML")
dmr_summary <- summarize_regions(dmrs_annot, "diff.Methy", "DMR")

cat("Differentially Methylated Loci (DML) Summary:\n")
print(dml_summary, row.names = FALSE)

cat("\nDifferentially Methylated Regions (DMR) Summary:\n")
print(dmr_summary, row.names = FALSE)

# 5.2 DML visualization
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

# 5.3 DMR visualization
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
# 6. SESSION INFO
# ======================

writeLines(capture.output(sessionInfo()), file.path(config$output_dir, "session_info.txt"))

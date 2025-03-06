library(tidyverse)
library(DESeq2)
library(rstudioapi)
library(ggplot2)
library(dplyr)
library(qqman)

#---------------------------------Loading Data----------------------------------
setwd(dirname(getActiveDocumentContext()$path))
sample_info <- read_tsv("../addresses/sample_metadata.tsv")
counts_with_locs <- read_tsv("../../data/human_binned_sequences/ctMatrices/geneCounts/150_180_gene_counts.tsv")
# bin_sizes <- "150_180"

#----------------------------Processing counts----------------------------------

# There shouldn't ever be a mixture of bins here
# bin_sizes <- unique(gsub("A\\d+_", "", colnames(counts_with_locs)))
# unique(gsub(".*_(\\d+_\\d+)$", "\\1", colnames(counts_with_locs), perl=TRUE))


# Fixing columns, removing extra information, gene_ids to row name (gene names has duplicates)
count_mtx <- counts_with_locs %>%
  select(-chr, -start, -end, -gene_name) %>%
  column_to_rownames("gene_id") %>%
  rename_with(~ gsub("_.*", "", .))
  
# Save gene metadata for later
gene_metadata <- counts_with_locs %>%
  select(gene_id, gene_name, chr, start, end)

#----------------------------Processing Metadata----------------------------------

sample_info <- sample_info %>%
  mutate(TGen_ID = as.character(TGen_ID)) %>%
  column_to_rownames("TGen_ID")
  
# Ensure columns in count_mat match row names in sample_info
all(colnames(count_mtx) %in% rownames(sample_info))  # Should return TRUE
all(rownames(sample_info) %in% colnames(count_mtx))  # Should return TRUE

# Reordering metadata to match matrix columns
sample_info <- sample_info[colnames(count_mtx), ]

#--------------------------Check overdispersion-------------------
# dds <- DESeqDataSetFromMatrix(
#   countData = round(as.matrix(count_mtx)),   # Ensure integer counts
#   colData   = sample_info,  # Sample metadata
#   design    = ~ Timepoint + Tube + Input  # Adjust as needed
# )
# 
# dds <- DESeq(dds)  # Estimate dispersions
# 
# png("../../results/BinAnalysis_R/dispersionPlot_150_180.png", width = 1600, height = 1200, res = 300)
# plotDispEsts(dds)  # Generate the plot
# dev.off()
# 
# poisson_model <- estimateSizeFactors(dds)
# poisson_variance <- rowMeans(counts(poisson_model, normalized=TRUE))
# observed_variance <- rowVars(counts(dds, normalized=TRUE))
# 
# dispersion_ratio <- observed_variance / poisson_variance
# summary(dispersion_ratio)

#-------------------- Tp2 vs all----------------------------------------

# Converting all timepoints that are not 2 to other for comparisons
sample_info$TP2_vs_All <- ifelse(sample_info$Timepoint == 2, "TP2", "Other")
sample_info$TP2_vs_All <- factor(sample_info$TP2_vs_All, levels = c("Other", "TP2"))

dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(count_mtx)),   # Ensure integer counts
  colData   = sample_info,  # Updated metadata
  design    = ~ TP2_vs_All + Tube + Input  # Adjust as needed
)

# Run deseq2
dds <- DESeq(dds)

# Extract Results for TP2 vs. All Other Timepoints
res <- results(dds, contrast = c("TP2_vs_All", "TP2", "Other"), alpha = 0.15)


# Selecting by padj
sig_genes_padj_0.15 <- res %>%
  as.data.frame() %>%
  rownames_to_column("gene_id") %>%
  filter(padj < 0.15) %>%
  arrange(padj)

# Selecting by pvalue
sig_genes_pval_0.05 <- res %>%
  as.data.frame() %>%
  rownames_to_column("gene_id") %>%
  filter(pvalue < 0.05) %>%
  arrange(pvalue)


# Merging location 
sig_genes_with_loc <- sig_genes_pval_0.05 %>%
  inner_join(gene_metadata, by="gene_id")

# Manhattan Plot?
ggplot(sig_genes_with_loc, aes(x = start, y = -log10(padj), color = chr)) +
  geom_point(alpha = 0.8) +
  labs(x = "Genomic Position", y = "-log10(Adjusted P-value)", title = "Manhattan Plot: TP2 vs. All") +
  theme_minimal()


#----------------------------- Manhattan plot----------------------------------

manhattan_df <- sig_genes_with_loc %>%
  dplyr::rename(SNP = gene_id, CHR = chr, BP = start, P = pvalue) %>%
  mutate(
    CHR = case_when(
      CHR == "X"  ~ 23,  # Convert X to 23
      CHR == "Y"  ~ 24,  # Convert Y to 24
      CHR == "MT" ~ 25,  # Convert MT to 25
      TRUE ~ as.numeric(CHR)  # Keep numeric chromosomes as they are
    )
  ) %>%
  drop_na()  # Remove any missing values

manhattan(manhattan_df, col = c("blue4", "orange3"), genomewideline = -log10(0.05), suggestiveline = -log10(0.1), main = "Manhattan Plot: TP2 vs. All")


#-------------------------- manhattan all results ---------------------------
res_df <- as.data.frame(res) %>%
  tibble::rownames_to_column(var = "gene_id")  # Ensure gene_id is a column

res_df <- res_df %>%
  inner_join(gene_metadata, by="gene_id")

manhattan_df_all <- res_df %>%
  dplyr::rename(SNP = gene_id, CHR = chr, BP = start, P = pvalue) %>%
  mutate(
    CHR = case_when(
      CHR == "X"  ~ 23,  # Convert X to 23
      CHR == "Y"  ~ 24,  # Convert Y to 24
      CHR == "MT" ~ 25,  # Convert MT to 25
      TRUE ~ as.numeric(CHR)  # Keep numeric chromosomes as they are
    )
  ) %>%
  drop_na()  # Remove any missing values

manhattan(manhattan_df_all, col = c("blue4", "orange3"), genomewideline = -log10(0.05), suggestiveline = -log10(0.1), main = "Manhattan Plot: TP2 vs. All")




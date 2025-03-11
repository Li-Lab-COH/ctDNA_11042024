library(tidyverse)
library(DESeq2)
library(rstudioapi)
library(ggplot2)
library(dplyr)
library(qqman)

#---------------------------------Loading Data----------------------------------
setwd(dirname(getActiveDocumentContext()$path))
sample_info <- read_tsv("../../addresses/sample_metadata.tsv")
counts_with_locs <- read_tsv("../../../data/human_binned_sequences/ctMatrices/nuclCounts/125_155_nucleosome_counts.tsv")

#----------------------------Processing counts----------------------------------

# Fixing columns, removing extra information, nucleosome_id to row names
count_mtx <- counts_with_locs %>%
  select(-chr, -start, -end) %>%  # Remove location columns for now
  column_to_rownames("nuc_id") %>%
  rename_with(~ gsub("_.*", "", .))  # Rename sample columns

# Save nucleosome metadata for later
nuc_metadata <- counts_with_locs %>%
  select(nuc_id, chr, start, end)

#Filderting low quality reads
count_mtx <- count_mtx[rowSums(count_mtx > 0) >= 3 & rowMeans(count_mtx) >= 1 & apply(count_mtx, 1, max) >= 5, ]

# Ensure nuc_metadata only contains nucleosome sites that remain in count_mtx
nuc_metadata <- nuc_metadata %>%
  filter(nuc_id %in% rownames(count_mtx))

rm(counts_with_locs)

#----------------------------Processing Metadata----------------------------------

sample_info <- sample_info %>%
  mutate(TGen_ID = as.character(TGen_ID)) %>%
  column_to_rownames("TGen_ID") %>%
  mutate(
    Timepoint = factor(Timepoint),
    Tube = factor(Tube),
    Input = factor(Input)
  )

# Ensure sample order matches count matrix
sample_info <- sample_info[colnames(count_mtx), ]

#--------------------------Check overdispersion-------------------
dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(count_mtx)),   # Ensure integer counts
  colData   = sample_info,  # Sample metadata
  design    = ~ Timepoint + Tube + Input  # Adjust as needed
)

dds <- DESeq(dds)  # Estimate dispersions

png("../../../results/BinAnalysis_R/dispersionPlot_Nucleosome.png", width = 1600, height = 1200, res = 300)
plotDispEsts(dds)  # Generate the plot
dev.off()

poisson_model <- estimateSizeFactors(dds)
poisson_variance <- rowMeans(counts(poisson_model, normalized=TRUE))
observed_variance <- rowVars(counts(dds, normalized=TRUE))

dispersion_ratio <- observed_variance / poisson_variance
summary(dispersion_ratio)


# Extract dispersion values from DESeq2
dispersion_values <- dispersions(dds)

# Convert to dataframe for visualization
dispersion_df <- data.frame(dispersion = dispersion_values)

# Plot histogram with log-scaled y-axis
ggplot(dispersion_df, aes(x = dispersion)) +
  geom_histogram(bins = 50, fill = "blue", color = "black", alpha = 0.7) +
  scale_y_log10() +  # Log scale for frequency
  labs(
    title = "Histogram of Dispersion Values (DESeq2)",
    x = "Dispersion",
    y = "Frequency (log scale)"
  ) +
  theme_minimal()

# Filtering low variable features
low_dispersion_threshold <- quantile(dispersion_values, 0.10)

# Remove nucleosome sites with extremely low variance
filtered_dds <- dds[dispersions(dds) > low_dispersion_threshold, ]

plotDispEsts(filtered_dds)
#-------------------- Tp2 vs all----------------------------------------

# Convert all timepoints that are not 2 to "Other"
sample_info$TP2_vs_All <- ifelse(sample_info$Timepoint == 2 , "TP2", "Other")
sample_info$TP2_vs_All <- factor(sample_info$TP2_vs_All, levels = c("Other", "TP2"))

dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(count_mtx)),   # Ensure integer counts
  colData   = sample_info,  # Updated metadata
  design    = ~ TP2_vs_All + Tube + Input  # Adjust as needed
)

# Run DESeq2
dds <- DESeq(dds)

# Extract Results for TP2 vs. All Other Timepoints
res <- results(dds, contrast = c("TP2_vs_All", "TP2", "Other"), alpha = 0.15)


# Selecting by padj
sig_nuc_padj_0.15 <- res %>%
  as.data.frame() %>%
  rownames_to_column("nuc_id") %>%
  filter(padj < 0.15) %>%
  arrange(padj)

# Selecting by pvalue
sig_nuc_pval_0.05 <- res %>%
  as.data.frame() %>%
  rownames_to_column("nuc_id") %>%
  filter(pvalue < 0.05) %>%
  arrange(pvalue)

# Merging location data
sig_nuc_with_loc <- sig_nuc_pval_0.05 %>%
  inner_join(nuc_metadata, by="nuc_id")

#----------------------------- Manhattan plot----------------------------------

# Prepare data for Manhattan plot
manhattan_df <- sig_nuc_with_loc %>%
  dplyr::rename(SNP = nuc_id, CHR = chr, BP = start, P = pvalue) %>%
  mutate(
    CHR = case_when(
      CHR == "X"  ~ 23,  # Convert X to 23
      CHR == "Y"  ~ 24,  # Convert Y to 24
      CHR == "MT" ~ 25,  # Convert MT to 25
      TRUE ~ as.numeric(CHR)  # Keep numeric chromosomes as they are
    )
  ) %>%
  drop_na()  # Remove any missing values

# Plot the Manhattan plot
manhattan(
  manhattan_df, 
  col = c("blue4", "orange3"), 
  genomewideline = -log10(0.05), 
  suggestiveline = -log10(0.1), 
  main = "Manhattan Plot: TP2 vs. All (Nucleosome Data)"
)

#-------------------------- Highlight significant nucleosome sites ---------------------------
res_df <- as.data.frame(res) %>%
  tibble::rownames_to_column(var = "nuc_id")  # Ensure nuc_id is a column

res_df <- res_df %>%
  inner_join(nuc_metadata, by="nuc_id")

manhattan_df_all <- res_df %>%
  dplyr::rename(SNP = nuc_id, CHR = chr, BP = start, P = pvalue) %>%
  mutate(
    CHR = case_when(
      CHR == "X"  ~ 23,  # Convert X to 23
      CHR == "Y"  ~ 24,  # Convert Y to 24
      CHR == "MT" ~ 25,  # Convert MT to 25
      TRUE ~ as.numeric(CHR)  # Keep numeric chromosomes as they are
    )
  ) %>%
  drop_na()  # Remove any missing values

# Identify significant nucleosome sites
sig_nuc_padj_0.15 <- res_df %>%
  filter(padj < 0.25) %>%
  arrange(padj)

highlight_nucs <- sig_nuc_padj_0.15 %>%
  select(nuc_id)

manhattan_df_highlight <- manhattan_df_all %>%
  inner_join(highlight_nucs, by = c("SNP" = "nuc_id"))

# Manhattan plot with highlighted nucleosome sites
png("../../../results/BinAnalysis_R/Nucleosome_Manhattan_TP2_vs_All.png", width = 1200, height = 800, res = 150)
manhattan(
  manhattan_df_all, 
  col = c("blue4", "orange3"), 
  genomewideline = -log10(0.05), 
  suggestiveline = -log10(0.1), 
  highlight = manhattan_df_highlight$SNP, 
  main = "Nucleosome - TP2 vs. All - padj: 0.25"
)
dev.off()


#-------------------------- some more analysis -----------------------------

manhattan_df_all <- manhattan_df_all %>% 
  mutate(difference = abs(BP - end))

summary(manhattan_df_all$difference)
hist(manhattan_df_all$difference)

# Calculate the 90th percentile threshold
threshold <- quantile(manhattan_df_all$difference, 0.99)

# Subset the data frame to keep rows below or equal to the threshold
manhattan_df_all_filtered <- manhattan_df_all[manhattan_df_all$difference <= 120, ]

hist(manhattan_df_all_filtered$difference)

summary(manhattan_df_all_filtered$difference)


boxplot(manhattan_df_all_filtered$difference, 
        main = "Boxplot of Difference",
        ylab = "Difference",
        col = "lightblue",
        border = "darkblue")




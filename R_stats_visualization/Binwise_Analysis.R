###############################################################################
# Load libraries
###############################################################################
library(dplyr)
library(ggplot2)
library(purrr)      # for functional programming (map, etc.)
library(tidyr)
library(forcats)    # for factor reordering, if needed
library(reshape2)   # if you want melt/cast for heatmaps
library(pheatmap)   # or ComplexHeatmap, whichever you prefer
library(stringr)

###############################################################################
# Read in Data
###############################################################################
# Replace 'combined_data.csv' with your actual data file
combined_data <- read.csv("../../data/combined_data_hist.csv", header = TRUE, stringsAsFactors = FALSE)

###############################################################################
# Data Preparation
###############################################################################
# 1) Define your bins. Example below uses 50-bp bin widths from 0 to 350 bp
#    Adjust to your desired bin scheme. 
bin_breaks <- seq(0, 350, by = 50)  # 0-50, 51-100, 101-150, etc.
bin_labels <- paste(head(bin_breaks, -1) + 1, bin_breaks[-1], sep = "-")

# Create a new column 'SizeBin'
combined_data <- combined_data %>%
  mutate(SizeBin = cut(InsertSize,
                       breaks = bin_breaks,
                       labels = bin_labels,
                       include.lowest = TRUE,
                       right = TRUE))

# 2) Calculate total frequency per sample (for CPM calculation)
#    We'll create a total read count (sum of all frequencies for each sample).
total_counts <- combined_data %>%
  group_by(TGen_ID) %>%
  summarise(TotalFrequency = sum(Frequency, na.rm = TRUE))

# 3) Bin frequencies
binned_data <- combined_data %>%
  group_by(TGen_ID, Patient_ID, Timepoint, Tube, Input, SizeBin) %>%
  summarise(BinFrequency = sum(Frequency, na.rm = TRUE)) %>%
  ungroup() %>%
  left_join(total_counts, by = "TGen_ID") %>%
  mutate(CPM = (BinFrequency / TotalFrequency) * 1e6)

###############################################################################
# Compute Z-scores
###############################################################################
# For a heatmap, you typically want each row to be a bin and each column a sample,
# or vice versa. Then you compute the Z-score across either samples or bins.
# We'll compute Z-scores across samples (for each bin separately).

# Reshape data into wide format for Z-score calculation:
binned_wide <- binned_data %>%
  select(TGen_ID, Tube, Input, Timepoint, SizeBin, BinFrequency) %>%
  pivot_wider(names_from = TGen_ID, values_from = BinFrequency, values_fill = 0)

# We'll keep track of the SizeBin and other metadata in separate objects
sizebins <- binned_wide$SizeBin
# The matrix of frequencies
freq_matrix <- as.matrix(binned_wide[, -1])  # remove the SizeBin column

# Function to compute row-wise z-scores
row_zscore <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

z_matrix <- t(apply(freq_matrix, 1, row_zscore))
colnames(z_matrix) <- colnames(freq_matrix)
rownames(z_matrix) <- sizebins

# If you want to compute z-scores by columns (i.e. by sample across bins),
# you can do:  z_matrix <- apply(freq_matrix, 2, row_zscore)

###############################################################################
# Directory Structure
###############################################################################
# You mentioned categories like:
# 1) EDTA (10 and 2.5)
# 2) Streck (10 and 2.5)
# 3) EDTA (10)
# 4) EDTA (2.5)
# 5) Streck (10)
# 6) Streck (2.5)
# 7) EDTA, Streck (2.5)
# 8) EDTA, Streck (10)
#
# We'll create a vector of these categories and then filter binned_data accordingly.

categories <- list(
  "EDTA_2.5_10"       = quote(Tube == "EDTA" & Input %in% c(2.5, 10)),
  "Streck_2.5_10"     = quote(Tube == "Streck" & Input %in% c(2.5, 10)),
  "EDTA_10"           = quote(Tube == "EDTA" & Input == 10),
  "EDTA_2.5"          = quote(Tube == "EDTA" & Input == 2.5),
  "Streck_10"         = quote(Tube == "Streck" & Input == 10),
  "Streck_2.5"        = quote(Tube == "Streck" & Input == 2.5),
  "EDTA_Streck_2.5"   = quote(Tube %in% c("EDTA","Streck") & Input == 2.5),
  "EDTA_Streck_10"    = quote(Tube %in% c("EDTA","Streck") & Input == 10)
)

# Create output directory if needed
base_dir <- "Fragment_Size_Analysis"
if(!dir.exists(base_dir)) dir.create(base_dir)

###############################################################################
# Helper functions to plot
###############################################################################

# 1) Z-score heatmap
plot_zscore_heatmap <- function(zscore_mat, output_file, main_title="Z-score Heatmap") {
  pheatmap::pheatmap(
    mat               = zscore_mat,
    cluster_rows      = TRUE,
    cluster_cols      = TRUE,
    scale             = "none",  # we already computed z-scores
    show_rownames     = TRUE,
    show_colnames     = TRUE,
    main              = main_title,
    filename          = output_file,
    width             = 8,
    height            = 6
  )
}

# 2) Barplot of raw frequencies
plot_raw_bar <- function(data_in, output_file, main_title="Raw Frequency Barplot") {
  p <- ggplot(data_in, aes(x = SizeBin, y = BinFrequency, fill = as.factor(Timepoint))) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_bw() +
    labs(title = main_title, x = "Size Bin", y = "Raw Frequency") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(output_file, plot = p, width = 8, height = 6)
}

# 3) Barplot of CPM
plot_cpm_bar <- function(data_in, output_file, main_title="CPM Barplot") {
  p <- ggplot(data_in, aes(x = SizeBin, y = CPM, fill = as.factor(Timepoint))) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_bw() +
    labs(title = main_title, x = "Size Bin", y = "CPM") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(output_file, plot = p, width = 8, height = 6)
}

###############################################################################
# Main Loop Over Categories
###############################################################################
for(cat_name in names(categories)) {
  
  # Make subdirectory for the category
  cat_dir <- file.path(base_dir, cat_name)
  if(!dir.exists(cat_dir)) dir.create(cat_dir)
  
  # Filter the data for the category
  # Using the '!!!' pattern plus rlang::parse_expr if using a string expression
  # We stored a quoted expression, so we can do:
  condition_expr <- categories[[cat_name]]
  
  category_data <- binned_data %>%
    filter(!!condition_expr)
  
  # -- 1) Z-score heatmap for this category --
  # Here we must subset the z_matrix or re-construct a subset matrix
  # if we want a "category-specific" heatmap (only for the samples in the category).
  # We'll do the simpler approach of re-constructing a frequency matrix from
  # just the subset category_data, then compute z-scores.
  
  # Pivot wide
  cat_wide <- category_data %>%
    select(TGen_ID, SizeBin, BinFrequency) %>%
    pivot_wider(names_from = TGen_ID, values_from = BinFrequency, values_fill = 0)
  
  # If there is no data or only 1 sample, skip plotting to avoid errors
  if(ncol(cat_wide) < 3) {
    message(paste("Skipping category", cat_name, "because it has too few samples."))
    next
  }
  
  cat_sizebins <- cat_wide$SizeBin
  cat_freq_mat <- as.matrix(cat_wide[, -1])
  rownames(cat_freq_mat) <- cat_sizebins
  
  # compute row-wise z-scores
  cat_z_matrix <- t(apply(cat_freq_mat, 1, row_zscore))
  colnames(cat_z_matrix) <- colnames(cat_freq_mat)
  rownames(cat_z_matrix) <- cat_sizebins
  
  # Plot and save z-score heatmap
  zscore_file <- file.path(cat_dir, paste0("ZScoreHeatmap_", cat_name, ".pdf"))
  plot_zscore_heatmap(cat_z_matrix, zscore_file,
                      main_title = paste("Z-score Heatmap -", cat_name))
  
  # -- 2) Barplot of raw frequencies for each sample/timepoint --
  # We'll produce one combined barplot showing the sum of frequencies by SizeBin
  # across the different TGen_ID or Timepoints. The easiest is a stacked or
  # dodged bar plot. We'll do a single plot for the entire category.
  
  # If you'd rather have separate plots per timepoint, loop again. For now, one plot:
  raw_file <- file.path(cat_dir, paste0("RawBarplot_", cat_name, ".pdf"))
  plot_raw_bar(category_data, raw_file,
               main_title = paste("Raw Frequency -", cat_name))
  
  # -- 3) Barplot of CPM for each sample/timepoint --
  cpm_file <- file.path(cat_dir, paste0("CPMBarplot_", cat_name, ".pdf"))
  plot_cpm_bar(category_data, cpm_file,
               main_title = paste("CPM -", cat_name))
  
}

message("All done!")

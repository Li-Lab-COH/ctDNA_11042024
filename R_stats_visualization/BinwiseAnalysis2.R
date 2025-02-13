library(dplyr)
library(ggplot2)
library(tidyr)
library(pheatmap)
library(stringr)

###############################################################################
# Read data
###############################################################################
# Replace 'combined_data.csv' with your actual data file path
combined_data <- read.csv("../../data/combined_data_hist.csv", header = TRUE, stringsAsFactors = FALSE)

###############################################################################
# Define bin ranges
###############################################################################

# Last element is removed added 1 to it
# First element is then removed. Begining and ends match in the list locations
bin_breaks <- seq(0, 350, by = 10)  
bin_labels <- paste(head(bin_breaks, -1) + 1, bin_breaks[-1], sep = "-")

combined_data <- combined_data %>%
  mutate(SizeBin = cut(InsertSize,
                       breaks = bin_breaks,
                       labels = bin_labels,
                       include.lowest = TRUE,
                       right = TRUE))

###############################################################################
# Summaries needed for raw frequency vs. CPM
###############################################################################
# We want the frequency aggregated for each (Sample, Bin), plus total for CPM.
binned_data <- combined_data %>%
  group_by(TGen_ID, Patient_ID, Timepoint, Tube, Input, SizeBin) %>%
  summarise(BinFrequency = sum(Frequency), .groups = "drop")

# total frequency per sample
sample_totals <- binned_data %>%
  group_by(TGen_ID) %>%
  summarise(TotalFreq = sum(BinFrequency), .groups = "drop")

# attach total frequency to binned data
binned_data <- binned_data %>%
  left_join(sample_totals, by = "TGen_ID") %>%
  mutate(CPM = (BinFrequency / TotalFreq) * 1e6)

# TODO: I will need to change how the CPM is plotted, it looks like this will
# only plot one bar for each bin, but I want to plot all the bars within the bin

###############################################################################
# Define your comparisons
###############################################################################
comparisons <- list(
  "EDTA_2.5_vs_Streck_2.5" = quote(Tube %in% c("EDTA","Streck") & Input == 2.5),
  "EDTA_10_vs_Streck_10"   = quote(Tube %in% c("EDTA","Streck") & Input == 10)
  # Add other comparisons as needed ...
)

###############################################################################
# Make a base folder
###############################################################################
base_dir <- "Fragment_Size_Analysis_BinWise"
dir.create(base_dir, showWarnings = FALSE)

###############################################################################
# 6. Helper functions to plot
###############################################################################
# A) Z-score matrix for "raw" or "CPM"
zscore_matrix <- function(data_long, value_col) {
  # data_long should have TGen_ID, BinFrequency or CPM, etc.
  # pivot to wide format: rows = TGen_ID, columns = ?
  # Actually, we want to compute Z-scores across *samples* or across *timepoints*?
  # Typically you'd do a heatmap with bins as rows and samples as columns.
  # But if you are only focusing on ONE bin at a time, there's only 1 row -> no real "heatmap".
  # Possibly you want a heatmap across timepoints or across TGen_ID.
  # The example below might break if there's only 1 row or 1 column.
  
  # We'll do: row = Timepoint, col = TGen_ID. That might be more interesting if you have multiple TGen_ID per bin
  wide <- data_long %>%
    pivot_wider(names_from = TGen_ID, values_from = all_of(value_col), values_fill = 0)
  
  # If there's only Timepoint in each row, you can do row_zscore across columns (TGen_ID)
  mat <- as.matrix(wide[, -1])  # remove the Timepoint or whatever your ID is
  rownames(mat) <- wide[[1]]    # the first column is Timepoint
  
  # compute row-wise z-score
  row_zscore <- function(x) (x - mean(x)) / sd(x)
  
  # If there's only 1 row, or 1 sample, standard deviation = 0 => you get NaN. 
  # Let's guard against that:
  if(nrow(mat) < 2 && ncol(mat) < 2) {
    # Not enough data for a meaningful z-score. Return a matrix with the same shape, but no transformations
    return(mat)
  }
  
  z_mat <- t(apply(mat, 1, row_zscore))  # row-wise
  colnames(z_mat) <- colnames(mat)
  rownames(z_mat) <- rownames(mat)
  
  return(z_mat)
}

# B) Heatmap function
plot_zscore_heatmap <- function(z_mat, out_file, main_title) {
  if(nrow(z_mat) <= 1 || ncol(z_mat) <= 1) {
    # If there's no dimension to cluster, skip
    message("Not enough dimensions for a meaningful heatmap. Saving anyway.")
  }
  pheatmap::pheatmap(z_mat,
                     cluster_rows = TRUE,
                     cluster_cols = TRUE,
                     main = main_title,
                     filename = out_file,
                     width = 6,
                     height = 5)
}

# C) Barplot function
#   We'll do a quick barplot across TGen_ID or Timepoint. 
#   Since you might have multiple TGen_ID in each bin, we can do:
plot_bar <- function(data_long, value_col, out_file, main_title) {
  # data_long has columns Timepoint, TGen_ID, Tube, Input, ...
  # We'll group by e.g. Timepoint or TGen_ID. Decide which x-axis is meaningful for you.
  # Example: x-axis = TGen_ID, fill = Timepoint
  p <- ggplot(data_long, aes(x = TGen_ID, y = .data[[value_col]], fill = factor(Timepoint))) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = main_title, x = "Sample (TGen_ID)", y = value_col) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(out_file, p, width = 6, height = 4)
}


###############################################################################
# 7. Main loops
###############################################################################
# We want:
#  For each bin label, create a folder.
#    For each comparison, create (or use) a subfolder.
#      Make 4 plots: 
#         Z-score heatmap of raw, Z-score heatmap of CPM,
#         barplot of raw, barplot of CPM.

for(one_bin_label in bin_labels) {
  
  # Create a folder for this bin
  bin_dir <- file.path(base_dir, paste0("Bin_", one_bin_label))
  dir.create(bin_dir, showWarnings = FALSE)
  
  # Filter data to this bin
  bin_data <- binned_data %>%
    filter(SizeBin == one_bin_label)
  
  # If there's no data for this bin, skip
  if(nrow(bin_data) == 0) {
    message("No data for bin ", one_bin_label, "; skipping.")
    next
  }
  
  # Now loop over each comparison
  for(comp_name in names(comparisons)) {
    
    subfolder <- file.path(bin_dir, comp_name)
    dir.create(subfolder, showWarnings = FALSE)
    
    # Filter data to this comparison
    comp_expr <- comparisons[[comp_name]]
    comp_data <- bin_data %>%
      filter(!!comp_expr)
    
    if(nrow(comp_data) == 0) {
      message("No data for comparison ", comp_name, " in bin ", one_bin_label, "; skipping.")
      next
    }
    
    ### Make the 4 plots ###
    
    # 1) Z-score heatmap of raw values
    raw_z_mat <- zscore_matrix(comp_data, "BinFrequency")
    z_raw_file <- file.path(subfolder, paste0("ZscoreHeatmapRaw_", one_bin_label, "_", comp_name, ".pdf"))
    plot_zscore_heatmap(raw_z_mat, z_raw_file, 
                        main_title = paste0("Z-score Heatmap (Raw) - Bin ", one_bin_label, " - ", comp_name))
    
    # 2) Z-score heatmap of CPM
    cpm_z_mat <- zscore_matrix(comp_data, "CPM")
    z_cpm_file <- file.path(subfolder, paste0("ZscoreHeatmapCPM_", one_bin_label, "_", comp_name, ".pdf"))
    plot_zscore_heatmap(cpm_z_mat, z_cpm_file, 
                        main_title = paste0("Z-score Heatmap (CPM) - Bin ", one_bin_label, " - ", comp_name))
    
    # 3) Barplot of raw values
    raw_bar_file <- file.path(subfolder, paste0("BarplotRaw_", one_bin_label, "_", comp_name, ".pdf"))
    plot_bar(comp_data, "BinFrequency", raw_bar_file, 
             main_title = paste0("Raw Frequency - Bin ", one_bin_label, " - ", comp_name))
    
    # 4) Barplot of CPM
    cpm_bar_file <- file.path(subfolder, paste0("BarplotCPM_", one_bin_label, "_", comp_name, ".pdf"))
    plot_bar(comp_data, "CPM", cpm_bar_file, 
             main_title = paste0("CPM - Bin ", one_bin_label, " - ", comp_name))
    
  } # end for(comp_name...)
} # end for(one_bin_label...)

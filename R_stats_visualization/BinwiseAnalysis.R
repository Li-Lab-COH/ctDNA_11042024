###############################################################################
# TODO: double check cpm calculation and results
# TODO: double check appropriate z-score values
# TODO: ADD A description ehre
###############################################################################



# Libraries
###############################################################################
library(dplyr)
library(tidyr)
library(pheatmap)

###############################################################################
# Read data
###############################################################################
# Replace 'combined_data.csv' with your actual data file path
# combined_data <- read.csv("../../data/combined_data_hist.csv", header = TRUE, stringsAsFactors = FALSE)
combined_data <- read.csv("../../data/combined_data_hist_mice.csv", header = TRUE, stringsAsFactors = FALSE)
###############################################################################
# Binning the Insert Size
###############################################################################
# Last element is removed added 1 to it
# First element is then removed. Begining and ends match in the list locations

bin_breaks <- seq(0, 350, by = 10)  # example
bin_labels <- paste(head(bin_breaks, -1) + 1, bin_breaks[-1], sep = "-")

combined_data <- combined_data %>%
  mutate(SizeBin = cut(InsertSize,
                       breaks = bin_breaks,
                       labels = bin_labels,
                       include.lowest = TRUE,
                       right = TRUE))

###############################################################################
# Aggregate frequencies by (Sample, Bin)
###############################################################################
# identify each sample by TGen_ID, plus keep Tube, Input, Timepoint, etc. for later labeling
# binned_data <- combined_data %>%
#   group_by(TGen_ID, Patient_ID, Timepoint, Tube, Input, SizeBin) %>%
#   summarise(BinFreq = sum(Frequency, na.rm = TRUE), .groups = "drop")

binned_data <- combined_data %>%
  group_by(TGen_ID, treatment, drug_dose, radiation, cell_line, SizeBin) %>%
  summarise(BinFreq = sum(Frequency, na.rm = TRUE), .groups = "drop")

###############################################################################
# Create a CPM column for each row
###############################################################################
# Calculate CPM within each sample, keep raw frequencies

binned_data <- binned_data %>%
  group_by(TGen_ID) %>%
  mutate(TotalFreq = sum(BinFreq, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(CPM = (BinFreq / TotalFreq) * 1e6)

###############################################################################
# Create Two Wide Matrices (Raw vs. CPM)
###############################################################################
# Create row labels
# binned_data <- binned_data %>%
#   mutate(SampleLabel = paste(Timepoint, Tube, Input, sep="|"))
binned_data <- binned_data %>%
    mutate(SampleLabel = paste(cell_line, treatment, 
                               paste0(drug_dose, "µM"), paste0(radiation, "GY"),  
                               sep="|"))

unique(binned_data$SampleLabel)


# For RAW frequencies
raw_wide <- binned_data %>%
  select(SampleLabel, SizeBin, BinFreq) %>%
  pivot_wider(names_from = SizeBin, values_from = BinFreq, values_fill = 0)

# For CPM
cpm_wide <- binned_data %>%
  select(SampleLabel, SizeBin, CPM) %>%
  pivot_wider(names_from = SizeBin, values_from = CPM, values_fill = 0)

# Turn them into matrices
raw_mat <- as.matrix(raw_wide[, -1])  # remove SampleLabel col
rownames(raw_mat) <- raw_wide$SampleLabel

cpm_mat <- as.matrix(cpm_wide[, -1])
rownames(cpm_mat) <- cpm_wide$SampleLabel

###############################################################################
# Compute Column-wise Z-scores (across samples) for each matrix
###############################################################################
# Which samples are “above/below average” for each bin.
col_zscore <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

z_raw_mat <- apply(raw_mat, 2, col_zscore)
z_cpm_mat <- apply(cpm_mat, 2, col_zscore)

###############################################################################
# F. Define Custom Row Orders
###############################################################################


custom_order_mice <- c(
  "PE|DMSO|0µM|0GY",
  "PE|DMSO|0µM|4GY",
  "PE|enzalutamide|5µM|0GY",
  "PE|enzalutamide|5µM|4GY",
  "PLV8|media|0µM|0GY",
  "PLV8|media|0µM|4GY",
  "PLV8|enzalutamide|10µM|0GY",
  "PLV8|enzalutamide|10µM|4GY"
)


custom_order_mice_1 <- c(
  "PE|DMSO|0µM|0GY",
  "PE|enzalutamide|5µM|0GY",
  "PLV8|media|0µM|0GY",
  "PLV8|enzalutamide|10µM|0GY",
  "PE|DMSO|0µM|4GY",
  "PE|enzalutamide|5µM|4GY",
  "PLV8|media|0µM|4GY",
  "PLV8|enzalutamide|10µM|4GY"
)


# custom_order_1 <- c(
#   "1|EDTA|2.5",
#   "2|EDTA|2.5",
#   "3|EDTA|2.5",
#   "4|EDTA|2.5",
#   "1|EDTA|10",
#   "2|EDTA|10",
#   "3|EDTA|10",
#   "4|EDTA|10",
#   "5|EDTA|10",
#   "1|Streck|2.5",
#   "2|Streck|2.5",
#   "4|Streck|2.5",
#   "1|Streck|10",
#   "2|Streck|10",
#   "4|Streck|10",
#   "5|Streck|10"
# )
# 
# custom_order_2 <- c(
#   "1|EDTA|2.5",
#   "1|EDTA|10",
#   "1|Streck|2.5",
#   "1|Streck|10",
#   "2|EDTA|2.5",
#   "2|EDTA|10",
#   "2|Streck|10",
#   "2|Streck|2.5",
#   "3|EDTA|2.5",
#   "3|EDTA|10",
#   "4|EDTA|2.5",
#   "4|EDTA|10",
#   "4|Streck|2.5",
#   "4|Streck|10",
#   "5|EDTA|10",
#   "5|Streck|10"
# )


# Ensure we keep only rows that exist in the matrix
custom_order_mice <- custom_order_mice[custom_order_mice %in% rownames(z_raw_mat)]
custom_order_mice_1 <- custom_order_mice_1[custom_order_mice_1 %in% rownames(z_raw_mat)]

# Subset and reorder the matrices
z_raw_mat_1 <- z_raw_mat[custom_order_mice, , drop=FALSE]
z_cpm_mat_1 <- z_cpm_mat[custom_order_mice, , drop=FALSE]

z_raw_mat_2 <- z_raw_mat[custom_order_mice_1, , drop=FALSE]
z_cpm_mat_2 <- z_cpm_mat[custom_order_mice_1, , drop=FALSE]


# # Ensure we keep only rows that exist in the matrix
# custom_order_1 <- custom_order_1[custom_order_1 %in% rownames(z_raw_mat)]
# custom_order_2 <- custom_order_2[custom_order_2 %in% rownames(z_raw_mat)]
# 
# # Subset and reorder the matrices
# z_raw_mat_1 <- z_raw_mat[custom_order_1, , drop=FALSE]
# z_cpm_mat_1 <- z_cpm_mat[custom_order_1, , drop=FALSE]
# 
# z_raw_mat_2 <- z_raw_mat[custom_order_2, , drop=FALSE]
# z_cpm_mat_2 <- z_cpm_mat[custom_order_2, , drop=FALSE]


###############################################################################
# Removing 
###############################################################################
# # Define the row to be removed
# row_to_remove <- "4|EDTA|10"
# 
# # Remove from all matrices
# z_raw_mat <- z_raw_mat[!rownames(z_raw_mat) %in% row_to_remove, , drop=FALSE]
# # z_cpm_mat <- z_cpm_mat[!rownames(z_cpm_mat) %in% row_to_remove, , drop=FALSE]
# 
# z_raw_mat_1 <- z_raw_mat_1[!rownames(z_raw_mat_1) %in% row_to_remove, , drop=FALSE]
# # z_cpm_mat_1 <- z_cpm_mat_1[!rownames(z_cpm_mat_1) %in% row_to_remove, , drop=FALSE]
# 
# z_raw_mat_2 <- z_raw_mat_2[!rownames(z_raw_mat_2) %in% row_to_remove, , drop=FALSE]
# # z_cpm_mat_2 <- z_cpm_mat_2[!rownames(z_cpm_mat_2) %in% row_to_remove, , drop=FALSE]


###############################################################################
# Plot the Heatmaps
###############################################################################

# # 1) Default clustering heatmaps
# pheatmap(z_raw_mat, main="Z-score Heatmap (Raw Frequencies)")
# pheatmap(z_cpm_mat, main="Z-score Heatmap (CPM)")
# 
# # 2) Manually ordered heatmaps (set 1)
pheatmap(z_raw_mat_1, cluster_rows=FALSE, cluster_cols = FALSE,
         main="Z-score Heatmap (Raw)")
pheatmap(z_cpm_mat_1, cluster_rows=FALSE, cluster_cols = FALSE,
         main="Z-score Heatmap (CPM) - Paired")


# 3) Manually ordered heatmaps (set 2)
pheatmap(z_raw_mat_2, cluster_rows=FALSE, cluster_cols = FALSE,
         main="Z-score Heatmap (Raw)")
pheatmap(z_cpm_mat_2, cluster_rows=FALSE, cluster_cols = FALSE,
         main="Z-score Heatmap (CPM, bin-wise)")


# Save locations:
# raw_image <- "../../figures/Heatmap_raw_bin10.png"
# CPM_image <- "../../figures/Heatmap_CPM_bin10.png"
# 
# png(raw_image, width = 3000, height = 2000, res = 300)
# pheatmap(z_raw_mat_2, cluster_rows=FALSE, cluster_cols=FALSE,
#          main="Z-score Heatmap (Raw)")
# dev.off()
# 
# png(CPM_image, width = 3000, height = 2000, res = 300)
# pheatmap(z_cpm_mat_2, cluster_rows=FALSE, cluster_cols=FALSE,
#          main="Z-score Heatmap (CPM)")
# dev.off()

# dev.cur()  # Check which graphics device is active
# dev.new()


#-------------------------rowwise heatmap-----------------------------------

# Compute row-wise z-scores (across bins) for each matrix
row_zscore <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

# For the raw frequency matrix:
z_raw_mat_samples <- t(apply(raw_mat, 1, row_zscore))
# For the CPM matrix:
z_cpm_mat_samples <- t(apply(cpm_mat, 1, row_zscore))

# re-order the data
z_raw_mat_samples_1 <- z_raw_mat_samples[custom_order_mice, , drop=FALSE]
z_cpm_mat_samples_1 <- z_cpm_mat_samples[custom_order_mice, , drop=FALSE]

z_raw_mat_samples_2 <- z_raw_mat_samples[custom_order_mice_1, , drop=FALSE]
z_cpm_mat_samples_2 <- z_cpm_mat_samples[custom_order_mice_1, , drop=FALSE]

rm(z_raw_mat_1)

# Plot the heatmaps using pheatmap:
pheatmap(z_raw_mat_samples_1, cluster_rows = FALSE, cluster_cols = FALSE,
         main = "Z-score Heatmap (Raw, Row-wise)")
pheatmap(z_cpm_mat_samples_1, cluster_rows = FALSE, cluster_cols = FALSE,
         main = "Z-score Heatmap (CPM, Row-wise)")

pheatmap(z_raw_mat_samples_2, cluster_rows = FALSE, cluster_cols = FALSE,
         main = "Z-score Heatmap (Raw Frequencies, Row-wise)")
pheatmap(z_cpm_mat_samples_2, cluster_rows = FALSE, cluster_cols = FALSE,
         main = "Z-score Heatmap (CPM, Row-wise)")



# ------------------------Mice Images to save------------------------------

save_location <- "../../figures/mice_bin_fragments/"

png(raw_image, width = 3000, height = 2000, res = 300)


# bin-wise
pheatmap(z_cpm_mat_1, cluster_rows=FALSE, cluster_cols = FALSE,
         main="Z-score Heatmap (CPM) - Paired - Bin-wise")
pheatmap(z_cpm_mat_2, cluster_rows=FALSE, cluster_cols = FALSE,
         main="Z-score Heatmap (CPM, Bin-wise)")


# Row wise heatmap
pheatmap(z_cpm_mat_samples_1, cluster_rows = FALSE, cluster_cols = FALSE,
         main = "Z-score Heatmap (CPM, Row-wise)")

pheatmap(z_cpm_mat_samples_2, cluster_rows = FALSE, cluster_cols = FALSE,
         main = "Z-score Heatmap (CPM, Row-wise)")





# Define the save location folder
save_location <- "../../figures/mice_bin_fragments/"

# Define file names based on the analytical results
file_binwise_paired <- paste0(save_location, "Zscore_CPM_Paired_BinWise.png")
file_binwise <- paste0(save_location, "Zscore_CPM_BinWise.png")
file_rowwise_paired <- paste0(save_location, "Zscore_CPM_Paired_RowWise.png")
file_rowwise <- paste0(save_location, "Zscore_CPM_RowWise.png")

#----------------- Save the Bin-wise Heatmaps -----------------

# Save the CPM heatmap for the paired custom order (bin-wise)
png(file_binwise_paired, width = 3000, height = 2000, res = 300)
pheatmap(z_cpm_mat_1, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         main = "Z-score Heatmap (CPM) - Paired - Bin-wise")
dev.off()

# Save the CPM heatmap for the alternative custom order (bin-wise)
png(file_binwise, width = 3000, height = 2000, res = 300)
pheatmap(z_cpm_mat_2, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         main = "Z-score Heatmap (CPM) - Bin-wise")
dev.off()

#----------------- Save the Row-wise Heatmaps -----------------

# Save the CPM heatmap for the paired custom order (row-wise)
png(file_rowwise_paired, width = 3000, height = 2000, res = 300)
pheatmap(z_cpm_mat_samples_1, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         main = "Z-score Heatmap (CPM) - Paired - Row-wise")
dev.off()

# Save the CPM heatmap for the alternative custom order (row-wise)
png(file_rowwise, width = 3000, height = 2000, res = 300)
pheatmap(z_cpm_mat_samples_2, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         main = "Z-score Heatmap (CPM) - Row-wise")
dev.off()






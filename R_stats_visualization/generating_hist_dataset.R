# Load necessary libraries
library(tidyverse)
library(data.table)
library(ggplot2)
library(gridExtra)
library(grid)
library(rstudioapi)

#----------------------Prepping -----------------

#setting current wd to files location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Define directories
metadata_file <- "/Users/janzules/Roselab/ctDNA_11042024/docs/samples_metadata.csv"
qualimap_dir <- "/Users/janzules/Roselab/ctDNA_11042024/results/alignment_metrics/qualimap"

# Load Metadata
metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)

# Get a list of all samples' raw data folders
sample_dirs <- list.dirs(qualimap_dir, recursive = FALSE)

# Create an empty list to hold data
insert_size_data <- list()

#-----------------loop-----------------

# Loop through each sample directory
for (dir in sample_dirs) {
  # Extract the full directory name
  full_name <- basename(dir)
  print(paste0("Processing: ", full_name))
  
  # Extract the TGen_ID from the directory name
  tgen_id <- sub("_qcReport$", "", full_name)
  print(paste0("Extracted TGen_ID: ", tgen_id))
  
  # Define the path to the insert size histogram file
  insert_size_file <- file.path(dir, "raw_data_qualimapReport/insert_size_histogram.txt")
  print(paste0("Insert Size File Path: ", insert_size_file))
  
  # Check if the file exists
  if (file.exists(insert_size_file)) {
    # Read the histogram data, skipping the comment row
    data <- fread(insert_size_file, skip = 1, header = FALSE)
    colnames(data) <- c("InsertSize", "Frequency")  # Rename columns for clarity
    
    # Add the TGen_ID to the data
    data$TGen_ID <- tgen_id
    print(head(data))  # Optional: Print the first few rows of the data for verification
    
    # Append the data to the list
    insert_size_data[[tgen_id]] <- data
  } else {
    print(paste0("File not found for TGen_ID: ", tgen_id))
  }
}

# Verify the structure of the combined data
# print(insert_size_data)


#------------Combining ---------------------

# Combine all histogram data into one dataframe
insert_size_df <- bind_rows(insert_size_data)

# Join metadata with the histogram data by TGen_ID
combined_data <- left_join(insert_size_df, metadata, by = "TGen_ID")

# write.csv(combined_data, "../../data/combined_data_hist.csv")

#---------------Figures----------------

# Combining the figures
# Ensure your dataset is a data.table for easier manipulation
setDT(combined_data)

# Define a color palette for the five timepoints
timepoint_colors <- scale_color_manual(
  values = c("1" = "red", "2" = "blue", "3" = "green", "4" = "orange", "5" = "purple")
)

# Create the plot
ggplot(combined_data, aes(x = InsertSize, y = Frequency, color = as.factor(Timepoint))) +
  geom_line(size = 1) +  # Line plot for frequency across insert sizes
  facet_wrap(~ Tube, ncol = 1, scales = "free_y") +  # One plot per tube type
  timepoint_colors +  # Apply custom colors
  theme_minimal() +  # Use a clean theme
  labs(
    title = "Insert Size Distribution by Tube Type and Timepoint",
    x = "Insert Size (bp)",
    y = "Frequency",
    color = "Timepoint"
  ) +
  theme(
    strip.text = element_text(size = 12, face = "bold"),  # Facet titles (tube types)
    axis.title = element_text(size = 12),  # Axis titles
    legend.title = element_text(size = 12),  # Legend title
    legend.text = element_text(size = 10)  # Legend text
  )


# Step 6: Perform Statistical Analysis (e.g., comparing means of insert sizes)
stats <- combined_data %>%
  group_by(Tube, Patient_ID) %>%
  summarize(mean_insert_size = weighted.mean(InsertSize, Frequency, na.rm = TRUE),
            total_reads = sum(Frequency, na.rm = TRUE))

print(stats)



#-------------------Normalizing Data proportional Scaling----------------

# Normalize by total counts
normalized_data <- combined_data[, .(
  Timepoint,
  InsertSize,
  Frequency,
  Tube,
  NormalizedFrequency = Frequency / sum(Frequency)  # Calculate NormalizedFrequency
), by = TGen_ID]

head(normalized_data)

# Ensure your dataset is a data.table for easier manipulation
setDT(normalized_data) 

# Define a color palette for the five timepoints
timepoint_colors <- scale_color_manual(
  values = c("1" = "red", "2" = "blue", "3" = "green", "4" = "orange", "5" = "purple")
)

# Create the plot
ggplot(normalized_data, aes(x = InsertSize, y = NormalizedFrequency, color = as.factor(Timepoint))) +
  geom_line(size = 1) +  # Line plot for frequency across insert sizes
  facet_wrap(~ Tube, ncol = 1, scales = "free_y") +  # One plot per tube type
  timepoint_colors +  # Apply custom colors
  theme_minimal() +  # Use a clean theme
  labs(
    title = "Insert Size Distribution (Frequency)",
    x = "Insert Size (bp)",
    y = "Frequency",
    color = "Timepoint"
  ) +
  theme(
    strip.text = element_text(size = 12, face = "bold"),  # Facet titles (tube types)
    axis.title = element_text(size = 12),  # Axis titles
    legend.title = element_text(size = 12),  # Legend title
    legend.text = element_text(size = 10)  # Legend text
  )


#----------------------CPM----------------------
# Calculate Counts Per Million (CPM)
cpm_data <- combined_data[, .(
  Timepoint,
  InsertSize,
  Frequency,
  Tube,
  CPM = (Frequency / sum(Frequency)) * 1e6  # Calculate CPM
), by = TGen_ID]

# Ensure the dataset is a data.table for easier manipulation
setDT(cpm_data)

# Define a color palette for the five timepoints
timepoint_colors <- scale_color_manual(
  values = c("1" = "red", "2" = "blue", "3" = "green", "4" = "orange", "5" = "purple")
)

# Create the plot
all_hist = ggplot(cpm_data, aes(x = InsertSize, y = CPM, color = as.factor(Timepoint), group = TGen_ID)) +
  geom_line(size = 0.5) +  # Line plot for CPM across insert sizes
  facet_wrap(~ Tube, ncol = 1, scales = "free_y") +  # One plot per tube type
  timepoint_colors +  # Apply custom colors
  theme_minimal() +  # Use a clean theme
  labs(
    title = "Insert Size Distribution (Counts Per Million)",
    x = "Insert Size (bp)",
    y = "Counts Per Million (CPM)",
    color = "Timepoint"
  ) +
  theme(
    strip.text = element_text(size = 12, face = "bold"),  # Facet titles (tube types)
    axis.title = element_text(size = 12),  # Axis titles
    legend.title = element_text(size = 12),  # Legend title
    legend.text = element_text(size = 10),  # Legend text
    plot.background = element_rect(fill = "white", color = NA),  # White plot background
    panel.background = element_rect(fill = "white", color = NA)  # White panel background
  )


# Ensure each sample has its own line
all_hist_50_150 = ggplot(cpm_data, aes(x = InsertSize, y = CPM, color = as.factor(Timepoint), group = TGen_ID)) +
  geom_line(size = 0.5) +  # Line plot for CPM across insert sizes
  facet_wrap(~ Tube, ncol = 1, scales = "free_y") +  # One plot per tube type
  timepoint_colors +  # Apply custom colors
  theme_minimal() +  # Use a clean theme
  labs(
    title = "Insert Size Distribution (Counts Per Million)",
    x = "Insert Size (bp)",
    y = "Counts Per Million (CPM)",
    color = "Timepoint"
  ) +
  theme(
    strip.text = element_text(size = 12, face = "bold"),  # Facet titles (tube types)
    axis.title = element_text(size = 12),  # Axis titles
    legend.title = element_text(size = 12),  # Legend title
    legend.text = element_text(size = 10),  # Legend text
    plot.background = element_rect(fill = "white", color = NA),  # White plot background
    panel.background = element_rect(fill = "white", color = NA)  # White panel background
  )+
  xlim(50,150)+
  ylim(0, 15000)


#################
# Saving

# Save Plot 1
ggsave(
  filename = "/Users/janzules/Roselab/ctDNA_11042024/figures/insert_size_distr.png",
  plot = all_hist,
  width = 6,
  height = 4,
  dpi = 300
)

# Save Plot 2
ggsave(
  filename = "/Users/janzules/Roselab/ctDNA_11042024/figures/insert_size_distr_50_150.png",
  plot = all_hist_50_150,
  width = 6,
  height = 4,
  dpi = 300
)


#------------New plot - dash, summary-------------------------

# Calculate Counts Per Million (CPM)
cpm_data <- combined_data[, .(
  Timepoint,
  InsertSize,
  Frequency,
  Tube,
  Input,
  CPM = (Frequency / sum(Frequency)) * 1e6  # Calculate CPM
), by = TGen_ID]

# Ensure the dataset is a data.table
setDT(cpm_data)

# Define a color palette for the five timepoints
timepoint_colors <- scale_color_manual(
  values = c("1" = "red", "2" = "blue", "3" = "green", "4" = "orange", "5" = "purple")
)

# Define linetype based on Input column
linetype_mapping <- scale_linetype_manual(
  values = c("2.5" = "dashed", "10" = "solid")
)

# Step 1: Calculate the Summary Table
summary_table <- combined_data[, .(
  MedianFrequency = median(Frequency)  # Calculate median Frequency
), by = .(Patient_ID_Input = paste0(Patient_ID, "_", Input))]  # Combine Patient_ID and Input

# Step 2: Convert Summary Table to a Grob
summary_table_grob <- tableGrob(
  summary_table, 
  rows = NULL,  # No row names
  theme = ttheme_minimal(
    core = list(fg_params = list(cex = 0.8)),  # Adjust font size
    colhead = list(fg_params = list(cex = 0.8, fontface = "bold"))
  )
)


# Create the plot
all_hist <- ggplot(cpm_data, aes(x = InsertSize, y = CPM, 
                                 color = as.factor(Timepoint), 
                                 group = TGen_ID, 
                                 linetype = as.factor(Input))) + 
  geom_line(size = 0.2) +  # Line plot for CPM across insert sizes
  facet_wrap(~ Tube, ncol = 1, scales = "free_y") +  # One plot per tube type
  timepoint_colors +  # Apply custom colors
  linetype_mapping +  # Apply custom linetypes
  theme_minimal() +  # Use a clean theme
  labs(
    title = "Insert Size Distribution (Counts Per Million)",
    x = "Insert Size (bp)",
    y = "Counts Per Million (CPM)",
    color = "Timepoint",
    linetype = "Input"  # Legend label for linetype
  ) +
  theme(
    strip.text = element_text(size = 12, face = "bold"),  # Facet titles (tube types)
    axis.title = element_text(size = 12),  # Axis titles
    legend.title = element_text(size = 12),  # Legend title
    legend.text = element_text(size = 10),  # Legend text
    plot.background = element_rect(fill = "white", color = NA),  # White plot background
    panel.background = element_rect(fill = "white", color = NA)  # White panel background
  )



all_hist_50_150 = ggplot(cpm_data, aes(x = InsertSize, y = CPM, 
                                       color = as.factor(Timepoint), 
                                       group = TGen_ID, 
                                       linetype = as.factor(Input))) + 
  geom_line(size = 0.2) +  # Line plot for CPM across insert sizes
  facet_wrap(~ Tube, ncol = 1, scales = "free_y") +  # One plot per tube type
  timepoint_colors +  # Apply custom colors
  linetype_mapping +  # Apply custom linetypes
  theme_minimal() +  # Use a clean theme
  labs(
    title = "Insert Size Distribution (Counts Per Million)",
    x = "Insert Size (bp)",
    y = "Counts Per Million (CPM)",
    color = "Timepoint",
    linetype = "Input"  # Legend label for linetype
  ) +
  theme(
    strip.text = element_text(size = 12, face = "bold"),  # Facet titles (tube types)
    axis.title = element_text(size = 12),  # Axis titles
    legend.title = element_text(size = 12),  # Legend title
    legend.text = element_text(size = 10),  # Legend text
    plot.background = element_rect(fill = "white", color = NA),  # White plot background
    panel.background = element_rect(fill = "white", color = NA)  # White panel background
  )+
  xlim(50,150)+
  ylim(0, 15000)

# Print the plot
print(all_hist)
print(all_hist_50_150)


#################
# Saving

# Save Plot 1 - High Resolution
ggsave(
  filename = "/Users/janzules/Roselab/ctDNA_11042024/figures/insert_size_distr.png",
  plot = all_hist,
  width = 12,     # Double the width (e.g., from 6 to 12)
  height = 8,     # Double the height (e.g., from 4 to 8)
  dpi = 2400,      # Increase DPI (300 to 600 or even 1200 for very high quality)
  units = "in"    # Ensure dimensions are in inches
)

# Save Plot 2 - High Resolution
ggsave(
  filename = "/Users/janzules/Roselab/ctDNA_11042024/figures/insert_size_distr_50_150.png",
  plot = all_hist_50_150,
  width = 12,     # Double the width
  height = 8,     # Double the height
  dpi = 2400,      # Increase DPI
  units = "in"    # Units in inches
)

















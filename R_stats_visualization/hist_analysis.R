library(matrixStats)
library(rstudioapi)
library(data.table)
library(writexl)

#setting current wd to files location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))



#----------------Preparing data------------------------
combined_data <- read.csv("../../data/combined_data_hist.csv")
setDT(combined_data)

# Normalize by total counts
# normalized_data <- combined_data[, .(
#   InsertSize,
#   Frequency,
#   Patient_ID,
#   Timepoint,
#   Tube,
#   Input,
#   NormalizedFrequency = Frequency / sum(Frequency)  # Calculate NormalizedFrequency
# ), by = TGen_ID]

cpm_data <- combined_data[, .(
  InsertSize,
  Frequency,
  Patient_ID,
  Timepoint,
  Tube,
  Input,
  CPM = (Frequency / sum(Frequency)) * 1e6  # Calculate CPM
), by = TGen_ID]



#----------------Summary stats--------------
# I used frequency counts to calculate these stats, cpm produces the same results (duh, but i had to double check)

# Ensure cpm_data is a data.table
setDT(cpm_data)

# Create the Patient_ID_Input column
cpm_data[, Patient_ID_Input := paste0(Patient_ID, "_", Input)]

# Calculate weighted mean, weighted median, and percentages in one go
summary_stats <- cpm_data[, .(
  WeightedMean = weighted.mean(InsertSize, w = Frequency, na.rm = TRUE),       # Weighted Mean
  WeightedMedian = weightedMedian(InsertSize, w = Frequency, na.rm = TRUE),    # Weighted Median
  Total_Reads = sum(Frequency, na.rm = TRUE),                                    # Total Frequency
  Below150_Percent = 100 * sum(Frequency[InsertSize < 150], na.rm = TRUE) / sum(Frequency, na.rm = TRUE), # % Below 150
  Above160_Percent = 100 * sum(Frequency[InsertSize > 160], na.rm = TRUE) / sum(Frequency, na.rm = TRUE)  # % Above 160
), by = Patient_ID_Input]

# Round values for readability
summary_stats[, `:=`(
  WeightedMean = round(WeightedMean),
  WeightedMedian = round(WeightedMedian),
  Total_Reads = formatC(Total_Reads, format = "e", digits = 2),
  Below150_Percent = round(Below150_Percent, 2),
  Above160_Percent = round(Above160_Percent, 2)
)]

# View the result
print(summary_stats)


write_xlsx(summary_stats, "../../results/alignment_metrics/mapping_QC_analysis/summary_hist_stats.xlsx")


write.csv(median_fragment_sizes, "../../results/alignment_metrics/mapping_QC_analysis/median_fragment_sizes")


#-----------------------Binning Outputs--------------------------------------




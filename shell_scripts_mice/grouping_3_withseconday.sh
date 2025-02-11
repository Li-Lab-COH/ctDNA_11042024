#!/bin/bash
#SBATCH --job-name=group_reads_3        # Job name
#SBATCH --output=/home/janzules/ctDNA_11042024/code/slurmOutput/groupReads_3.out  # Standard output log
#SBATCH --error=/home/janzules/ctDNA_11042024/code/slurmOutput/groupReads_3.err   # Standard error log
#SBATCH --ntasks=1                     # Number of tasks (processes)
#SBATCH --cpus-per-task=8              # Number of CPU cores per task
#SBATCH --mem=32G                      # Memory allocation
#SBATCH --time=12:00:00                # Time limit (hh:mm:ss)

# Load necessary modules
module load Mamba

# Activate the Mamba environment
mamba activate my_fgbio_env

# Confirm the environment is active (optional, for debugging)
echo "Activated environment: $(conda info --envs | grep '*' | awk '{print $1}')"

# Define directories and files
mapped_bam_dir="/home/janzules/ctDNA_11042024/data/consensusPipeline/intermediate_files/mapped_bam"
grouped_bam_dir="/home/janzules/ctDNA_11042024/data/consensusPipeline/intermediate_files/grouped_bam_3"
histograms_dir="/home/janzules/ctDNA_11042024/data/consensusPipeline/output/histograms_3"
logs_dir="/home/janzules/ctDNA_11042024/data/consensusPipeline/output/logs/grouping_3"
anumber_file="/home/janzules/ctDNA_11042024/code/addresses/human_Anumbers.txt"

# Create necessary directories if they do not exist
mkdir -p "$grouped_bam_dir" "$histograms_dir" "$logs_dir"

# Process each anumber
while read -r anumber; do
    mapped_bam_file="$mapped_bam_dir/${anumber}.mapped.bam"
    grouped_bam_file="$grouped_bam_dir/${anumber}.grouped.bam"
    histogram_file="$histograms_dir/${anumber}.family-size-histogram.txt"
    log_file="$logs_dir/group_reads.${anumber}.log"

    if [ -f "$mapped_bam_file" ]; then
        echo "Processing: $anumber"

        # Run GroupReadsByUmi
        fgbio -Xmx8g --compression 1 --async-io GroupReadsByUmi \
            --input "$mapped_bam_file" \
            --strategy Adjacency \
            --edits 2 \
            --include-secondary TRUE\
            --include-supplementary  TRUE\
            --threads 8 \
            --output "$grouped_bam_file" \
            --family-size-histogram "$histogram_file" &> "$log_file"

        if [ $? -eq 0 ]; then
            echo "Grouping completed successfully for: $anumber"
        else
            echo "Error during grouping for: $anumber. Check $log_file for details."
        fi
    else
        echo "Warning: Missing mapped BAM file for $anumber. Skipping."
    fi
done < "$anumber_file"

echo "GroupReadsByUmi pipeline complete."

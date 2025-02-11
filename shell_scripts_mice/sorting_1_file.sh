#!/bin/bash
#SBATCH --job-name=sort_bam_single         # Job name
#SBATCH --output=/home/janzules/ctDNA_11042024/code/slurmOutput/sort_bam_single.out  # Standard output log
#SBATCH --error=/home/janzules/ctDNA_11042024/code/slurmOutput/sort_bam_single.err   # Standard error log
#SBATCH --ntasks=1                         # Number of tasks (processes)
#SBATCH --cpus-per-task=8                  # Number of CPU cores per task
#SBATCH --mem=300G                          # Memory allocation
#SBATCH --time=12:00:00                    # Time limit (hh:mm:ss)

# Load samtools module
module load samtools

# Set the specific anumber to process
anumber="A14903"  # Change this to your specific anumber

# Define directories and files
mapped_bam_dir="/home/janzules/ctDNA_11042024/data/consensusPipeline/intermediate_files/mapped_bam"
sorted_bam_dir="/home/janzules/ctDNA_11042024/data/consensusPipeline/intermediate_files/sorted_bam"

# Define input and output BAM file names
mapped_bam_file="$mapped_bam_dir/${anumber}.mapped.bam"
sorted_bam_file="$sorted_bam_dir/${anumber}.sorted.bam"

# Check if the input BAM file exists
if [ -f "$mapped_bam_file" ]; then
    echo "Sorting BAM file for: $anumber"
    samtools sort -@ 8 -m 32G -o "$sorted_bam_file" "$mapped_bam_file"

    if [ $? -eq 0 ]; then
        echo "Sorting successful for: $anumber"
    else
        echo "Error during sorting for: $anumber"
    fi
else
    echo "Error: Missing mapped BAM file for $anumber. Exiting."
    exit 1
fi

echo "Sorting complete for: $anumber"

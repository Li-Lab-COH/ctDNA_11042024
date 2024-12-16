#!/bin/bash
#SBATCH --job-name=sort_bam_array         # Job name
#SBATCH --output=/home/janzules/ctDNA_11042024/code/slurmOutput/sort_jobarray/sort_bam_array_%A_%a.out  # Standard output log
#SBATCH --error=/home/janzules/ctDNA_11042024/code/slurmOutput/sort_jobarray/sort_bam_array_%A_%a.err   # Standard error log
#SBATCH --ntasks=1                        # Number of tasks (processes)
#SBATCH --cpus-per-task=8                 # Number of CPU cores per task
#SBATCH --mem=96G                         # Memory allocation per job
#SBATCH --time=12:00:00                   # Time limit (hh:mm:ss)
#SBATCH --array=0-7                       # Job array index (8 jobs for 16 files, 2 files per job)

# Load samtools module
module load samtools

# Define directories and files
mapped_bam_dir="/home/janzules/ctDNA_11042024/data/consensusPipeline/intermediate_files/mapped_bam"
sorted_bam_dir="/home/janzules/ctDNA_11042024/data/consensusPipeline/output/sorted_bam"
anumber_file="/home/janzules/ctDNA_11042024/code/addresses/human_Anumbers.txt"

# Create output directory if it does not exist
mkdir -p "$sorted_bam_dir"

# Read anumber_file into an array
mapfile -t anumbers < "$anumber_file"

# Calculate the range of anumbers to process for this job array task
start_index=$(( SLURM_ARRAY_TASK_ID * 2 ))
end_index=$(( start_index + 1 ))

# Ensure end_index does not exceed the number of available anumbers
if (( end_index >= ${#anumbers[@]} )); then
    end_index=$(( ${#anumbers[@]} - 1 ))
fi

# Process each anumber in the range for this job
for i in $(seq $start_index $end_index); do
    anumber="${anumbers[$i]}"
    mapped_bam_file="$mapped_bam_dir/${anumber}.mapped.bam"
    sorted_bam_file="$sorted_bam_dir/${anumber}.sorted.bam"

    if [ -f "$mapped_bam_file" ]; then
        echo "Sorting BAM file for: $anumber"
        samtools sort -@ 8 -m 12G -o "$sorted_bam_file" "$mapped_bam_file"

        if [ $? -eq 0 ]; then
            echo "Sorting successful for: $anumber"
        else
            echo "Error during sorting for: $anumber"
        fi
    else
        echo "Warning: Missing mapped BAM file for $anumber. Skipping."
    fi
done

echo "Job $SLURM_ARRAY_TASK_ID complete"

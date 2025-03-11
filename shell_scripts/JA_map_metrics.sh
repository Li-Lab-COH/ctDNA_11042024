#!/bin/bash
#SBATCH --job-name=alignment_metrics_array  # Job name
#SBATCH --output=/home/janzules/ctDNA_11042024/code/slurmOutput/mapping_qc/alignment_metrics_%A_%a.out  # Standard output log
#SBATCH --error=/home/janzules/ctDNA_11042024/code/slurmOutput/mapping_qc/alignment_metrics_%A_%a.err   # Standard error log
#SBATCH --ntasks=1                          # Number of tasks (processes)
#SBATCH --cpus-per-task=16                  # Number of CPU cores per task
#SBATCH --mem=128G                          # Memory allocation
#SBATCH --time=12:00:00                     # Time limit (hh:mm:ss)
#SBATCH --array=0-7                         # Job array index (8 jobs for 16 files, 2 files per job)

# Loading modules
module load R
module load Qualimap

# Define directories and files
sorted_bam_dir="/home/janzules/ctDNA_11042024/data/consensusPipeline/intermediate_files/sorted_bam"
anumber_file="/home/janzules/ctDNA_11042024/code/addresses/human_Anumbers.txt"
output_dir="/home/janzules/ctDNA_11042024/data/consensusPipeline/output/alignment_metrics"

# Additional folders
insert_size_dir="$output_dir/insert_size"
qualmap_dir="$output_dir/qualimap"

# Create necessary directories if they do not exist
mkdir -p "$output_dir" "$insert_size_dir" "$qualmap_dir"

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
    sorted_bam_file="$sorted_bam_dir/${anumber}.sorted.bam"
    IS_txt="$insert_size_dir/${anumber}.insert_size_metrics.txt"
    IS_hist="$insert_size_dir/${anumber}.insert_size_histogram.pdf"
    qualimap_out="$qualmap_dir/${anumber}_qcReport"

    if [ -f "$sorted_bam_file" ]; then
        echo "Processing alignment metrics for: $anumber"

        # # Run Picard CollectInsertSizeMetrics
        # java -jar /opt/picard/2.21.1/picard.jar CollectInsertSizeMetrics \
        #     I=$sorted_bam_file \
        #     O=$IS_txt \
        #     H=$IS_hist \
        #     M=0.5

        # if [ $? -eq 0 ]; then
        #     echo "Insert size metrics successful for: $anumber"
        # else
        #     echo "Error during Insert size metrics collection for: $anumber"
        # fi

        # Run Qualimap QC
        echo "Starting quality control for: $anumber"
        qualimap bamqc \
            -bam $sorted_bam_file \
            -outdir $qualimap_out \
            -nt 16 \
            --java-mem-size=75G

        if [ $? -eq 0 ]; then
            echo "Quality control successful for: $anumber"
        else
            echo "Error during quality control for: $anumber"
        fi

    else
        echo "Warning: Missing sorted BAM file for $anumber. Skipping."
    fi
done

echo "Job $SLURM_ARRAY_TASK_ID complete"

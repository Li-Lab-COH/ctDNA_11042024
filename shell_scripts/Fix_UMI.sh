#!/bin/bash
#SBATCH --job-name=umi_correction_job    # Job name
#SBATCH --output=/home/janzules/ctDNA_11042024/code/slurmOutput/umiCorrection.out   # Standard output log
#SBATCH --error=/home/janzules/ctDNA_11042024/code/slurmOutput/umiCorrection.err    # Standard error log
#SBATCH --ntasks=1                      # Number of tasks (processes)
#SBATCH --cpus-per-task=4               # Number of CPU cores per task
#SBATCH --mem=16G                       # Memory allocation
#SBATCH --time=24:00:00                 # Time limit (hh:mm:ss)

# Load the Mamba module
module load Mamba

# Activate the Mamba environment
mamba activate my_fgbio_env

# Confirm the environment is active (optional)
echo "Activated environment: $(conda info --envs | grep '*' | awk '{print $1}')"

# Define directories and files
ubam_folder="/home/janzules/ctDNA_11042024/data/consensusPipeline/intermediate_files/ubam"
umi_corrected_folder="/home/janzules/ctDNA_11042024/data/consensusPipeline/intermediate_files/umi_corrected"
metrics_folder="/home/janzules/ctDNA_11042024/data/consensusPipeline/output/umi_metrics"
umi_file="/path/to/file-of-expected-umis-one-per-line.txt"

# Create output directories if they do not exist
mkdir -p "$umi_corrected_folder"
mkdir -p "$metrics_folder"

# Iterate through each uBAM file
for ubam_file in "$ubam_folder"/*.unmapped.bam; do
    # Extract sample name from file
    sample_name=$(basename "$ubam_file" .unmapped.bam)
    
    # Define output files
    corrected_ubam="$umi_corrected_folder/${sample_name}.corrected.unmapped.bam"
    metrics_file="$metrics_folder/${sample_name}_umi_metrics.txt"

    echo "Processing UMI correction for sample: $sample_name"

    fgbio -Xmx8g --compression 1 --async-io CorrectUmis \
        --input "$ubam_file" \
        --output "$corrected_ubam" \
        --max-mismatches 1 \
        --min-distance 2 \
        --metrics "$metrics_file"

    echo "UMI correction completed for sample: $sample_name"
done

echo "All samples processed."

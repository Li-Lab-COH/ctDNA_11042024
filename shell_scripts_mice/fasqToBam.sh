#!/bin/bash
#SBATCH --job-name=fgbio_job           # Job name
#SBATCH --output=/home/janzules/ctDNA_11042024/code/slurmOutput/fastqToBAM.out      # Standard output log (%j will be replaced with the job ID)
#SBATCH --error=/home/janzules/ctDNA_11042024/code/slurmOutput/fastqToBAM.err       # Standard error log
#SBATCH --ntasks=1                     # Number of tasks (processes)
#SBATCH --cpus-per-task=4              # Number of CPU cores per task
#SBATCH --mem=8G                      # Memory allocation
#SBATCH --time=12:00:00                # Time limit (hh:mm:ss)

# Load the Mamba module
module load Mamba

# Activate the Mamba environment
mamba activate my_fgbio_env

# Confirm the environment is active (optional, for debugging)
echo "Activated environment: $(conda info --envs | grep '*' | awk '{print $1}')"

# Define directories and files
file_locations="/home/janzules/ctDNA_11042024/data/fastp_cleaned/cleaned"
#anumber_file="/home/janzules/ctDNA_11042024/code/addresses/human_Anumbers.txt"
anumber_file="/home/janzules/ctDNA_11042024/code/addresses/Step1_leftOver_human.txt"
ubam_loc="/home/janzules/ctDNA_11042024/data/consensusPipeline/intermediate_files/ubam"
logs_dir="/home/janzules/ctDNA_11042024/data/consensusPipeline/output/logs"

# anumber_report_dir="$fasp_reports_base/$anumber"

while read -r anumber; do
    # Define input files (R1 and R2 for each sample)
    R1_strand="${file_locations}/${anumber}_R1.fastq.gz"
    R2_strand="${file_locations}/${anumber}_R2.fastq.gz"

    if [ -f "$R1_strand" ] && [ -f "$R2_strand" ]; then
        echo "Processing: $anumber"
        
        ubam_output="$ubam_loc/${anumber}.unmapped.bam"
        fgbio -Xmx8g --compression 1 --async-io FastqToBam \
            --input "$R1_strand" "$R2_strand" \
            --read-structures 5M2S+T 5M2S+T \
            --sample "$anumber" \
            --library "$anumber" \
            --output "$ubam_output" &> "$logs_dir/fastq_to_ubam.${anumber}.log"
        
        echo "Completed processing for $anumber"
    else
        echo "Warning: Missing files for $anumber. Skipping."
    fi
done < "$anumber_file"

echo "fastq to BAM complete."

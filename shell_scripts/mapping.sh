#!/bin/bash
#SBATCH --job-name=bwa_alignment        # Job name
#SBATCH --output=/home/janzules/ctDNA_11042024/code/slurmOutput/bwaAlignment.out  # Standard output log
#SBATCH --error=/home/janzules/ctDNA_11042024/code/slurmOutput/bwaAlignment.err   # Standard error log
#SBATCH --ntasks=1                     # Number of tasks (processes)
#SBATCH --cpus-per-task=16             # Number of CPU cores per task
#SBATCH --mem=32G                      # Memory allocation
#SBATCH --time=24:00:00                # Time limit (hh:mm:ss)

# Load necessary modules
module load Mamba
module load SAMtools/1.9-foss-2018a
module load BWA/0.7.17-foss-2018b

# Activate the Mamba environment
mamba activate my_fgbio_env

# Confirm the environment is active (optional, for debugging)
echo "Activated environment: $(conda info --envs | grep '*' | awk '{print $1}')"

# Define directories and files
ubam_dir="/home/janzules/ctDNA_11042024/data/consensusPipeline/intermediate_files/ubam"
mapped_bam_dir="/home/janzules/ctDNA_11042024/data/consensusPipeline/intermediate_files/mapped_bam"
logs_dir="/home/janzules/ctDNA_11042024/data/consensusPipeline/output/logs/mapping_logs"
reference="/home/janzules/bwa_index/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
anumber_file="/home/janzules/ctDNA_11042024/code/addresses/human_Anumbers.txt"

# Make sure output directories exist
# mkdir -p "$mapped_bam_dir" "$logs_dir"

# Process each anumber
while read -r anumber; do
    ubam_file="$ubam_dir/${anumber}.unmapped.bam"
    mapped_bam_file="$mapped_bam_dir/${anumber}.mapped.bam"
    log_file="$logs_dir/bwa_alignment.${anumber}.log"

    if [ -f "$ubam_file" ]; then
        echo "Processing: $anumber"
        
        # Run the alignment and mapping pipeline
        samtools fastq "$ubam_file" | \
        bwa mem -t 16 -p -K 150000000 -Y "$reference" - | \
        fgbio -Xmx28g --compression 1 --async-io ZipperBams \
            --unmapped "$ubam_file" \
            --ref "$reference" \
            --output "$mapped_bam_file" &> "$log_file"
        
        if [ $? -eq 0 ]; then
            echo "Alignment and mapping completed for: $anumber"
        else
            echo "Error during processing for: $anumber. Check $log_file for details."
        fi
    else
        echo "Warning: Missing uBAM file for $anumber. Skipping."
    fi
done < "$anumber_file"

echo "Alignment pipeline complete."

#!/bin/bash
#SBATCH --job-name=fgbio_job           # Job name
#SBATCH --output=fgbio_job.%j.out      # Standard output log (%j will be replaced with the job ID)
#SBATCH --error=fgbio_job.%j.err       # Standard error log
#SBATCH --ntasks=1                     # Number of tasks (processes)
#SBATCH --cpus-per-task=4              # Number of CPU cores per task
#SBATCH --mem=16G                      # Memory allocation
#SBATCH --time=02:00:00                # Time limit (hh:mm:ss)

# Load the Mamba module
module load Mamba

# Activate the Mamba environment
source activate my_fgbio_env

# Confirm the environment is active (optional, for debugging)
echo "Activated environment: $(conda info --envs | grep '*' | awk '{print $1}')"


#!/bin/bash

# Define directories and files
file_locations="/home/janzules/ctDNA_11042024/data/fastp_cleaned/cleaned"
anumber_file="/home/janzules/ctDNA_11042024/code/addresses/Anumbers.txt"
genome="/path/to/hs38DH.fa"  # Replace with actual genome file path

# Create log directory
log_dir="./logs"
mkdir -p "$log_dir"

# Process each Anumber
while read -r anumber; do
    # Define input files (R1 and R2 for each sample)
    R1_strand="${file_locations}/${anumber}_R1.fastq.gz"
    R2_strand="${file_locations}/${anumber}_R2.fastq.gz"

    if [ -f "$R1_strand" ] && [ -f "$R2_strand" ]; then
        echo "Processing: $anumber"

        # 1. Generate uBAM from R1 and R2
        ubam_output="${anumber}.unmapped.bam"
        echo "Starting command: Generate uBAM for $anumber"
        fgbio -Xmx1g --compression 1 --async-io FastqToBam \
            --input "$R1_strand" "$R2_strand" \
            --read-structures "5M2S+T" "5M2S+T" \
            --sample "$anumber" \
            --output "$ubam_output" &> "$log_dir/fastq_to_ubam.${anumber}.log"
        echo "Finished command: Generate uBAM for $anumber"

        # 2. Align BAM using bwa and ZipperBams
        mapped_bam_output="${anumber}.mapped.bam"
        echo "Starting command: Align BAM for $anumber"
        samtools fastq "$ubam_output" \
            | bwa mem -t 16 -p -K 150000000 -Y "$genome" - \
            | fgbio -Xmx4g --compression 1 --async-io ZipperBams \
                --unmapped "$ubam_output" \
                --ref "$genome" \
                --output "$mapped_bam_output" &> "$log_dir/align_bam.${anumber}.log"
        echo "Finished command: Align BAM for $anumber"

        # 3. Group reads by UMI and position
        grouped_bam_output="${anumber}.grouped.bam"
        grouped_stats_output="${anumber}.grouped-family-sizes.txt"
        echo "Starting command: Group reads for $anumber"
        fgbio -Xmx8g --compression 1 --async-io GroupReadsByUmi \
            --input "$mapped_bam_output" \
            --strategy Adjacency \
            --edits 1 \
            --output "$grouped_bam_output" \
            --family-size-histogram "$grouped_stats_output" &> "$log_dir/group_reads.${anumber}.log"
        echo "Finished command: Group reads for $anumber"

        # 4. Call consensus reads
        consensus_ubam_output="${anumber}.cons.unmapped.bam"
        echo "Starting command: Call consensus reads for $anumber"
        fgbio -Xmx4g --compression 1 CallMolecularConsensusReads \
            --input "$grouped_bam_output" \
            --output "$consensus_ubam_output" \
            --min-reads 1 \
            --min-input-base-quality 20 \
            --threads 4 &> "$log_dir/call_consensus_reads.${anumber}.log"
        echo "Finished command: Call consensus reads for $anumber"

        # 5. Align consensus reads to genome
        consensus_mapped_bam_output="${anumber}.cons.mapped.bam"
        echo "Starting command: Align consensus reads for $anumber"
        samtools fastq "$consensus_ubam_output" \
            | bwa mem -t 16 -p -K 150000000 -Y "$genome" - \
            | fgbio -Xmx4g --compression 1 --async-io ZipperBams \
                --unmapped "$consensus_ubam_output" \
                --ref "$genome" \
                --tags-to-reverse Consensus \
                --tags-to-revcomp Consensus \
                --output "$consensus_mapped_bam_output" &> "$log_dir/align_consensus.${anumber}.log"
        echo "Finished command: Align consensus reads for $anumber"

        # 6. Filter and sort consensus reads
        filtered_bam_output="${anumber}.cons.filtered.bam"
        echo "Starting command: Filter and sort consensus reads for $anumber"
        fgbio -Xmx8g --compression 0 FilterConsensusReads \
            --input "$consensus_mapped_bam_output" \
            --output /dev/stdout \
            --ref "$genome" \
            --min-reads 3 \
            --min-base-quality 40 \
            --max-base-error-rate 0.2 \
            | samtools sort --threads 8 -o "$filtered_bam_output"##idx##"$filtered_bam_output.bai" --write-index &> "$log_dir/filter_consensus.${anumber}.log"
        echo "Finished command: Filter and sort consensus reads for $anumber"

        echo "Completed processing for $anumber"
    else
        echo "Warning: Missing files for $anumber. Skipping."
    fi
done < "$anumber_file"

echo "Pipeline completed."


while read -r anumber; do
    # Define input files (R1 and R2 for each sample)
    R1_strand="${file_locations}/${anumber}_R1.fastq.gz"
    R2_strand="${file_locations}/${anumber}_R2.fastq.gz"

    if [ -f "$R1_strand" ] && [ -f "$R2_strand" ]; then
        echo "Processing: $anumber"

        
        echo "Completed processing for $anumber"
    else
        echo "Warning: Missing files for $anumber. Skipping."
    fi
done < "$anumber_file"

echo "Pipeline completed."
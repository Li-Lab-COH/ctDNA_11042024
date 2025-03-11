#!/bin/bash
#SBATCH --job-name=sort_and_qualimap          # Job name
#SBATCH --output=/home/janzules/ctDNA_11042024/slurmOutput/sort_and_qualimap_%A_%a.out
#SBATCH --error=/home/janzules/ctDNA_11042024/slurmOutput/sort_and_qualimap_%A_%a.err
#SBATCH --ntasks=1                            # Number of tasks (processes)
#SBATCH --cpus-per-task=16                    # CPU cores per task (up to you; 16 recommended)
#SBATCH --mem=64G                            # Memory allocation (adjust as needed)
#SBATCH --time=12:00:00                       # Time limit (hh:mm:ss)
#SBATCH --array=0-7                           # One job per sample, adjust range for your total samples

# -----------------------------
#     LOAD MODULES
# -----------------------------
module load samtools
module load Qualimap

# -----------------------------
#    PATHS AND DIRECTORIES
# -----------------------------

# MOUSE
mapped_bam_dir="/home/janzules/ctDNA_11042024/data/output_mice/consensusPipeline/intermediate_files/mapped_bam" 
sorted_bam_dir="/home/janzules/ctDNA_11042024/data/output_mice/consensusPipeline/intermediate_files/sorted_bam"
qualimap_dir="/home/janzules/ctDNA_11042024/data/output_mice/consensusPipeline/output/alignment_metrics/qualimap"
anumber_file="/home/janzules/ctDNA_11042024/code/addresses/mice_Anumbers.txt"


# HUMAN
# mapped_bam_dir="/home/janzules/ctDNA_11042024/data/consensusPipeline/intermediate_files/mapped_bam"
# sorted_bam_dir="/home/janzules/ctDNA_11042024/data/consensusPipeline/intermediate_files/sorted_bam"
# qualimap_dir="/home/janzules/ctDNA_11042024/data/consensusPipeline/output/alignment_metrics/qualimap"
# anumber_file="/home/janzules/ctDNA_11042024/code/addresses/human_Anumbers.txt"

# -----------------------------
#       SETUP OUTPUT DIRS
# -----------------------------
mkdir -p "$sorted_bam_dir"
mkdir -p "$qualimap_dir"

# -----------------------------
#      READ SAMPLE LIST
# -----------------------------
mapfile -t anumbers < "$anumber_file"

# The current sample index is given by the array task ID
ANUMBER_INDEX=$SLURM_ARRAY_TASK_ID

# Basic error checking if index is out of range
if [[ $ANUMBER_INDEX -ge ${#anumbers[@]} ]]; then
    echo "[ERROR] SLURM_ARRAY_TASK_ID ($ANUMBER_INDEX) exceeds number of samples (${#anumbers[@]})."
    exit 1
fi

# Retrieve the A-number
ANUMBER="${anumbers[$ANUMBER_INDEX]}"

# -----------------------------
# 1) SORT THE BAM FILE
# -----------------------------
mapped_bam_file="${mapped_bam_dir}/${ANUMBER}.mapped.bam"
sorted_bam_file="${sorted_bam_dir}/${ANUMBER}.sorted.bam"

if [[ -f "$mapped_bam_file" ]]; then
    echo "=== Sorting BAM for: $ANUMBER ==="
    samtools sort \
        -@ 8 \
        -m 5G \
        -o "$sorted_bam_file" \
        "$mapped_bam_file"

    if [[ $? -eq 0 ]]; then
        echo "[Samtools sort] Completed for $ANUMBER"
    else
        echo "[Samtools sort] ERROR for $ANUMBER"
        exit 1
    fi
else
    echo "[Warning] Missing mapped BAM file for $ANUMBER. Skipping."
    exit 0
fi

# -----------------------------
# 2) RUN QUALIMAP
# -----------------------------
if [[ -f "$sorted_bam_file" ]]; then
    echo "=== Running Qualimap for: $ANUMBER ==="
    qualimap bamqc \
        -bam "$sorted_bam_file" \
        -outdir "${qualimap_dir}/${ANUMBER}_qcReport" \
        -nt 16 \
        --java-mem-size=64G

    if [[ $? -eq 0 ]]; then
        echo "[Qualimap] Completed for $ANUMBER"
    else
        echo "[Qualimap] ERROR for $ANUMBER - Check log"
    fi
else
    echo "[Warning] Sorted BAM not found for $ANUMBER. Skipping Qualimap."
fi

echo "=== All done for $ANUMBER ==="

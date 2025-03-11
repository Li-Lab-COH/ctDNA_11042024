#!/bin/bash
#SBATCH --job-name=fastqToBam_BWA_pipeline    # Job name
#SBATCH --output=/home/janzules/ctDNA_11042024/slurmOutput/ubam_n_mapping.out
#SBATCH --error=/home/janzules/ctDNA_11042024/slurmOutput/ubam_n_mapping.err
#SBATCH --ntasks=1                            # Number of tasks (processes)
#SBATCH --cpus-per-task=16                    # Number of CPU cores per task
#SBATCH --mem=48G                             # Memory allocation
#SBATCH --time=24:00:00                       # Time limit (hh:mm:ss)

# -----------------------------
#      USER CONFIGURATION
# -----------------------------

# MICE

READS_DIR="/home/janzules/ctDNA_11042024/data/output_mice/fastp_cleaned/cleaned"  #
ANUMBER_FILE="/home/janzules/ctDNA_11042024/code/addresses/mice_Anumbers.txt" #

UBAM_DIR="/home/janzules/ctDNA_11042024/data/output_mice/consensusPipeline/intermediate_files/ubam" #
MAPPED_BAM_DIR="/home/janzules/ctDNA_11042024/data/output_mice/consensusPipeline/intermediate_files/mapped_bam" #
LOGS_DIR="/home/janzules/ctDNA_11042024/data/output_mice/consensusPipeline/output/logs/ubam_n_mappings" #

# Reference genome used for alignment
REFERENCE="/home/janzules/reference_genomes/mouse_GRCh39/GRCm39.primary_assembly.genome.fa" #

# Number of threads for alignment
THREADS=16



# HUMAN
# READS_DIR="/home/janzules/ctDNA_11042024/data/fastp_cleaned/cleaned"
# ANUMBER_FILE="/home/janzules/ctDNA_11042024/code/addresses/Step1_leftOver_human.txt"

# UBAM_DIR="/home/janzules/ctDNA_11042024/data/consensusPipeline/intermediate_files/ubam"
# MAPPED_BAM_DIR="/home/janzules/ctDNA_11042024/data/consensusPipeline/intermediate_files/mapped_bam"
# LOGS_DIR="/home/janzules/ctDNA_11042024/data/consensusPipeline/output/logs"

# # Reference genome used for alignment
# REFERENCE="/home/janzules/bwa_index/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

# # Number of threads for alignment
# THREADS=16

# -----------------------------
#    LOAD MODULES / ENV
# -----------------------------

module load Mamba
# SAMtools and BWA modules (used in the mapping step)
module load SAMtools/1.9-foss-2018a
module load BWA/0.7.17-foss-2018b

# Activate the Mamba environment
mamba activate my_fgbio_env

echo "Activated environment: $(conda info --envs | grep '*' | awk '{print $1}')"

# Make sure output directories exist (no errors if they already exist)
mkdir -p "${UBAM_DIR}" "${MAPPED_BAM_DIR}" "${LOGS_DIR}"

# -----------------------------
#       MAIN PIPELINE
# -----------------------------
while read -r anumber; do
    
    # Define the expected reads
    R1="${READS_DIR}/${anumber}_R1.fastq.gz"
    R2="${READS_DIR}/${anumber}_R2.fastq.gz"
    UBAM_OUTPUT="${UBAM_DIR}/${anumber}.unmapped.bam"
    MAPPED_BAM="${MAPPED_BAM_DIR}/${anumber}.mapped.bam"
    
    LOG_FASTQ2UBAM="${LOGS_DIR}/fastq_to_ubam.${anumber}.log"
    LOG_BWA="${LOGS_DIR}/bwa_alignment.${anumber}.log"

    # -----------------------------
    #   FASTQ -> unmapped BAM
    # -----------------------------
    if [[ -f "$R1" && -f "$R2" ]]; then
        echo "=== [FastqToBam] Processing: ${anumber} ==="
        fgbio -Xmx8g --compression 1 --async-io FastqToBam \
            --input "$R1" "$R2" \
            --read-structures "5M2S+T" "5M2S+T" \
            --sample "$anumber" \
            --library "$anumber" \
            --output "$UBAM_OUTPUT" &> "$LOG_FASTQ2UBAM"

        if [[ $? -eq 0 ]]; then
            echo "[FastqToBam] Completed for ${anumber}"
        else
            echo "[FastqToBam] ERROR for ${anumber} - see log: ${LOG_FASTQ2UBAM}"
            continue  # Move to the next sample
        fi
    else
        echo "[FastqToBam] Warning: Missing files for $anumber (R1 or R2). Skipping."
        continue
    fi

    # -----------------------------
    #        MAPPING STEP
    # -----------------------------
    if [[ -f "$UBAM_OUTPUT" ]]; then
        echo "=== [Mapping] Processing: ${anumber} ==="
        samtools fastq "$UBAM_OUTPUT" | \
        bwa mem -t "$THREADS" -p -K 150000000 -Y "$REFERENCE" - | \
        fgbio -Xmx28g --compression 1 --async-io ZipperBams \
            --unmapped "$UBAM_OUTPUT" \
            --ref "$REFERENCE" \
            --output "$MAPPED_BAM" &> "$LOG_BWA"

        if [[ $? -eq 0 ]]; then
            echo "[Mapping] Completed for ${anumber}"
        else
            echo "[Mapping] ERROR for ${anumber} - see log: ${LOG_BWA}"
        fi
    else
        echo "[Mapping] Warning: Missing uBAM file for $anumber. Skipping."
        continue
    fi

done < "$ANUMBER_FILE"

echo "=== Pipeline complete. ==="

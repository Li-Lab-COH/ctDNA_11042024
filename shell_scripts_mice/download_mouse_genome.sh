#!/bin/bash
#SBATCH --job-name=download_gencode_mouse
#SBATCH --output=/home/janzules/ctDNA_11042024/code/shell_scripts_mice/slurmOutput/download_gencode_mouse.log
#SBATCH --error=/home/janzules/ctDNA_11042024/code/shell_scripts_mice/slurmOutput/download_gencode_mouse.err
#SBATCH --time=01:00:00  # Adjust time as needed
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G

# Define variables
GENCODE_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/latest_release/"
DEST_DIR="/home/janzules/bwa_index/mouse_GRCh39/"  # Adjust to your preferred directory

# Ensure destination directory exists
mkdir -p "$DEST_DIR"

# Navigate to the directory
cd "$DEST_DIR"

# Download the primary assembly genome
wget -O GRCm39.primary_assembly.genome.fa.gz "$GENCODE_URL"

# Verify successful download
if [ -f "GRCm39.primary_assembly.genome.fa.gz" ]; then
    echo "Download completed successfully."
else
    echo "Download failed!"
    exit 1
fi

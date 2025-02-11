#!/bin/bash
#SBATCH --job-name=bwa_index
#SBATCH --output=slurmOutput/bwa_index.out      # Output file name (%j will be replaced with job ID)
#SBATCH --error=slurmOutput/bwa_index.err       # Error file name
#SBATCH --ntasks=1                     # Number of tasks (usually 1 for BWA indexing)
#SBATCH --cpus-per-task=8              # Number of CPU cores
#SBATCH --mem=16G                      # Memory requested
#SBATCH --time=04:00:00                # Time limit (adjust as necessary)


# Load BWA module if available (adjust according to your environment)
module load BWA/0.7.17-foss-2018b

# Path to your reference genome file
GENOME_PATH="/home/janzules/bwa_index/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

# Run BWA indexing
bwa index $GENOME_PATH

echo "BWA indexing completed for $GENOME_PATH"

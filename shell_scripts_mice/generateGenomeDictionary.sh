#!/bin/bash -l
#SBATCH --job-name=create_dict            # Job name
#SBATCH --output=/home/janzules/ctDNA_11042024/code/slurmOutput/create_dict.out  # Standard output
#SBATCH --error=/home/janzules/ctDNA_11042024/code/slurmOutput/create_dict.err   # Standard error
#SBATCH --ntasks=1                        # Number of tasks
#SBATCH --cpus-per-task=8                 # Number of CPU cores
#SBATCH --mem=16G                          # Memory allocation
#SBATCH --time=01:00:00                   # Time limit (hh:mm:ss)

# Load Picard
module load picard

# Define reference genome and output dictionary paths
REFERENCE="/home/janzules/bwa_index/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
OUTPUT="/home/janzules/bwa_index/Homo_sapiens.GRCh38.dna.primary_assembly.dict"

# Create the sequence dictionary
java -jar /opt/picard/2.21.1/picard.jar CreateSequenceDictionary \
    R="$REFERENCE" \
    O="$OUTPUT"

# Confirm successful dictionary creation
if [ -f "$OUTPUT" ]; then
    echo "Sequence dictionary created: $OUTPUT"
else
    echo "Error: Sequence dictionary creation failed!"
fi

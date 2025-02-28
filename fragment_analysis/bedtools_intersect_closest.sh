#!/bin/bash
#SBATCH --job-name=process_bam
#SBATCH --output=logs/%x_%A_%a.out  # Unique log files per job
#SBATCH --error=logs/%x_%A_%a.err   # Unique error log files per job
#SBATCH --array=0-7                 # 8 jobs for 16 files (2 per job)
#SBATCH --time=01:00:00             # Adjust as needed
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G                     # Adjust memory as needed

# Load required modules
module load bedtools  # Ensure bedtools is available

# Ensure the FOLDER_PATH variable is set
if [ -z "$FOLDER_PATH" ]; then
    echo "Error: FOLDER_PATH is not set. Exiting."
    exit 1
fi

# Defining locations
BAM_FOLDER="$FOLDER_PATH/bam"
OUTPUT_DIR="$FOLDER_PATH/bedtools_out/geneIntersect"
GFF3_File="/home/janzules/reference_genomes/Human_GRCh38_ensembl/annotation_files/Homo_sapiens.GRCh38.113.chr.gff3"

# Create the output directories if they don't exist
mkdir -p "$OUTPUT_DIR"

# Get the list of BAM files
BAM_FILES=($(ls "$BAM_FOLDER"/*.bam))

# Calculate the range of BAM files to process in this array job
START_INDEX=$(( SLURM_ARRAY_TASK_ID * 2 ))
END_INDEX=$(( START_INDEX + 1 ))

# Ensure the end index does not exceed the total number of BAM files
if (( END_INDEX >= ${#BAM_FILES[@]} )); then
    END_INDEX=$(( ${#BAM_FILES[@]} - 1 ))
fi

# Loop over the two BAM files assigned to this task
for (( i=START_INDEX; i<=END_INDEX; i++ )); do
    BAM_FILE="${BAM_FILES[$i]}"
    
    # Ensure BAM file is valid
    if [[ -f "$BAM_FILE" ]]; then
        echo "Processing: $BAM_FILE"
        
        # Extract the filename without the path and extension
        BAM_BASENAME=$(basename "$BAM_FILE" .bam)

        # Define the output file path
        OUTPUT_FILE="$OUTPUT_DIR/${BAM_BASENAME}_geneIntersect.bed"

        # Run bedtools intersect
        bedtools intersect -a "$BAM_FILE" -b "$GFF3_File" -wa -wb -bed > "$OUTPUT_FILE"

        echo "Output saved to: $OUTPUT_FILE"
    else
        echo "Warning: BAM file not found - $BAM_FILE"
    fi
done


bedtools intersect -a 

# For each bin + sample, you’ll produce:
# A14891_40_150_genesIntersect.bed
# A14891_40_150_nucIntersect.bed
# A14891_40_150_nucClosest.bed
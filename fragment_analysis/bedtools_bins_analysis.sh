#!/bin/bash
#SBATCH --job-name=BinBedtools
#SBATCH --output=/home/janzules/ctDNA_11042024/slurmOutput/binningAnalysis/bedtools_JA/bin/%x_%A.out  
#SBATCH --error=/home/janzules/ctDNA_11042024/slurmOutput/binningAnalysis/bedtools_JA/bin/%x_%A.err   
#SBATCH --array=0-7                 # 8 jobs for 16 files (2 per job)
#SBATCH --time=012:00:00             
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --mail-user=janzules@coh.org
#SBATCH --mail-type=FAIL

# Load required modules
module load bedtools  # Ensure bedtools is available

# Ensure the FOLDER_PATH variable is set
# TODO: change this for the nonbinned files
if [ -z "$FOLDER_PATH" ]; then
    echo "[ERROR] FOLDER_PATH is not set. Exiting." >&2
    exit 1
fi

# Define locations
# TODO: change this for the nonbinned files
BED_FOLDER="$FOLDER_PATH/bed"


OUTPUT_DIR_GENE="$FOLDER_PATH/bedtools_out/geneIntersect"
OUTPUT_DIR_NUCL_INTERSECT="$FOLDER_PATH/bedtools_out/nuclIntersect"
# OUTPUT_DIR_NUCL_CLOSEST="$FOLDER_PATH/bedtools_out/nuclClosest"

# This is removed because we don't to find matches to EVERYTHING, but keeping it here as a record for the change
# GFF3_FILE="/home/janzules/reference_genomes/Human_GRCh38_ensembl/annotation_files/Homo_sapiens.GRCh38.113.chr.gff3"
GFF3_FILE="/home/janzules/reference_genomes/Human_GRCh38_ensembl/annotation_files/genes_only.gff3"
NUC_FILE="/home/janzules/reference_genomes/Human_GRCh38_ensembl/nucleosomes/GSE71378_nuc_with_IDs.bed"

# Create output directories if they don't exist
mkdir -p "$OUTPUT_DIR_GENE"
mkdir -p "$OUTPUT_DIR_NUCL_INTERSECT"
# mkdir -p "$OUTPUT_DIR_NUCL_CLOSEST"


# Get the list of BED files (sorted for consistency)
BED_FILES=($(ls "$BED_FOLDER"/*.bed | sort -V))

# Ensure there are BED files
if [[ ${#BED_FILES[@]} -eq 0 ]]; then
    echo "[ERROR] No BED files found in $BED_FOLDER. Exiting." >&2
    exit 1
fi

# Calculate the range of BED files to process in this array job
START_INDEX=$(( SLURM_ARRAY_TASK_ID * 2 ))
END_INDEX=$(( START_INDEX + 1 ))

# Ensure the end index does not exceed the total number of BED files
if (( END_INDEX >= ${#BED_FILES[@]} )); then
    END_INDEX=$(( ${#BED_FILES[@]} - 1 ))
fi

# Loop over the two BED files assigned to this task
for (( i=START_INDEX; i<=END_INDEX; i++ )); do
    BED_FILE="${BED_FILES[$i]}"

    # Ensure BED file is valid
    if [[ -f "$BED_FILE" ]]; then
        echo "[INFO] Processing: $BED_FILE"
        
        # Extract filename without the path and extension
        BED_BASENAME=$(basename "$BED_FILE" .bed)

        # Define output file paths
        GENES_OUTPUT="$OUTPUT_DIR_GENE/${BED_BASENAME}_genesIntersect.bed"
        NUC_OUTPUT="$OUTPUT_DIR_NUCL_INTERSECT/${BED_BASENAME}_nucIntersect.bed"
        NUC_CLOSEST_OUTPUT="$OUTPUT_DIR_NUCL_CLOSEST/${BED_BASENAME}_nucClosest.bed"

        # Intersect BED with GFF3 (genes)
        bedtools intersect -a "$GFF3_FILE" -b "$BED_FILE" -c > "$GENES_OUTPUT"
        # bedtools intersect -a "$BED_FILE" -b "$GFF3_FILE" -wa -wb > "$GENES_OUTPUT"
        echo "[INFO] Gene intersection saved to: $GENES_OUTPUT"

        # Intersect BED with nucleosome BED
        # bedtools intersect -a "$BED_FILE" -b "$NUC_FILE" -wa -wb > "$NUC_OUTPUT"
        bedtools intersect -a "$NUC_FILE" -b "$BED_FILE" -c > "$NUC_OUTPUT"
        echo "[INFO] Nucleosome intersection saved to: $NUC_OUTPUT"

        # Find closest nucleosome (excluding overlaps, reporting distances)
        # bedtools closest -a "$BED_FILE" -b "$NUC_FILE" -D ref -io > "$NUC_CLOSEST_OUTPUT"
        # echo "[INFO] Closest nucleosome positions saved to: $NUC_CLOSEST_OUTPUT"

    else
        echo "[WARNING] BED file not found - $BED_FILE" >&2
    fi
done

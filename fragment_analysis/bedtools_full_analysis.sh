#!/bin/bash
#SBATCH --job-name=BedtoolsFull
#SBATCH --output=/home/janzules/ctDNA_11042024/slurmOutput/binningAnalysis/bedtools_JA/full/%x_%A.out  
#SBATCH --error=/home/janzules/ctDNA_11042024/slurmOutput/binningAnalysis/bedtools_JA/full/%x_%A.err   
#SBATCH --array=0-15  # 16 jobs for 16 BAM files (1 per job)
#SBATCH --time=12:00:00             
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --mail-user=janzules@coh.org
#SBATCH --mail-type=FAIL

# Load required modules
module load bedtools  # Ensure bedtools is available

# Ensure the FOLDER_PATH variable is set
FOLDER_PATH="/home/janzules/ctDNA_11042024/data/consensusPipeline/intermediate_files/"
if [ -z "$FOLDER_PATH" ]; then
    echo "[ERROR] FOLDER_PATH is not set. Exiting." >&2
    exit 1
fi

# Define locations
BAM_FOLDER="$FOLDER_PATH/sorted_bam"
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
mkdir -p "$OUTPUT_DIR_NUCL_CLOSEST"

# Get the list of BAM files (sorted for consistency)
BAM_FILES=($(ls "$BAM_FOLDER"/*.sorted.bam | sort -V))

# Ensure there are BAM files
if [[ ${#BAM_FILES[@]} -eq 0 ]]; then
    echo "[ERROR] No BAM files found in $BAM_FOLDER. Exiting." >&2
    exit 1
fi

# Get the BAM file corresponding to this job array index
BAM_FILE="${BAM_FILES[$SLURM_ARRAY_TASK_ID]}"

# Ensure BAM file is valid
if [[ -f "$BAM_FILE" ]]; then
    echo "[INFO] Processing: $BAM_FILE"
    
    # Extract filename without the path and extension
    BAM_BASENAME=$(basename "$BAM_FILE" .sorted.bam)

    # Define output file paths
    GENES_OUTPUT="$OUTPUT_DIR_GENE/${BAM_BASENAME}_genesIntersect.bed"
    NUC_OUTPUT="$OUTPUT_DIR_NUCL_INTERSECT/${BAM_BASENAME}_nucIntersect.bed"
    NUC_CLOSEST_OUTPUT="$OUTPUT_DIR_NUCL_CLOSEST/${BAM_BASENAME}_nucClosest.bed"

    # Intersect BAM with GFF3 (genes)
    # bedtools intersect -abam "$BAM_FILE" -b "$GFF3_FILE" -wa -wb -bed > "$GENES_OUTPUT" # this generated terrabytes of data (too much man)
    bedtools intersect -a "$GFF3_FILE" -b "$BAM_FILE" -c > "$GENES_OUTPUT"
    echo "[INFO] Gene intersection saved to: $GENES_OUTPUT"
    
    # Intersect BAM with nucleosome BED
    # bedtools intersect -abam "$BAM_FILE" -b "$NUC_FILE" -wa -wb -bed > "$NUC_OUTPUT"
    bedtools intersect -a "$NUC_FILE" -b "$BAM_FILE" -c > "$NUC_OUTPUT"
    echo "[INFO] Nucleosome intersection saved to: $NUC_OUTPUT"

    # echo "[INFO] Beginning Nucleosome closest for $BAM_FILE "
    # # Find closest nucleosome (excluding overlaps, reporting distances)
    # bedtools bamtobed -i "$BAM_FILE" | \
    #     bedtools closest -a stdin -b "$NUC_FILE" -D ref -io > "$NUC_CLOSEST_OUTPUT"
    # echo "[INFO] Closest nucleosome positions saved to: $NUC_CLOSEST_OUTPUT"

else
    echo "[WARNING] BAM file not found - $BAM_FILE" >&2
fi

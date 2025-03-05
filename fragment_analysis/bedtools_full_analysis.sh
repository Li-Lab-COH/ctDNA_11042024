#!/bin/bash
#SBATCH --job-name=BedtoolsFull
#SBATCH --output=/home/janzules/ctDNA_11042024/slurmOutput/binningAnalysis/bedtools_JA/full/%x_%A.out  
#SBATCH --error=/home/janzules/ctDNA_11042024/slurmOutput/binningAnalysis/bedtools_JA/full/%x_%A.err   
#SBATCH --array=0-7 #8 jobs
#SBATCH --time=12:00:00             
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=128G
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


# This is removed because we don't need to find matches to EVERYTHING, but keeping it here as a record for the change
# GFF3_FILE="/home/janzules/reference_genomes/Human_GRCh38_ensembl/annotation_files/Homo_sapiens.GRCh38.113.chr.gff3"
GFF3_FILE="/home/janzules/reference_genomes/Human_GRCh38_ensembl/annotation_files/genes_only.gff3"
NUC_FILE="/home/janzules/reference_genomes/Human_GRCh38_ensembl/nucleosomes/GSE71378_nuc_with_IDs.bed"

# Create output directories if they don't exist
mkdir -p "$OUTPUT_DIR_GENE"
mkdir -p "$OUTPUT_DIR_NUCL_INTERSECT"
# mkdir -p "$OUTPUT_DIR_NUCL_CLOSEST"

# Get the list of BAM files (sorted for consistency)
BAM_FILES=($(ls "$BAM_FOLDER"/*.sorted.bam | sort -V))

# Ensure there are BAM files - sanity check
if [[ ${#BAM_FILES[@]} -eq 0 ]]; then
    echo "[ERROR] No BAM files found in $BAM_FOLDER. Exiting." >&2
    exit 1
fi

# Determine which 2 files this array task should process
# Each array index will handle two BAM files.
START_INDEX=$(( SLURM_ARRAY_TASK_ID * 2 ))
END_INDEX=$(( START_INDEX + 1 ))

# Make sure end_index doesn't exceed the array length
if (( END_INDEX >= ${#BAM_FILES[@]} )); then
    END_INDEX=$(( ${#BAM_FILES[@]} - 1 ))
fi

# Get the BAM file corresponding to this job array index
# BAM_FILE="${BAM_FILES[$SLURM_ARRAY_TASK_ID]}


for (( i=START_INDEX; i<=END_INDEX; i++ )); do
    # Safety check
    if (( i >= ${#BAM_FILES[@]} )); then
        break
    fi

    BAM_FILE="${BAM_FILES[$i]}"
    if [[ ! -f "$BAM_FILE" ]]; then
        echo "[WARNING] Could not find $BAM_FILE"
        continue
    fi
    
    BAM_BASENAME=$(basename "$BAM_FILE" .sorted.bam)
    echo "[INFO] Processing: $BAM_FILE"

    # Output file paths
    GENES_OUTPUT="$OUTPUT_DIR_GENE/${BAM_BASENAME}_genesIntersect.bed"
    NUC_OUTPUT="$OUTPUT_DIR_NUCL_INTERSECT/${BAM_BASENAME}_nucIntersect.bed"


    # - Count how many alignments overlap each gene
    bedtools bamtobed -i "$BAM_FILE" | bedtools intersect -a "$GFF3_FILE" -b stdin -c > "$GENES_OUTPUT"
    echo "[INFO] Gene intersection saved to: $GENES_OUTPUT"


    # - Count how many alignments overlap each nucleosome region
    bedtools bamtobed -i "$BAM_FILE" | bedtools intersect -a "$NUC_FILE" -b stdin -c > "$NUC_OUTPUT"
    echo "[INFO] Nucleosome intersection saved to: $NUC_OUTPUT"
    
    # echo "[INFO] Beginning Nucleosome closest for $BAM_FILE "
    # # Find closest nucleosome (excluding overlaps, reporting distances)
    # bedtools bamtobed -i "$BAM_FILE" | \
    #     bedtools closest -a stdin -b "$NUC_FILE" -D ref -io > "$NUC_CLOSEST_OUTPUT"
    # echo "[INFO] Closest nucleosome positions saved to: $NUC_CLOSEST_OUTPUT"
done

echo "[INFO] Done with array task $SLURM_ARRAY_TASK_ID"

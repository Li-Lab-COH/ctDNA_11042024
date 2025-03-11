#!/bin/bash
#SBATCH --job-name=/home/janzules/ctDNA_11042024/slurmOutput/extract_rx
#SBATCH --output=/home/janzules/ctDNA_11042024/slurmOutput/rx_counts.out
#SBATCH --error=rx_counts.err
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2

module load samtools  # Load samtools module if needed

# Define input and output files
BAM_FILE="/home/janzules/ctDNA_11042024/data/consensusPipeline/intermediate_files/grouped_bam/A14897.grouped.bam"
RX_COUNT_FILE="/home/janzules/ctDNA_11042024/data/consensusPipeline/intermediate_files/grouped_bam/A14897_rx_counts.txt"
OUTPUT_BAM="/home/janzules/ctDNA_11042024/data/consensusPipeline/intermediate_files/grouped_bam/A14897_top3_rx.bam"

# Step 1: Extract RX values and their counts, sorted
samtools view "$BAM_FILE" | awk '{for(i=12;i<=NF;i++) if($i ~ /^RX:Z:/) print $i}' | sort | uniq -c > "$RX_COUNT_FILE"

# Step 2: Extract the first 3 unique RX values
TOP3_RX=$(awk '{print $2}' "$RX_COUNT_FILE" | head -3)

# Step 3: Filter BAM file for reads matching the first 3 RX values
samtools view -h "$BAM_FILE" | awk -v rx1=$(echo $TOP3_RX | awk '{print $1}') -v rx2=$(echo $TOP3_RX | awk '{print $2}') -v rx3=$(echo $TOP3_RX | awk '{print $3}') '
    BEGIN {OFS="\t"}
    /^@/ {print $0; next} 
    {
        for(i=12; i<=NF; i++) {
            if($i ~ ("RX:Z:" rx1) || $i ~ ("RX:Z:" rx2) || $i ~ ("RX:Z:" rx3)) {
                print $0;
                break;
            }
        }
    }' | samtools view -bS -o "$OUTPUT_BAM"

# Step 4: Sort and index the output BAM file for IGV
samtools sort -o "${OUTPUT_BAM%.bam}.sorted.bam" "$OUTPUT_BAM"
samtools index "${OUTPUT_BAM%.bam}.sorted.bam"

echo "Filtered BAM saved as ${OUTPUT_BAM%.bam}.sorted.bam"

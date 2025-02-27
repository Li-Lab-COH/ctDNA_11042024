#!/usr/bin/env bash

# Absolute path to your GTF file
GTF_FILE="/Users/janzules/Roselab/references/hg38_gtf_n_bed/Homo_sapiens.GRCh38.113.chr.gtf"

# extracts the folder path from the GTF filename.
OUT_DIR="$(dirname "$GTF_FILE")"
OUT_FILE="${OUT_DIR}/Homo_sapiens.GRCh38.113.chr.gene.bed"

echo "[INFO] Converting GTF to BED..."
echo "[INFO] GTF:       $GTF_FILE"
echo "[INFO] Output:    $OUT_FILE"

ggrep -P "\tgene\t" "$GTF_FILE" \
  | cut -f1,4,5,7,9 \
  | sed 's/[[:space:]]/\t/g' \
  | sed 's/[;|"]//g' \
  | awk -F $'\t' 'BEGIN { OFS=FS } {
      # $1 = chromosome
      # $2 = start (1-based in GTF)
      # $3 = end
      # $4 = strand
      # $5 = attribute string (gene_id, gene_version, gene_name, etc.)

      # Convert from 1-based to 0-based start:
      print $1, ($2 - 1), $3, $6, ".", $4, $10, $12, $14
    }' \
  | sort -k1,1 -k2,2n > "$OUT_FILE"

echo "[INFO] Finished creating gene-level BED."
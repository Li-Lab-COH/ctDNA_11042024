#!/bin/bash

# Path to the file containing the list of folders
# FOLDER_LIST="/home/janzules/ctDNA_11042024/code/addresses/mapped.n.bin.bam.files.txt"
FOLDER_LIST="/home/janzules/ctDNA_11042024/code/addresses/fixingCorruption.mapped.n.bin.bam.files.txt"

# Read the folder paths line by line
while read -r FOLDER_PATH; do
    # Submit a job for each folder
    sbatch --export=FOLDER_PATH="$FOLDER_PATH" bedtools_bins_analysis.sh

done < "$FOLDER_LIST"

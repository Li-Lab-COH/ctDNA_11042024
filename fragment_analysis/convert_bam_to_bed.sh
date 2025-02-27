#!/bin/bash
#SBATCH --job-name=convert_bam_to_bed      # Job name
#SBATCH --output=/home/janzules/ctDNA_11042024/slurmOutput/binningAnalysis/convert_bam_to_bed.out
#SBATCH --error=/home/janzules/ctDNA_11042024/slurmOutput/binningAnalysis/convert_bam_to_bed.err
#SBATCH --ntasks=1                         
#SBATCH --cpus-per-task=2                  
#SBATCH --mem=8G                          
#SBATCH --time=4:00:00                     

module load bedtools


# Define the base directory that contains your binned subfolders
BINNED_DIR="/home/janzules/ctDNA_11042024/data/human_binned_sequences"

# TODO: for future analysis with more custom folder, you can take the name of the bin folders in the folder of interest
bins=("40_150" "125_155" "150_180" "170_220" "220_345")


# Loop over each bin folder
for bin_label in "${bins[@]}"; do
    
    # Paths for this bin
    bam_dir="${BINNED_DIR}/${bin_label}/bam"
    bed_dir="${BINNED_DIR}/${bin_label}/bed"

    # Create the bed_dir if it doesn't exist
    mkdir -p "$bed_dir"


    # Convert each .bam in bam_dir to .bed
    for bam_file in "${bam_dir}"/*.bam; do
        
        # check if file exists
        if [ ! -e "$bam_file" ]; then
        echo "[ERROR] No BAM files found in $bam_dir" >&2
        continue
        fi

        # Extract the base filename (e.g. A14891_40_150) without the .bam
        base_name=$(basename "$bam_file" .bam)

        # Construct the output bed path
        bed_file="${bed_dir}/${base_name}.bed"

        echo "[INFO] Converting $bam_file -> $bed_file"

        # Run bedtools
        bedtools bamtobed -i "$bam_file" > "$bed_file"
    done
done

echo "[INFO] BAM to BED conversion complete."

#!/bin/bash
#SBATCH --job-name=alignment_metrics        # Job name
#SBATCH --output=/home/janzules/ctDNA_11042024/code/slurmOutput/alignment_metrics.out  # Standard output log
#SBATCH --error=/home/janzules/ctDNA_11042024/code/slurmOutput/alignment_metrics.err   # Standard error log
#SBATCH --ntasks=1                     # Number of tasks (processes)
#SBATCH --cpus-per-task=16              # Number of CPU cores per task
#SBATCH --mem=188G                      # Memory allocation
#SBATCH --time=12:00:00                # Time limit (hh:mm:ss)

# Loading modules
module load R
module load Qualimap
module load samtools


sorted_bam_dir="/home/janzules/ctDNA_11042024/data/consensusPipeline/intermediate_files/sorted_bam"
# anumber_file="/home/janzules/ctDNA_11042024/code/addresses/Step1_leftOver_human.txt"
anumber_file="/home/janzules/ctDNA_11042024/code/addresses/test_1_human.txt"
output_dir="/home/janzules/ctDNA_11042024/data/consensusPipeline/output/alignment_metrics"
intermediate_dir="/home/janzules/ctDNA_11042024/data/consensusPipeline/intermediate_files"
# name_mapping_file="/home/janzules/ctDNA_11042024/code/addresses/human_samples_name_mapping.csv"

# Additional folders
insert_size_dir="$output_dir/insert_size"
qualmap_dir="$output_dir/qualimap"
dedup_dir="$intermediate_dir/dedup_bam"
dup_dir_metrics="$output_dir/duplicates"





# Create necessary directories if they do not exist"
mkdir -p "$output_dir" "$insert_size_dir" "$qualmap_dir" "$dedup_dir" "$dup_dir_metrics"

# Process each anumber
while read -r anumber; do

    # Customizing inputs and outputs
    sorted_bam_file="$sorted_bam_dir/${anumber}.sorted.bam"
    
    # Insert size (IS) metrics Picard
    IS_txt="$insert_size_dir/${anumber}.insert_size_metrics.txt"
    IS_hist="$insert_size_dir/${anumber}.inser_size_histogram.pdf"
    
    # Quality map
    qualimap_out="$qualmap_dir/${anumber}.qcReport"

    # Duplication report
    dedup_bam_out="$dedup_dir/${anumber}.dedup.bam"
    dup_metrics_out="$dup_dir_metrics/${anumber}.dup_metrics.txt" 

    # Analyzing

    if [ -f "$sorted_bam_file" ]; then
        echo "Collecting insert size for: $anumber"

        # Run Insert Size metrics
        java -jar /opt/picard/2.21.1/picard.jar CollectInsertSizeMetrics \
            I=$sorted_bam_file \
            O=$IS_txt \
            H=$IS_hist \
            M=0.5 # 50% of the reads must have an insert size to include that size in the calculations.

        if [ $? -eq 0 ]; then
            echo "Insert size metrics successful: $anumber"
        else
            echo "Error during Insert size collection"
        fi
        echo "starting quality control collection for: $anumber"
        # Quality Control 
        qualimap bamqc \
            -bam $sorted_bam_file \
            -outdir $qualimap_out \
            -nt 16

        if [ $? -eq 0 ]; then
            echo "Quality map generation successful for: $anumber"
        else
            echo "Error during Quality map generation"
        fi
        echo "Starting Duplication analysis: $anumber"
        # # Duplicate metrics
        # java -jar /opt/picard/2.21.1/picard.jar MarkDuplicates \
        #     I=$sorted_bam_file \
        #     O=$dedup_bam_out \
        #     M=$dup_metrics_out \
        #     REMOVE_DUPLICATES=false

        # if [ $? -eq 0 ]; then
        #     echo "Duplication analysis complete for for: $anumber"
        # else
        #     echo "Error during Duplication analysis"
        # fi   



    else
        echo "Warning: Missing sorted BAM file for $anumber. Skipping."
    fi
done < "$anumber_file"

echo "Metrics collection complete"

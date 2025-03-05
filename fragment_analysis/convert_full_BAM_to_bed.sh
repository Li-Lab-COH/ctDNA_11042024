#!/bin/bash
#SBATCH --job-name=BamToBed      # Job name
#SBATCH --output=/home/janzules/ctDNA_11042024/slurmOutput/binningAnalysis/convert_bam_to_bed.out
#SBATCH --error=/home/janzules/ctDNA_11042024/slurmOutput/binningAnalysis/convert_bam_to_bed.err
#SBATCH --ntasks=1                         
#SBATCH --cpus-per-task=8                  
#SBATCH --mem=32G                          
#SBATCH --time=8:00:00
#SBATCH --array=0-7                      # Job array index (8 jobs total)     

anumber_file="/home/janzules/ctDNA_11042024/code/addresses/human_Anumbers.txt"

# 2) Read the lines of anumber_file into an array
mapfile -t anumbers < "$anumber_file"

# 3) Calculate the range of anumbers to process for this job array task
#    Each task processes two anumbers: indices (start_index) and (end_index).
start_index=$(( SLURM_ARRAY_TASK_ID * 2 ))
end_index=$(( start_index + 1 ))

# Ensure end_index doesn't exceed the array length
if (( end_index >= ${#anumbers[@]} )); then
    end_index=$(( ${#anumbers[@]} - 1 ))
fi


# We nare going to try the intersect of the bam file against the .gtf files of the human genome

# This is the command I'm using:

bedtools intersect -a mapped_reads.bam -b human_annotation.gff3 -bed

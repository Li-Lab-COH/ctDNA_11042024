#!/bin/bash
#SBATCH --job-name=bins_array            # Job name
#SBATCH --output=/home/janzules/ctDNA_11042024/slurmOutput/binningAnalysis/bins_%a.out
#SBATCH --error=/home/janzules/ctDNA_11042024/slurmOutput/binningAnalysis/bins_%a.err
#SBATCH --ntasks=1                       
#SBATCH --cpus-per-task=8               
#SBATCH --mem=32G                       
#SBATCH --time=12:00:00                 
#SBATCH --array=0-7                      # Job array index (8 jobs total)

# 1) File with 16 lines, each line = an A-number (e.g., A14891)
#    E.g.: /home/janzules/ctDNA_11042024/code/slurmOutput/binning/anumber_file.txt
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

# 4) Load Mamba/conda environment (adjust as needed for your cluster)
module load Mamba  # or module load Anaconda3, etc.
mamba activate pysam_env

# 5) Build the comma-separated string for the Python script
anumber_list_for_python=""
for i in $(seq "$start_index" "$end_index"); do
    anumber_list_for_python+="${anumbers[$i]},"
done
# Remove trailing comma
anumber_list_for_python="${anumber_list_for_python%,}"

# 6) Define input & output directories
INPUT_DIR="/home/janzules/ctDNA_11042024/data/consensusPipeline/intermediate_files/mapped_bam"
OUTPUT_DIR="/home/janzules/ctDNA_11042024/data/human_binned_sequences"

# 7) Run the Python script
# Make sure separate_by_fragments.py is accessible, e.g., in the same dir or in your PATH
python /home/janzules/ctDNA_11042024/code/fragment_analysis/separate_by_fragments.py \
    "$anumber_list_for_python" \
    "$INPUT_DIR" \
    "$OUTPUT_DIR"


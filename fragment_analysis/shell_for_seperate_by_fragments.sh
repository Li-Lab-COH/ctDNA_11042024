#!/bin/bash

##SBATCH --array=0-7                         # Job array index (8 jobs for 16 files, 2 files per job)
#
#
## Calculate the range of anumbers to process for this job array task
#start_index=$(( SLURM_ARRAY_TASK_ID * 2 ))
#end_index=$(( start_index + 1 ))
#
#
## Ensure end_index does not exceed the number of available anumbers
#if (( end_index >= ${#anumbers[@]} )); then
#   end_index=$(( ${#anumbers[@]} - 1 ))
#Fi


#!/bin/bash

# Usage:
#  1) Make this script executable:
#       chmod +x test_run_separate.sh
#  2) Run it:
#       ./test_run_separate.sh

# Parameters for your test:
ANUMBERS="A14891,A14892,A14893"
INPUT_DIR="/Users/janzules/Roselab/ctDNA_11042024/data/testing_binned_analysis/mapped_bam_subset"
OUTPUT_DIR="/Users/janzules/Roselab/ctDNA_11042024/data/testing_binned_analysis/binned_sequences"

# Just a safety check to see if the Python script is in the current directory.
# If it's somewhere else, you can provide the full path:
# PY_SCRIPT="/path/to/separate_by_fragments.py"
PY_SCRIPT="./separate_by_fragments.py"

# Run the Python script with these parameters:
python "$PY_SCRIPT" "$ANUMBERS" "$INPUT_DIR" "$OUTPUT_DIR"

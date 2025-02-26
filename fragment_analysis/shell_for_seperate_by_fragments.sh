#!/bin/zsh

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


bam_file_loc=""
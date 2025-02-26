#!/bin/bash
#SBATCH --job-name=fragment_split
#SBATCH --output=fragment_split_%j.out    # Save output log (with job ID)
#SBATCH --error=fragment_split_%j.err     # Save error log (with job ID)
#SBATCH --time=2:00:00                    # Set max time (hh:mm:ss)
#SBATCH --partition=compute               # Adjust to your cluster's partition
#SBATCH --ntasks=1                        # Single task
#SBATCH --cpus-per-task=4                 # Adjust based on BAM file size
#SBATCH --mem=8G                          # Adjust memory as needed
#SBATCH --mail-user=your_email@example.com # (Optional) Email notifications
#SBATCH --mail-type=END,FAIL              # Email on job completion/failure

# Load Mamba (if necessary)
module load Mamba   # Only if required by your HPC system

# Activate the Mamba environment
mamba activate pysam_env

# Run the Python script
python separate_by_fragments.py input.bam output_directory/

# Deactivate the environment (optional)
mamba deactivate

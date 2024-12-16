#!/bin/bash
#SBATCH --job-name=copy_ctData
#SBATCH --output=/home/janzules/managing_homedirectory/slurmout/copy_ctData.log
#SBATCH --error=/home/janzules/managing_homedirectory/slurmout/copy_ctData.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=10G
#SBATCH --time=24:00:00

cp -r /labs/yunroseli_grp/Jon/ctDNA_11042024/data/20241101_LH00295_0150_B22VNMHLT3/ /home/janzules/ctDNA_11042024/data/

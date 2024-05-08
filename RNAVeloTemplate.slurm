#!/bin/bash
#SBATCH -t 12:00:00
#SBATCH -n 16
#SBATCH --job-name=SCVelo   # Job name
#SBATCH --output=./slurmOut/SCV.%j.out   # Output file
#SBATCH --error=./slurmOut/SCV.%j.err    # Error file

cd /home/gridsan/elu/20.440Project

source /state/partition1/llgrid/pkg/anaconda/anaconda3-2023a/bin/activate
conda init bash
conda activate rnavelo

python3 RNAVelo.py %%SRAID%% 12

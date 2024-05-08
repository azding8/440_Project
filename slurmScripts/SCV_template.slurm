#!/bin/bash
#SBATCH -t 12:00:00
#SBATCH -n 8
#SBATCH --job-name=SCVelo   # Job name
#SBATCH --output=/home/gridsan/elu/20.440Project/slurmOut/SCV.%j.out   # Output file
#SBATCH --error=/home/gridsan/elu/20.440Project/slurmOut/SCV.%j.err    # Error file

source /state/partition1/llgrid/pkg/anaconda/anaconda3-2023a/bin/activate
conda init bash
conda activate rnavelo

cd /home/gridsan/elu/20.440Project

python3 RNAVelo.py %%SRAID%% 8
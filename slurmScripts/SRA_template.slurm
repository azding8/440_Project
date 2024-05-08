#!/bin/bash
#SBATCH -t 6:00:00         
#SBATCH -n 4
#SBATCH --job-name SRA   # Job name
#SBATCH --output /home/gridsan/elu/20.440Project/slurmOut/SRA.%j.out   # Output file
#SBATCH --error /home/gridsan/elu/20.440Project/slurmOut/SRA.%j.err    # Error file
# Go to working directory

cd /home/gridsan/elu/20.440Project/

fasterq-dump /home/gridsan/elu/20.440Project/data/sra/%%SRAID%%.sra --include-technical -S

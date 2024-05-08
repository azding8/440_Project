#!/bin/bash
#SBATCH -t 16:00:00
#SBATCH -n 16
#SBATCH --job-name=AF
#SBATCH --output=/home/gridsan/elu/20.440Project/slurmOut/AF.%j.out
#SBATCH --error=/home/gridsan/elu/20.440Project/slurmOut/AF.%j.err


### Setting up environment  -----------------
export ALEVIN_FRY_HOME="/home/gridsan/elu/20.440Project/AlevinFry"
cd $ALEVIN_FRY_HOME

source /state/partition1/llgrid/pkg/anaconda/anaconda3-2023a/bin/activate
conda init bash
conda activate rnavelo

simpleaf set-paths
ulimit -n 2048


### Setting up paths  -----------------

# Set directory paths
AF_SAMPLE_DIR="$ALEVIN_FRY_HOME/quants"
FASTQ_DIR="/home/gridsan/elu/20.440Project/data/fastq"
IDX_DIR="$ALEVIN_FRY_HOME/%%INDEX%%/human-2020-A_splici3"

# Identify reads
r1_pat="%%SRAID%%_3.fastq"
r2_pat="%%SRAID%%_4.fastq"
reads1="$(find -L ${FASTQ_DIR} -name "*$r1_pat*" -type f | sort | awk -v OFS=, '{$1=$1;print}' | paste -sd, -)"
reads2="$(find -L ${FASTQ_DIR} -name "*$r2_pat*" -type f | sort | awk -v OFS=, '{$1=$1;print}' | paste -sd, -)"


### Run simpleaf  -----------------

simpleaf quant \
--reads1 $reads1 \
--reads2 $reads2 \
--threads 16 \
--index $IDX_DIR/index \
--chemistry 10xv3 --resolution cr-like \
--expected-ori fw --unfiltered-pl \
--t2g-map $IDX_DIR/index/t2g_3col.tsv \
--output $AF_SAMPLE_DIR/%%SRAID%%_quant_v3


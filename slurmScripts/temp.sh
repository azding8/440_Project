for filename in ./slurmScripts/*.slurm; do
    sbatch $filename
done

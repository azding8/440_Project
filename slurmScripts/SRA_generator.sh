SRAIDS=("SRR26383992" "SRR26383994" "SRR26383999")

mkdir "./readySlurm/"
# loop through each SRA ID and create a new slurm file for each
for SRAID in ${SRAIDS[@]}
do
    sed "s/%%SRAID%%/$SRAID/g" SRA_template.slurm > "./readySlurm/SRA_$SRAID.slurm"
done

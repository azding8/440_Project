SRAIDS=("SRR26383970" "SRR26383973" "SRR26383988" "SRR26383989" "SRR26383990" "SRR26383996" "SRR26384000" "SRR26383987" "SRR26383992" "SRR26383994" "SRR26383999")

mkdir "./readySlurm/"
# loop through each SRA ID and create a new slurm file for each
for SRAID in ${SRAIDS[@]}
do
    sed "s/%%SRAID%%/$SRAID/g" SCV_template.slurm > "./readySlurm/SCV_$SRAID.slurm"
done

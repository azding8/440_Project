# Manually input Index directory that contains human-2020-A_splici3, and which SRA files to use
SRAIDS=("SRR26383992" "SRR26383994" "SRR26383999")

# Must follow $AF_SAMPLE_DIR/Index_X_X/human-2020-A_splici
INDEX="Index150_150"

mkdir "./readySlurm/"

# loop through each SRA ID and create a new slurm file for each
for SRAID in ${SRAIDS[@]}
do
    sed "s/%%SRAID%%/$SRAID/g" SAF_template.slurm > "./readySlurm/SAF_$SRAID.slurm"
    sed -i "s/%%INDEX%%/$INDEX/g" "./readySlurm/SAF_$SRAID.slurm"
done

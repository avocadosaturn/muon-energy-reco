#!/bin/bash

#SBATCH --partition=turing
#SBATCH --account=mli:gampix
#
#SBATCH --job-name=gmpix-MPV
#SBATCH --output=logs/output-%j.txt
#SBATCH --error=logs/error-%j.txt
#
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=10g
#SBATCH --gpus=1
#
#SBATCH --time=12:00:00


OUTDIR="/sdf/data/neutrino/summer25/seohyeon/gampix_raw_redo"


EDEPDIR="/sdf/data/neutrino/summer25/seohyeon/edep-sim_54k_raw"
FILES=$(ls $EDEPDIR | sed -n '49,54p')
# echo "$FILES"

GAMPIXROOT=${HOME}/GAMPixPy
SINGULARITY_IMAGE_PATH=/sdf/group/neutrino/images/larcv2_ub22.04-cuda12.1-pytorch2.2.1-larndsim.sif


i=48
for FILE in $FILES; do
    echo "----"
    echo "$FILE"


    OUTPUT_FILE="$OUTDIR/muon1k_0-1gev_gampix_raw_run$i.h5"

    COMMAND="python3 ${GAMPIXROOT}/examples/batch_sim2.py $EDEPDIR/${FILE} -i root -o ${OUTPUT_FILE} -r ${GAMPIXROOT}/gampixpy/readout_config/GAMPixD_notruth.yaml"

    singularity exec --nv -B /sdf,/lscratch ${SINGULARITY_IMAGE_PATH} ${COMMAND}

    ((i++))
done









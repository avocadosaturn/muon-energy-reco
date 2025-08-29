#!/bin/bash

#SBATCH --partition=turing
#SBATCH --account=mli:gampix
#
#SBATCH --job-name=gmpix-analyzer
#SBATCH --output=logs/output-%j.txt
#SBATCH --error=logs/error-%j.txt
#
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=10g
#SBATCH --gpus=1
#
#SBATCH --time=10:00:00

IMAGE_PATH=/sdf/group/neutrino/images/develop.sif
COMMAND="python directory_gampix_analyzer.py"

singularity exec --nv -B /sdf,/lscratch ${IMAGE_PATH} ${COMMAND}




#!/bin/bash

#SBATCH --partition=turing
#SBATCH --account=mli:gampix
#
#SBATCH --job-name=muon-analysis
#SBATCH --output=logs/output-%j.txt
#SBATCH --error=logs/output-%j.txt
#
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=10g
#SBATCH --gpus=1
#
#SBATCH --time=10:00:00
#
#SBATCH --array=0

SINGULARITY_IMAGE_PATH=/sdf/group/neutrino/images/larcv2_ub20.04-cuda11.3-cudnn8-pytorch1.10.0-larndsim-2022-11-03.sif
COMMAND="python track_length_analysis.py"


singularity exec -B /sdf,/lscratch ${SINGULARITY_IMAGE_PATH} ${COMMAND}
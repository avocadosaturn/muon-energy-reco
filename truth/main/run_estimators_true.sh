#!/bin/bash

#SBATCH --partition=turing
#SBATCH --account=mli:gampix
#
#SBATCH --job-name=edep-estimators
#SBATCH --output=logs/output-%j.txt
#SBATCH --error=logs/error-%j.txt
#
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=10g
#SBATCH --gpus=1
#
#SBATCH --time=10:00:00


singularity exec --nv -B /sdf,/lscratch /sdf/group/neutrino/images/develop.sif python3 estimators_true.py

#!/bin/bash

#SBATCH --partition=turing
#SBATCH --account=mli:gampix
#
#SBATCH --job-name=gmpix-MPV
#SBATCH --output=logs/output-%j.txt
#SBATCH --error=logs/output-%j.txt
#
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=10g
#SBATCH --gpus=1
#
#SBATCH --time=10:00:00
#
#SBATCH --array=0-9

SINGULARITY_IMAGE_PATH=/sdf/group/neutrino/images/larcv2_ub20.04-cuda11.3-cudnn8-pytorch1.10.0-larndsim-2022-11-03.sif

OUTDIR=$1

# I've been using a UUID from the OS to make unique filenames
# these are very long and can be quite ugly
# you can also just use the slurm batch index
ID=$(cat /proc/sys/kernel/random/uuid) 

G4MACRO=$2

EDEPDIR=${OUTDIR}/edep-sim2
mkdir -p $EDEPDIR
EDEPOUTPUT=${EDEPDIR}/muon_0-1gev_1000ev_run0.root
NEVENTS=1000

COMMAND="edep-sim -g SeaOfArgon.gdml -e ${NEVENTS} -o ${EDEPOUTPUT} ${G4MACRO}"

singularity exec -B /sdf,/lscratch ${SINGULARITY_IMAGE_PATH} ${COMMAND}
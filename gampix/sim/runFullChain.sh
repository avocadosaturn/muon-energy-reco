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
#SBATCH --time=10:00:00
#SBATCH --array=0-19

SINGULARITY_IMAGE_PATH=/sdf/group/neutrino/images/larcv2_ub20.04-cuda11.3-cudnn8-pytorch1.10.0-larndsim-2022-11-03.sif

OUTDIR=$1

# I've been using a UUID from the OS to make unique filenames
ID=$(cat /proc/sys/kernel/random/uuid) 

G4MACRO=$2

EDEPDIR=${OUTDIR}/edep-sim
mkdir -p $EDEPDIR
EDEPOUTPUT=${EDEPDIR}/muon1k_edep_test.root
NEVENTS=1000

COMMAND="edep-sim -g SeaOfArgon.gdml -e ${NEVENTS} -o ${EDEPOUTPUT} ${G4MACRO}"

singularity exec -B /sdf,/lscratch ${SINGULARITY_IMAGE_PATH} ${COMMAND}

# #-----------------------

SINGULARITY_IMAGE_PATH=/sdf/group/neutrino/images/larcv2_ub22.04-cuda12.1-pytorch2.2.1-larndsim.sif

GAMPIXROOT=${HOME}/GAMPixPy
GAMPIXDDIR=${OUTDIR}/gampix
mkdir -p $GAMPIXDDIR
GAMPIXDOUTPUT=${GAMPIXDDIR}/muon1k_gampix_test.h5

COMMAND="python3 ${GAMPIXROOT}/examples/batch_sim.py -i root ${EDEPOUTPUT} -o ${GAMPIXDOUTPUT} -r ${GAMPIXROOT}/gampixpy/readout_config/GAMPixD.yaml"

singularity exec --nv -B /sdf,/lscratch ${SINGULARITY_IMAGE_PATH} ${COMMAND}

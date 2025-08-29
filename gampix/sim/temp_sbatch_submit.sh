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


singularity exec --nv -B /sdf,/lscratch /sdf/group/neutrino/images/larcv2_ub22.04-cuda12.1-pytorch2.2.1-larndsim.sif python3 ~/GAMPixPy/examples/batch_sim2.py /sdf/data/neutrino/summer25/seohyeon/edep-sim_54k_raw/edep_single_particle_7568eea7-38aa-46e0-becd-b5f9b79521cf.root -i root -o /sdf/data/neutrino/summer25/seohyeon/gampix_raw_redo/muon1k_0-1gev_gampix_raw_run29.h5 -r ~/GAMPixPy/gampixpy/readout_config/GAMPixD_notruth.yaml

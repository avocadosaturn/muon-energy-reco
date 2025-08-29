import h5py
import sys
import numpy as np
import matplotlib.pyplot as plt
import os
from tools import bin20

m_mu = 105.7 #MeV
energies=[]
deposited_Es=[]
ranges=[]
path_lengths=[]
evs=[]

directory = '/sdf/data/neutrino/summer25/seohyeon/edep-sim_h5_corrected/'

for filename in os.listdir(directory):
    f = h5py.File(directory + filename)

    print(filename)
    
    trajs = f['trajectories']
    segs = f['segments']

    
    primary_mask_trajs = trajs['traj_id']== 0
    muon_trajs = trajs[primary_mask_trajs] # primary muon trajectories only
    
    for traj in muon_trajs:
        # getting event number
        ev_n = traj['event_id']
        evs.append(ev_n)
        # print('----------')
        # print(f'event: {ev_n}')
    
        # getting energies
        KE = traj['E_start'] - m_mu
        energies.append(KE)
    
        # getting range ---------------
        start_vector = traj['xyz_start']
        end_vector = traj['xyz_end']
        distance_vector = np.subtract(end_vector, start_vector)
        length = np.linalg.norm(distance_vector) *10 #mm
        ranges.append(length)
    
        # getting path length -------------
        # getting segments in muon track only
        ev_mask = segs['event_id'] == ev_n #ev_mask is a boolean array, true if that segment is in that trajectory
        ev_segs = segs[ev_mask] #only segments in that event
        primary_mask = ev_segs['traj_id'] == 0 
        muon_segs = ev_segs[primary_mask] #only segments in primary muon track
    
        path_length = 0
        for seg in muon_segs:
            path_length += seg['dx'] * 10
        path_lengths.append(path_length)

        # getting deposited energy
        deposited_E = 0
        for seg in muon_segs:
            deposited_E += seg['dE']
        deposited_Es.append(deposited_E)



#pdg data
#https://pdg.lbl.gov/2011/AtomicNuclearProperties/MUON_ELOSS_TABLES/muonloss_289.pdf
Ar_density = 1.396 #g/cm^3
pdg_data_energies=np.array([10, 14, 20, 30, 40, 80, 100, 140, 200, 300, 400, 800, 1000]) #MeV
pdg_data_CSDA_range=np.array([0.9833, 1.786,  3.321, 6.598, 10.58, 30.84, 42.50, 67.32, 106.3, 172.5, 238.5, 493.4, 616.3]) #g/cm^2
pdg_data_range= (pdg_data_CSDA_range / Ar_density) *10 #mm



total_evs = np.max(evs) + 1
# range plots ------------------------------------------------------
# range vs. KE
plt.figure()
plt.plot(energies, ranges, 'o')
plt.xlabel('muon kinetic energy (MeV)')
plt.ylabel('range (cm)')
plt.title(f'range vs. muon KE, {total_evs} events')
plt.savefig(f'Plots/range_{total_evs}evs.png')

# # profile
range_bin_midpts, range_avgs, range_rms, range_relative_rms = bin20(energies, ranges)
plt.figure()
plt.errorbar(range_bin_midpts, range_avgs, range_rms, fmt='o')
plt.xlabel('muon kinetic energy (MeV)')
plt.ylabel('range (cm)')
plt.title(f'range vs. muon KE profile, {total_evs} events')
plt.savefig(f'Plots/range_profile_{total_evs}evs.png')

#relative rms
plt.figure()
plt.plot(range_bin_midpts, range_relative_rms, 'o')
plt.xlabel('muon KE (MeV)')
plt.ylabel('relative rms')
plt.title(f'relative variance in range vs. muon KE, {total_evs} events')
plt.savefig(f'Plots/range_rms_{total_evs}evs.png')


# path length plots---------------------------------------
# path length vs. KE
plt.figure()
plt.plot(energies, path_lengths, 'o')
plt.xlabel('muon kinetic energy (MeV)')
plt.ylabel('path length (cm)')
plt.title(f'path length vs. muon KE, {total_evs} events')
plt.savefig(f'Plots/path_length_{total_evs}evs.png')

# #profile
pl_bin_midpts, pl_avgs, pl_rms, pl_relative_rms = bin20(energies, path_lengths)
plt.figure()
plt.errorbar(pl_bin_midpts, pl_avgs, pl_rms, fmt='o')
plt.xlabel('muon kinetic energy (MeV)')
plt.ylabel('path length (cm)')
plt.title(f'path length vs. muon KE profile, {total_evs} events')
plt.savefig(f'Plots/path_length_profile_{total_evs}evs.png')

# relative rms
# relative rms is calculated relative to average path length of that bin. but that is ~proportional to muon KE?
plt.figure()
plt.plot(pl_bin_midpts, pl_relative_rms, 'o')
plt.xlabel('muon KE (MeV)')
plt.ylabel('relative rms')
plt.title(f'relative variance in path length vs. muon KE, {total_evs} events')
plt.savefig(f'Plots/path_length_rms_{total_evs}evs.png')


# energy deposition plots----------------------------------------------
#energy deposited vs. initial KE
plt.figure()
plt.plot(energies, deposited_Es, 'o')
plt.xlabel('muon kinetic energy (MeV)')
plt.ylabel('deposited E (MeV)')
plt.title(f'deposited energy vs. muon KE, {total_evs} events')
plt.savefig(f'Plots/energy_deposited_KE_{total_evs}evs.png')

#profile
KE_bin_midpts, KE_avgs, KE_rms, KE_relative_rms = bin20(initial_KEs, deposited_Es)
plt.figure()
plt.errorbar(KE_bin_midpts, KE_avgs, KE_rms, fmt='o')
plt.xlabel('muon kinetic energy (MeV)')
plt.ylabel('deposited E (MeV)')
plt.title(f'deposited energy vs. muon KE profile, {total_evs} events')
plt.savefig(f'Plots/energy_deposited_KE_profile_{total_evs}evs.png')

#relative rms
plt.figure()
plt.plot(KE_bin_midpts, KE_relative_rms, 'o')
plt.xlabel('muon KE (MeV)')
plt.ylabel('relative rms')
plt.title(f'relative variance in deposted energy vs. muon KE, {total_evs} events')
plt.savefig(f'Plots/energy_deposited_KE_relativerms_{total_evs}evs.png')


# comparison plots----------------------------------------
# range & path length vs. KE
plt.figure()
plt.errorbar(pl_bin_midpts, pl_avgs, pl_rms, fmt='o', color='blue', alpha=0.9, label='path length')
plt.errorbar(range_bin_midpts, range_avgs, range_rms, fmt='o', color='orange', alpha=0.65, label='range')
plt.xlabel('muon KE (MeV)')
plt.ylabel('length (cm)')
plt.title(f'length estimates vs. muon KE, {total_evs} events')
plt.legend()
plt.savefig(f'Plots/length_estimates_{total_evs}evs.png')

# # range vs. path length
# plt.figure()
# plt.plot(path_lengths, ranges, 'o')
# plt.xlabel('path length (cm)')
# plt.ylabel('range (cm)')
# plt.title(f'range vs. path length, {total_evs} events')
# plt.savefig(f'Plots/range_vs_path_length_{total_evs}evs.png')

# # range with pdg data 
# plt.figure()
# plt.errorbar(range_bin_midpts, range_avgs, range_rms, fmt='o', label='edepsim')
# plt.plot(pdg_data_energies, pdg_data_range, 'o', color='orange', label='pdg')
# plt.xlabel('muon kinetic energy (MeV)')
# plt.ylabel('range (cm)')
# plt.title(f'range vs. muon KE, {total_evs} events')
# plt.legend()
# plt.savefig(f'Plots/range_with_pdg_{total_evs}evs.png')


# rms comparison
KE_bin_midpts, KE_avgs, KE_rms, KE_relative_rms = bin20(energies, deposited_Es)
plt.figure()
plt.plot(range_bin_midpts, range_relative_rms, 'o', label='range', markerfacecolor='none')
plt.plot(pl_bin_midpts, pl_relative_rms, 'o', label='path length', markerfacecolor='none')
plt.plot(KE_bin_midpts, KE_relative_rms, 'o', label='energy deposited', markerfacecolor='none') 
plt.xlabel('muon KE (MeV)')
plt.ylabel('relative rms')
plt.legend()
plt.savefig(f'Plots/rms_comparison_{total_evs}evs.png')



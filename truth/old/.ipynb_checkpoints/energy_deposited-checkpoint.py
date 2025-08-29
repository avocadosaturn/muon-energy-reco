import h5py
import sys
import numpy as np
import matplotlib.pyplot as plt
from tools import bin20 

if len(sys.argv) != 2:
    print("usage: python [program].py [file.h5]")
    quit()

m_mu = 105.7 #MeV
deposited_Es=[]
initial_KEs=[]
ranges=[]
evs=[]

f = h5py.File('/sdf/home/s/seohyeon/summer25/data/' + sys.argv[1])
trajs = f['trajectories']
segs = f['segments']

primary_mask_trajs = trajs['traj_id']== 0
muon_trajs = trajs[primary_mask_trajs] # primary muon trajectories only
        
for traj in muon_trajs:
    # getting event number 
    ev_n = traj['event_id'] #these are in order: 0, 1, 2, ...
    evs.append(ev_n)
    print('----------')
    print(f'event: {ev_n}')

    # getting initial KE
    KE = traj['E_start'] - m_mu
    initial_KEs.append(KE)

    # getting range
    start_vector = traj['xyz_start']
    end_vector = traj['xyz_end']
    distance_vector = np.subtract(end_vector, start_vector)
    length = np.linalg.norm(distance_vector) 
    ranges.append(length)

    # getting segments in that muon track only
    ev_mask = segs['event_id'] == ev_n #ev_mask is a boolean array, true if that segment is in that trajectory
    ev_segs = segs[ev_mask] #only segments in that event
    primary_mask = ev_segs['traj_id'] == 0 
    muon_segs = ev_segs[primary_mask] #only segments in primary muon track

    # getting deposited energy
    deposited_E = 0
    for seg in muon_segs:
        deposited_E += seg['dE']
    deposited_Es.append(deposited_E)
    


total_evs = np.max(evs) + 1
#----------------------------------------------
#energy deposited vs. initial KE
plt.figure()
plt.plot(initial_KEs, deposited_Es, 'o')
plt.xlabel('muon kinetic energy (MeV)')
plt.ylabel('deposited E (MeV)')
plt.title(f'deposited energy vs. muon initial KE, {total_evs} events')
plt.savefig(f'Plots/energy_deposited_KE_{total_evs}evs.png')

#profile
KE_bin_midpts, KE_avgs, KE_rms, KE_relative_rms = bin20(initial_KEs, deposited_Es)
plt.figure()
plt.errorbar(KE_bin_midpts, KE_avgs, KE_rms, fmt='o')
plt.xlabel('muon kinetic energy (MeV)')
plt.ylabel('deposited E (MeV)')
plt.title(f'deposited energy vs. muon initial KE, {total_evs} events')
plt.savefig(f'Plots/energy_deposited_KE_profile_{total_evs}evs.png')

#relative rms
# relative rms is calculated relative to average path length of that bin. but that is ~proportional to muon KE 
plt.figure()
plt.plot(KE_bin_midpts, KE_relative_rms, 'o')
plt.xlabel('muon KE (MeV)')
plt.ylabel('relative rms')
plt.title(f'relative variance in deposted energy vs. muon KE, {total_evs} events')
plt.savefig(f'Plots/energy_deposited_KE_relativerms_{total_evs}evs.png')

#----------------------------------------------
# energy deposited vs. range
plt.figure()
plt.plot(ranges, deposited_Es, 'o')
plt.xlabel('range (cm)')
plt.ylabel('deposited E (MeV)')
plt.title(f'deposited energy vs. range, {total_evs} events')
plt.savefig(f'Plots/energy_deposited_range_{total_evs}evs.png')

#profile
range_bin_midpts, range_avgs, range_rms, range_relative_rms = bin20(ranges, deposited_Es)
plt.figure()
plt.errorbar(range_bin_midpts, range_avgs, range_rms, fmt='o')
plt.xlabel('range (cm)')
plt.ylabel('deposited E (MeV)')
plt.title(f'deposited energy vs. range, {total_evs} events')
plt.savefig(f'Plots/energy_deposited_range_profile_{total_evs}evs.png')

#relative rms
# relative rms is calculated relative to average path length of that bin. but that is ~proportional to muon KE 
plt.figure()
plt.plot(range_bin_midpts, range_relative_rms, 'o')
plt.xlabel('range (cm)')
plt.ylabel('relative rms')
plt.title(f'relative variance in deposited energy vs. range, {total_evs} events')
plt.savefig(f'Plots/energy_deposited_range_relativerms_{total_evs}evs.png')



import h5py
import sys
import numpy as np
import matplotlib.pyplot as plt
import os



m_mu = 105.7 #MeV


directory = '/sdf/data/neutrino/summer25/seohyeon/edep-sim_h5_54k_raw/'
# files = [f for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f))]

# print(files)

file_ns = range(54)

for file_n in file_ns:
# for filename in files:
    energy_deposited=[]
    
    filename = f'muon_0-1gev_run{file_n}.h5'
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
        print('----------')
        print(f'event: {ev_n}')

        # getting energies
        KE = traj['E_start'] - m_mu
        energies.append(KE)
    
        # getting range ---------------
        start = traj['xyz_start']
        stop = traj['xyz_end']
        distance_vector = np.subtract(stop, start)
        length = np.linalg.norm(distance_vector) * 10 #mm
        ranges.append(length)
    
        # getting path length -------------
        # getting segments in muon track only
        ev_mask = segs['event_id'] == ev_n #ev_mask is a boolean array, true if that segment is in that trajectory
        ev_segs = segs[ev_mask] #only segments in that event
        primary_mask = ev_segs['traj_id'] == 0 
        muon_segs = ev_segs[primary_mask] #only segments in primary muon track

        seg_starts = np.array([[seg['x_start'], seg['y_start'], seg['z_start']] for seg in muon_segs])
        seg_ends = np.array([[seg['x_end'], seg['y_end'], seg['z_end']] for seg in muon_segs])

        path_length = 0

        if np.any(start != seg_starts[0]):
            path_length += np.linalg.norm(start - seg_starts[0])
            print('mismatch')
            print(start)
            print(seg_starts[0])
            print(len(seg_starts))
        elif np.any(stop != seg_ends[-1]):
            path_length += np.linalg.norm(stop - seg_ends[-1])
            print('mismatch')
            print(start)
            print(seg_starts[0])
            print(len(seg_starts))

        for i, seg in enumerate(muon_segs):
            if i+1 == len(muon_segs):
                path_length += np.linalg.norm(seg_ends[i] - seg_starts[i]) 
            else:
                path_length += np.linalg.norm(seg_ends[i] - seg_starts[i]) 
                path_length += np.linalg.norm(seg_starts[i+1] - seg_ends[i]) 
        path_lengths.append(path_length * 10)


        # getting deposited energy
        ev_mask = segs['event_id'] == ev_n #ev_mask is a boolean array, true if that segment is in that trajectory
        ev_segs = segs[ev_mask] #only segments in that event
        primary_mask = ev_segs['traj_id'] == 0 
        muon_segs = ev_segs[primary_mask] #only segments in primary muon track

        deposited_E = 0
        for seg in muon_segs:
            deposited_E += seg['dE']
        deposited_Es.append(deposited_E)

    np.savez(f'/sdf/data/neutrino/summer25/seohyeon/temp/muon1k_0-1gev_estimators_edep_run{file_n}.npz', energies=energies, naive=ranges, detailed=path_lengths, edeps=deposited_Es, evs=evs)
import h5py
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from matplotlib.widgets import Slider
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import binned_statistic_2d
from scipy.interpolate import make_splprep


directory = '/sdf/data/neutrino/summer25/seohyeon/gampix_raw_redo/'
# files = [f for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f))]
# print(files)

file_ns = [19]
# counter = 0

for file_n in file_ns:
# for filename in files:
    # initialization
    energies=[]
    naive_lengths=[]
    detailed_lengths=[]
    charge_collected=[]
    evs_passed = []
    
    
    # parameters
    noise = 50 #e-
    threshold_sigma = 4
    e_lifetime = 10 #ms
    work_fxn = 2.36e-5 #MeV


    # read file
    # f = h5py.File(directory + filename)
    filename = f'muon1k_0-1gev_gampix_raw_run{file_n}.h5'
    f = h5py.File(directory + filename)
    print(filename)

    
    pixels = f['pixel_hits']
    coarse = f['coarse_hits']
    meta = f['meta']
    evs = f['meta']['event id']
    
    
    for ev_n in evs:
        ev_mask_pixels = pixels['event id'] == ev_n
        ev_pixels = pixels[ev_mask_pixels]
        # print('-------')
        # print(ev_n)
        
        points = np.array([[hit['pixel x'], hit['pixel y'], hit['hit z']] for hit in ev_pixels])
        charge = ev_pixels['hit charge']
    
        # detector effects ---------------------------------------------
        # noise_mask = charge - (threshold_sigma * noise) > 0
    
        
        # charge = charge[noise_mask] # pixel noise cut
        # points = points[noise_mask]
        
        # charge = (charge * np.exp(ev_pixels['hit t'])) / e_lifetime # attenuation
        
        # #TODO: recombination
        
        # charge = charge * work_fxn # work fxn
    x
        # event cuts ---------------------------------------------------
        # no hits in event
        if len(points) == 0:
            # print('no hits in this event')
            continue
        
        # not enough hits for a pca analysis
        min_hits = 40
        if len(points) < min_hits:
            # print('not enough hits in this event')
            continue
    
        # track runs too close too anode
        if min(points[:, 2]) < 5:
            # print('skipped, too close to anode')
            continue
    
        
        evs_passed.append(ev_n)
        energies.append(meta['primary energy'][ev_n])
        charge_collected.append(np.sum(charge))
        
        # naive length ----------------------------------------------------------
        # pca
        pca = PCA(n_components = 3)
        pca_on_hits = pca.fit(points)
        principal=pca.components_[0] #this is a unit vector
    
        projections ={}
        for pt in points:
            mag_proj = np.dot(principal, pt) #signed projection scalar
            projections[mag_proj] = pt
        
        
        start = projections[min(projections)]
        stop = projections[max(projections)]
        start = np.array(start)
        stop = np.array(stop)
    
        naive_length = max(projections) - min(projections)
        naive_vector = (stop - start) * naive_length / np.linalg.norm(stop - start)
    
        naive_lengths.append(naive_length)
        # print(f'naive: {naive_length}')
    
        # detailed length -------------------------------------------------------
        # 2d binning  
        bins1d = 5
        statistic, x_edges, y_edges, binnumber = binned_statistic_2d(
            points[:, 0], points[:, 1], points[:, 2], statistic='mean', bins=(bins1d, bins1d))
    
        
        # calculating the weighted avg of each bin
        weighted_avgs = {}
        for bin_n in range(max(binnumber)+1):
            indices = np.where(binnumber==bin_n)
            indices = indices[0]
            
            bin_points = points[indices]
            bin_charges = charge[indices]
        
            if len(bin_points) == 0:
                continue
                
            weighted_avg = 0
            num = 0
            denom = 0
        
            
            for i, pt in enumerate(bin_points):
                num += pt[2] * bin_charges[i]
            for val in bin_charges:
                denom += val
                
            weighted_avg = num / denom
        
            if np.isnan(weighted_avg) == True:
                continue
            else:
                weighted_avgs[bin_n] = weighted_avg
    
    
        # calculating the middle (x, y) of each bin
        spl_pts=[]
        for bin_n in weighted_avgs.keys():
            indices = np.where(binnumber==bin_n)
            indices = indices[0]
            
            bin_points = points[indices]
        
            x_avg = 0
            y_avg = 0
            n = len(bin_points)
        
            for pt in bin_points:
                x_avg += pt[0]
                y_avg += pt[1]
                
            x_avg = x_avg / n
            y_avg = y_avg / n
        
        
            spl_pt = [x_avg, y_avg, weighted_avgs[bin_n]]
            spl_pts.append(spl_pt)
    
        
        # sort binned points for spline
        spl_pts = np.array(spl_pts)
        
        dist_from_start = []
        for pt in spl_pts:
            dist_from_start.append(np.linalg.norm(pt - start))
        
        sorted_indices = np.argsort(dist_from_start)
        spl_pts = spl_pts[sorted_indices]
        
        
        # add start and stop points to spline points
        start_index = np.where(np.all(points == start, axis = 1))[0]
        stop_index = np.where(np.all(points == stop, axis = 1))[0]
        spl_pts = np.insert(spl_pts, 0, start, axis=0)
        spl_pts = np.vstack([spl_pts, stop])
    
        
        # weight start and stop points to clamp spline
        clamp_weights = np.ones(len(spl_pts) - 2)
        clamp_weights = np.insert(clamp_weights, 0, 1000)
        clamp_weights = np.append(clamp_weights, 1000)
        
        
        # spline 
        spline, u = make_splprep([spl_pts[:, 0], spl_pts[:, 1], spl_pts[:, 2]], s=25, w=clamp_weights)
        u_fine = np.linspace(0, 1, 1000)
        x_fine, y_fine, z_fine = spline(u_fine)
        
        
        # arc length of spline
        arc = 0
        u_ultrafine = np.linspace(0, 1, 5000)
        for i, u in enumerate(u_ultrafine):
            if i == 0 :
                continue
            else:
                u_prev = u_ultrafine[i-1]        
                arc += np.linalg.norm(spline(u) - spline(u_prev))
        detailed_lengths.append(arc)
    
        
        # sanity check 
        if naive_length > arc:
            print('naive > detailed!')
            print(f'naive length: {naive_length}')
            print(f'detailed length: {arc}')

    
    np.savez(f'/sdf/data/neutrino/summer25/seohyeon/gampix_estimators_redo/muon1k_0-1gev_estimators_gampix_run{file_n}.npz', naive = naive_lengths, detailed=detailed_lengths, charge=charge_collected, evs=evs_passed, energies=energies)


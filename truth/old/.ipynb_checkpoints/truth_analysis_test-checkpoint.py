#!/usr/bin/env python
# coding: utf-8

# In[2]:


import h5py
import sys
import numpy as np
import matplotlib.pyplot as plt
import os
from tools import bin20


# In[3]:


m_mu = 105.7 #MeV
energies=[]
deposited_Es=[]
ranges=[]
path_lengths=[]
evs=[]

f = h5py.File('/sdf/home/s/seohyeon/summer25/test_data/muon/muon1000ev_0-1gev_edep.h5')

trajs = f['trajectories']
segs = f['segments']


primary_mask_trajs = trajs['traj_id'] == 0
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



# In[112]:


import importlib
import tools
from scipy.optimize import curve_fit
importlib.reload(tools)

from tools import bin20  


# In[93]:


total_evs=1000

# # method 1: fit with polynomial
# # range plot fitting
# range_coeffs = np.polyfit(ranges, energies, deg=2)
# range_poly = np.poly1d(range_coeffs)
# range_yfit = range_poly(ranges)
# range_residuals = range_yfit-energies

# # edep plot fitting
# edep_coeffs = np.polyfit(deposited_Es, energies, deg=2)
# edep_poly = np.poly1d(edep_coeffs)
# edep_yfit = edep_poly(deposited_Es)
# edep_residuals = edep_yfit - energies

range_bmin = -np.min(ranges) + 1e-6  # so x + b > 0 always

range_lbs = [0, 0, range_bmin, -np.inf, -np.inf, -np.inf]  # [x0, a, b, c, m, d]
range_ubs = [np.inf, np.inf, np.inf, np.inf, np.inf, np.inf]

# method 2: piecewise fit
def piecewise(x, x0, a, b, c, m, d):
    return np.piecewise(x, [x < x0, x >= x0],
                        [lambda x: a * np.sqrt(x + b) + c,
                         lambda x: m*x + d])

# Initial guess: [x0, a, b, c, m, d]
range_p0 = [110, 12, 130, -150, 0.22, 30]

range_params, _ = curve_fit(piecewise, ranges, energies, range_p0, bounds=(range_lbs, range_ubs))
print(f"Fitted breakpoint at x0 = {range_params[0]:.3f}")

range_xfit_display = np.linspace(min(ranges), max(ranges), 500)
range_yfit_display = piecewise(range_xfit_display, *range_params)


# # this plot was for figuring out intial guess for the sqrt part
# plt.figure()
# plt.xlim(0, 1000)
# plt.plot(ranges, energies, 'o')
# x1 = np.linspace(0.01, 1000, 100)
# plt.plot(x1, 12*np.sqrt(x1+130) - 150)


# range plots ------------------------------------------------------
# KE vs. range

plt.figure()
plt.plot(ranges, energies, 'o')

# this was for figuring out initial guess for linear part
# x2 = np.linspace(0.01, 4500, 500)
# plt.plot(x2, 0.22 * x2 + 30)

# plt.plot(ranges, range_yfit) # quadratic fit

plt.plot(range_xfit_display, range_yfit_display, label="Piecewise Fit", color="red", linewidth=2)
plt.xlabel(r'$r_n$ (mm)', fontsize=12)
plt.ylabel(r'$T_\mu$ (MeV)', fontsize=12)
plt.title(f'true muon KE vs. range, {total_evs} events', fontsize=16)
plt.savefig(f'Plots/range_{total_evs}evs.png')
plt.close()

range_yfit=(piecewise(ranges, *params))
range_res = (range_yfit - energies)
plt.figure()
plt.plot(energies, range_res, 'o')
plt.close()


# In[115]:


# edep_bmin = -np.min(deposited_Es) + 1e-6  # so x + b > 0 always

# edep_lbs = [0, 0, edep_bmin, -np.inf, -np.inf, -np.inf]  # [x0, a, b, c, m, d]
# edep_ubs = [np.inf, np.inf, np.inf, np.inf, np.inf, np.inf]

# # method 2: piecewise fit

# # Initial guess: [x0, a, b, c, m, d]
# edep_0 = [110, 12, 130, -150, 0.22, 30]

# edep_params, _ = curve_fit(piecewise, deposited_Es, energies, edep_p0, bounds=(range_lbs, range_ubs))
# print(f"Fitted breakpoint at x0 = {edep_params[0]:.3f}")

# edep_xfit_display = np.linspace(min(ranges), max(ranges), 500)
# edep_yfit_display = piecewise(edep_xfit_display, *edep_params)
deposited_Es = np.array(deposited_Es)

edep_coeffs = np.polyfit(deposited_Es, energies, 1)  # returns [slope, intercept]
edep_slope, edep_intercept = edep_coeffs

edep_xfit_display = np.linspace(0, 900, 100)
edep_yfit_display = (edep_slope * edep_xfit_display) + edep_intercept

edep_yfit = (edep_slope * deposited_Es) + edep_intercept
edep_res = edep_yfit - energies

plt.figure()
plt.plot(deposited_Es, energies, 'o')
plt.plot(edep_xfit_display, edep_yfit_display)
plt.close()

res_bin_midpts, res_avgs, res_rms, res_relative_rms = bin20(edep_residuals, range_residuals)
plt.figure()
plt.plot(edep_res, range_res, 'o')

plt.figure()
plt.errorbar(res_bin_midpts, res_avgs, res_rms, fmt='o')

# plt.figure()
# plt.hist2d(edep_res, range_res, bins=60)


# In[ ]:





# In[ ]:





#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#%%
"""
Created on Thu Jun 1 11:32:34 2023
@author: HYUNHEE KWAK
Purpose: get autocorrelation curve from binned photon arrival data
        Paul Müller (2012) Python multiple-tau algorithm (Version 0.3.3). Available at https://pypi.python.org/pypi/multipletau/ (Accessed 01 JUN 2023)

 ---variables---
 below 5 variables are user input variables 
     filename : Strings, filepath or filename in current directory
     bin_width_sec: float, binning width in second, desired input value 
     resolution: int, number of points on each level of lag-time, when to scale up
     time_gate_op: boolean, if time-gate option applied
     time_gate_valid: float, desired range of valid time when time_gate_op is True
 
 
 nanotimes = photon_data.nanotimes.read()
 tcspc_unit = photon_data.nanotimes_specs.tcspc_unit.read()
 tcspc_num_bins = photon_data.nanotimes_specs.tcspc_num_bins.read()
 tcspc_range = photon_data.nanotimes_specs.tcspc_range.read()
 
 ts:            array-like, photon arrival time stamps in integer number 
 ts_unit:       float, unit time of timestamp in second, fixed information from fluoroscence spectroscopy
 bin_width_n:   int, binning width in timestamp unit, calculated from above two
 binned_data:   array-like, binned count of photon arrival
             
 autocorr_array: array-like, autocorrelation value along multiple tau (lag time in second) in logarithmic scale
        m = even integer, resolution (number of points on one level, when to scale up)
        deltat = bin width in second for multiple tau autocorrelate
"""

#%% User input values

filename = "test.hdf5"
bin_width_sec = 200* 10**-9  #second
resolution = 8
time_gate_op = True
time_gate_valid = None


#%% Data read
#%matplotlib inline #This line is for Jupyter Notebook environment
import os
import numpy as np
import tables
import matplotlib.pyplot as plt
import multipletau as mt
from binning import photon_bin


if os.path.exists(filename):
    print('File found, it will be proceed.')
else:
    print('ATTENTION: File not found, please check the directory name.\n'
          '           (current value "%s")' % filename)


h5file = tables.open_file(filename)
photon_data = h5file.root.photon_data
photon_data.measurement_specs.measurement_type.read().decode()

nanotimes = photon_data.nanotimes.read()
tcspc_unit = photon_data.nanotimes_specs.tcspc_unit.read()
tcspc_num_bins = photon_data.nanotimes_specs.tcspc_num_bins.read()
tcspc_range = photon_data.nanotimes_specs.tcspc_range.read()
ts = photon_data.timestamps.read()
ts_unit = photon_data.timestamps_specs.timestamps_unit.read() # ~50ns

#%% Figure 1: Tcspc histogram 

# figure 1 nanotimes histogram
fig1, ax = plt.subplots()
ax.hist(nanotimes, bins = tcspc_num_bins)

# figure 1 x-axis (nanotimes in unit count) 
ax.tick_params(top = True, labeltop = True, bottom = False, labelbottom = False, labelsize=15)
font_x = {'size':15}
ax.set_title('Time', fontdict = font_x) # used as Xaxis(top) label

# figure 1 second x-axis (nanotimes in nano-second)
def nt2sec(x):
    return x * tcspc_unit * 1E9 # nano second
def sec2nt(x):
    return x / tcspc_unit / 1E9

secax = ax.secondary_xaxis('bottom', functions=(nt2sec, sec2nt))
secax.set_xlabel('Time (ns)', fontdict = font_x)
secax.xaxis.set_tick_params(labelsize=15)

# figure 1 y-axis (photone count)
ax.set_ylabel('Counts in log scale', fontdict = font_x)
ax.set_yscale("log")
ax.yaxis.set_tick_params(labelsize=15)


#%% Figure 2: Autocorrelation graph

binwidth_n = int(bin_width_sec / ts_unit) # binwidth in numbers : how many cycles per bin
binned_data = photon_bin(ts, binwidth_n).astype(np.float64)
autocorr_array = mt.autocorrelate(binned_data, m = resolution, deltat = bin_width_sec, normalize=True, copy=True, dtype=None, compress='average', ret_sum=False)


# figure 2 autocorrelation graph
plt.rcParams.update({'figure.figsize':(7,5), 'figure.dpi':100})

plt.figure(2)
plt.plot(autocorr_array[:,0], autocorr_array[:, 1])

# figure 2 x-axis (lag-time in second) to Logarithmic scale 
plt.xticks(fontsize=15)
plt.xlabel('lag time (s)', fontdict = font_x)
plt.xscale("log")

# figure 2 y-axis (extent of correlation)
font_y = {'size':15}
plt.yticks(fontsize=15)
plt.ylim(0, 1)
plt.ylabel('g(τ)', fontdict = font_y)



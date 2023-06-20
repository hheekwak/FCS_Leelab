#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#%%
"""
Created on Thu Jun 1 11:32:34 2023
@author: HYUNHEE KWAK
Purpose: get autocorrelation curve from binned photon arrival data
        Paul Müller (2012) Python multiple-tau algorithm (Version 0.3.3). Available at https://pypi.python.org/pypi/multipletau/ (Accessed 01 JUN 2023)

 ---variables---
 bin_width_sec: float, binning width in second, desired input value 
 resolution: int, number of points on each level of lag-time, when to scale up
 
 ts_unit:       float, unit time of timestamp in second, fixed information from fluoroscence spectroscopy
 bin_width_n:   int, binning width in timestamp unit, calculated from above two
 
 binned_data:   array-like, binned count of photon arrival
             ts : array-like, photon arrival time stamps in integer number
             
 autocorr_array: array-like, autocorrelation value along multiple tau (lag time in second) in logarithmic scale
        m = even integer, resolution (number of points on one level, when to scale up)
        deltat = bin width in second for multiple tau autocorrelate
"""

#%%

# input values below 3

filename = "test.hdf5"
bin_width_sec = 200* 10**-9  #second
resolution = 8


#%%
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


ts = photon_data.timestamps.read()
ts_unit = photon_data.timestamps_specs.timestamps_unit.read() # ~50ns

binwidth_n = int(bin_width_sec / ts_unit) # binwidth in numbers : how many cycles per bin

binned_data = photon_bin(ts, binwidth_n).astype(np.float64)

autocorr_array = mt.autocorrelate(binned_data, m = resolution, deltat = bin_width_sec, normalize=True, copy=True, dtype=None, compress='average', ret_sum=False)



plt.rcParams.update({'figure.figsize':(7,5), 'figure.dpi':100})

# x-axis (lag-time in second) to Logarithmic scale 
plt.xscale("log")
plt.plot(autocorr_array[:,0], autocorr_array[:, 1])
# y-axis (extent of correlation)
plt.ylim(-0.2, 1)


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
     filename: Strings, filepath or filename in current directory
     bin_width_ns: int, binning width in nano second, desired input value 
     resolution: int, number of points on each level of lag-time, when to scale up
     time_gate_op: boolean, if time-gate option applied
     time_gate_valid: dic ('start': int, 'end': int), desired range of valid time in nano second when time_gate_op is True
 
 
 nanotimes: array-like, photon nanotimes in integer number
 tcspc_unit: float, unit time of tcspc in second
 tcspc_num_bins: int, split number between laser pulse
 tcspc_range: float, range in time between pulse
 bins: array-like, bin edges array from tcspc_num_bins for histogram
 
 ts:            array-like, photon arrival time stamps in integer number 
 ts_unit:       float, unit time of timestamp in second, fixed information from fluoroscence spectroscopy
 bin_width_n:   int, binning width in timestamp unit, calculated from above two
 binned_data:   array-like, binned count of photon arrival
       
 bin_width_sec: float, binning width in second, convert from bin_width_ns       
 autocorr_array: array-like, autocorrelation value along multiple tau (lag time in second) in logarithmic scale
        m = even integer, resolution (number of points on one level, when to scale up)
        deltat = bin width in second for multiple tau autocorrelate
"""

#%matplotlib inline #This line is for Jupyter Notebook environment
import os
import numpy as np
import math
import tables
import matplotlib.pyplot as plt
import multipletau as mt
from binning import photon_bin


#%% Tcspc histogram 

def tcspc_hist(nanotimes, tcspc_num_bins, tcspc_unit, start, end):
    
    plt.rcParams.update({'figure.figsize':(7,5), 'figure.dpi':100})
    fig, ax = plt.subplots()
    bins = np.arange(tcspc_num_bins)
    ax.hist(nanotimes, bins = bins)
    
    # x-axis (nanotimes in unit count) 
    ax.tick_params(top = True, labeltop = True, bottom = False, labelbottom = False, labelsize=15)
    font_x = {'size':15}
    ax.set_title('Time', fontdict = font_x) # used as Xaxis(top) label
    
    # second x-axis (nanotimes in nano-second)
    def nt2sec(x):
        return x * tcspc_unit * 1E9 # nano second
    def sec2nt(x):
        return x / tcspc_unit / 1E9
    
    secax = ax.secondary_xaxis('bottom', functions=(nt2sec, sec2nt))
    secax.set_xlabel('Time (ns)', fontdict = font_x)
    secax.xaxis.set_tick_params(labelsize=15)
    
    # y-axis (photone count)
    ax.set_ylabel('Counts in log scale', fontdict = font_x)
    ax.set_yscale("log")
    ax.yaxis.set_tick_params(labelsize=15)
    
    # valid range when time gated
    ax.axvspan(start, end, color='g', alpha=0.1) # color green


#%% Autocorrelation graph

def autocorr (ts, ts_unit, bin_width_sec, resolution):
    binwidth_n = int(bin_width_sec / ts_unit) # binwidth in numbers : how many cycles per bin
    binned_data = photon_bin(ts, binwidth_n).astype(np.float64)
    autocorr_array = mt.autocorrelate(binned_data, m = resolution, deltat = bin_width_sec, normalize=True, copy=True, dtype=None, compress='average', ret_sum=False)
    
    plt.rcParams.update({'figure.figsize':(7,5), 'figure.dpi':100})
    plt.figure()
    plt.plot(autocorr_array[:,0], autocorr_array[:, 1])
    
    # x-axis (lag-time in second) to Logarithmic scale 
    plt.xscale("log")
    font_x = {'size':15}
    plt.xlabel('lag time (s)', fontdict = font_x)
    plt.xticks(fontsize=15)
    
    # y-axis (extent of correlation)
    plt.ylim(0, 1)
    font_y = {'size':15}
    plt.ylabel('g(τ)', fontdict = font_y)
    plt.yticks(fontsize=15)
    
#%%

def main():
   # User input data

   filename = "test.hdf5"
   bin_width_ns = 200 #nano second
   resolution = 8
   time_gate_op = True
   time_gate_valid = {'start': 25, 'end': 50} # nano second


   #Data read
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
   
   if (time_gate_op):
       start = math.floor(time_gate_valid['start'] *1E-9 / tcspc_unit)  # covert to unit number
       end = math.ceil(time_gate_valid['end'] *1E-9 / tcspc_unit) # covert to unit number
       
   else:                    # whole range
       start = 0
       end = math.ceil(ts_unit/ tcspc_unit)
       
   
   tcspc_hist(nanotimes, tcspc_num_bins, tcspc_unit, start, end)
   
   if (time_gate_op):
       ts = ts[((nanotimes >= start) & (nanotimes <= end))]   # filter valid photons
    
   bin_width_sec = bin_width_ns * 1E-9  # convert unit from ns to second
   autocorr (ts, ts_unit, bin_width_sec, resolution)


if __name__ == "__main__":
    main()

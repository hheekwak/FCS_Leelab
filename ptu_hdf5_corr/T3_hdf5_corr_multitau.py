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
     time_gate_op: Strings, user input Y/N
         tg_op: boolean, if time-gate option applied
     time_gate_valid: dic ('start': int, 'end': int), desired range of valid time in nano second when time_gate_op is True
 
 
 nanotimes: array-like, photon nanotimes in integer number
 tcspc_unit: float, unit time of tcspc in second
 tcspc_num_bins: int, split number between laser pulse
 tcspc_range: float, range in time between pulse
 bins: array-like, bin edges array from tcspc_num_bins for histogram
 
 ts:            array-like, photon arrival time stamps in integer number 
 num_channel:
 ts_unit:       float, unit time of timestamp in second, fixed information from fluoroscence spectroscopy
 bin_width_n:   int, binning width in timestamp unit, calculated from above two
 binned_data:   array-like, binned count of photon arrival
       
 bin_width_sec: float, binning width in second, convert from bin_width_ns       
 autocorr_array: array-like, autocorrelation value along multiple tau (lag time in second) in logarithmic scale
        m = even integer, resolution (number of points on one level, when to scale up)
        deltat = bin width in second for multiple tau autocorrelate
"""

#%matplotlib inline #This line is for Jupyter Notebook environment
import csv
import numpy as np
import math
import tables
import matplotlib.pyplot as plt
import multipletau as mt
import user_input
from binning import photon_bin
from itertools import zip_longest




#%% Tcspc histogram 

def tcspc_hist(nanotimes, tcspc_num_bins, tcspc_unit, ch_detectors, channels, time_gate):
    
    plt.rcParams.update({'figure.figsize':(7,5), 'figure.dpi':100})
    fig, ax = plt.subplots()
    bins = np.arange(tcspc_num_bins)
    
    clr = 'g' # first channel span color green
    nanotimes_d = {}
    for ch in channels: 
        nanotimes_d['ch{0}'.format(ch)] = nanotimes[ch_detectors == ch]
        ax.hist(nanotimes_d['ch{0}'.format(ch)], bins = bins, color=clr)
        # valid range when time gated
        ax.axvspan(time_gate['ch{0}_st'.format(ch)], time_gate['ch{0}_end'.format(ch)], color=clr, alpha=0.1)
        clr = 'r' # second channel span color red
        
    
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
    
    return nanotimes_d


#%% Correlation graph

def corr(ts_gated_d, ts_unit, bin_width_sec, resolution, channels):
    """
    

    Parameters
    ----------
    ts_gated_d : TYPE
        DESCRIPTION.
    ts_unit : TYPE
        DESCRIPTION.
    bin_width_sec : TYPE
        DESCRIPTION.
    resolution : TYPE
        DESCRIPTION.
    channels : TYPE
        DESCRIPTION.

    Raises
    ------
    ValueError
        DESCRIPTION.

    Returns
    -------
    corr_array_d : TYPE
        DESCRIPTION.

    """
    
    plt.rcParams.update({'figure.figsize':(7,5), 'figure.dpi':100})
    plt.figure()
    
    binwidth_n = int(bin_width_sec / ts_unit) # binwidth in numbers : how many cycles per bin
    
    corr_array_d = {}
    if len(channels) == 1: 
        binned_data = photon_bin(ts_gated_d['ch{0}'.format(channels[0])], binwidth_n, ts_gated_d['ch{0}'.format(channels[0])][0], ts_gated_d['ch{0}'.format(channels[0])][-1] ).astype(np.float64)
        corr_array_d['auto_ch{0}'.format(channels[0])] = mt.autocorrelate(binned_data, m = resolution, deltat = bin_width_sec, normalize=True, copy=True, dtype=None, compress='average', ret_sum=False)
        plt.plot(corr_array_d['auto_ch{0}'.format(channels[0])][:, 0], corr_array_d['auto_ch{0}'.format(channels[0])][:, 1])
    
    elif len(channels) == 2:
        # set time range (start/end) of binning to match length for all 4 correlations, 
        start = min(ts_gated_d['ch{0}'.format(channels[0])][0], ts_gated_d['ch{0}'.format(channels[1])][0])
        end = max(ts_gated_d['ch{0}'.format(channels[0])][-1], ts_gated_d['ch{0}'.format(channels[1])][-1])
        
        binneddata_d = {} 
        clr = 'b'   # 1st ch autocorrelation graph color is blue
        for ch in channels: 
            binneddata_d['ch{0}'.format(ch)] = photon_bin(ts_gated_d['ch{0}'.format(ch)], binwidth_n, start, end).astype(np.float64)
            corr_array_d['auto_ch{0}'.format(ch)] = mt.autocorrelate(binneddata_d['ch{0}'.format(ch)], m = resolution, deltat = bin_width_sec, normalize=True, copy=True, dtype=None, compress='average', ret_sum=False)
            lbl = 'Ch{0}'.format(ch)
            plt.plot(corr_array_d['auto_ch{0}'.format(ch)][:,0], corr_array_d['auto_ch{0}'.format(ch)][:, 1], color = clr, label = lbl)
            clr = 'g' # 2nd ch autocorrelation graph color is green
            
        corr_array_d['cross_pos'] = mt.correlate(binneddata_d['ch{0}'.format(channels[0])], binneddata_d['ch{0}'.format(channels[1])], m = resolution, deltat = bin_width_sec, normalize=True, copy=True, dtype=None, compress="average", ret_sum=False)
        corr_array_d['cross_neg'] = mt.correlate(binneddata_d['ch{0}'.format(channels[1])], binneddata_d['ch{0}'.format(channels[0])], m = resolution, deltat = bin_width_sec, normalize=True, copy=True, dtype=None, compress="average", ret_sum=False)
        
        plt.plot(corr_array_d['cross_pos'][:,0], corr_array_d['cross_pos'][:, 1], color = 'k', label = 'Pos.Cross-corr.') # positive crosscorrelation color is black
        plt.plot(corr_array_d['cross_neg'][:,0], corr_array_d['cross_neg'][:, 1], color = 'm', label = 'Neg.Cross-corr.') # negative crosscorrelation color is magenta
        
    else:
        raise ValueError('ATTENTION: more than 2 channels detected. Check detectors data.\n')
    
    
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
    
    plt.legend()

    return corr_array_d
    
#%%

def main():
    # User input data
    filename = input("Enter file name: ") # test.hdf5
    a = user_input.UserInput(filename)
    b_w = input("Enter bin width in ns (default = 200 ns) : ") # nano second
    
    if b_w !="":
        a.bin_width_ns = int(b_w)
        
    res = input("Enter desired resolution in int (default = 8): ")    
    if res !="":
        a.resolution = res 
 
    tg_confirm = False
    while (tg_confirm == False):
        tg_op = input("Will you use time gating? (Y/N): ")
        if tg_op == "Y" or tg_op =="y":
            a.time_gate_op = True
            tg_confirm = True
        elif tg_op == "N" or  tg_op == "n":
            a.time_gate_op = False
            tg_confirm = True
        else:
            pass
    
    h5file = tables.open_file(a.filename)
    photon_data = h5file.root.photon_data
    photon_data.measurement_specs.measurement_type.read().decode()
 
    nanotimes = photon_data.nanotimes.read()
    tcspc_unit = photon_data.nanotimes_specs.tcspc_unit.read()
    tcspc_num_bins = photon_data.nanotimes_specs.tcspc_num_bins.read()
    tcspc_range = photon_data.nanotimes_specs.tcspc_range.read()
 
    ts = photon_data.timestamps.read()
    ch_detectors = photon_data.detectors.read()
    channels = np.unique(ch_detectors)
    ts_unit = photon_data.timestamps_specs.timestamps_unit.read() # ~50ns
 
    tg = {}
    if a.time_gate_op == True:
        for ch in channels: 
            tg['ch{0}_st'.format(ch)] = math.floor(int(input ("time gate START time of CH{0} in nano second: ".format(ch)))*1E-9 / tcspc_unit) # convert to count number
            tg['ch{0}_end'.format(ch)] = math.ceil(int(input ("time gate END time of CH{0} in nano second: ".format(ch)))*1E-9 / tcspc_unit) # convert to count number
    else:
        for ch in channels:  
            tg['ch{0}_st'.format(ch)] = 0
            tg['ch{0}_end'.format(ch)] = math.ceil(ts_unit/ tcspc_unit)
    
    nanotimes_d = tcspc_hist(nanotimes, tcspc_num_bins, tcspc_unit, ch_detectors, channels, tg)

    # split each channel's timestamps and apply timegating
    ts_gated_d = {} # gated timestamps 
    for ch in channels: 
        ts_gated_d['ch{0}'.format(ch)] = ts[ch_detectors == ch] 
        ts_gated_d['ch{0}'.format(ch)] = ts_gated_d['ch{0}'.format(ch)][(nanotimes_d['ch{0}'.format(ch)] >= tg['ch{0}_st'.format(ch)]) & (nanotimes_d['ch{0}'.format(ch)] <= tg['ch{0}_end'.format(ch)])]
 
    bin_width_sec = a.bin_width_ns * 1E-9  # convert unit from ns to second
    correlation_d = corr(ts_gated_d, ts_unit, bin_width_sec, a.resolution, channels)
   
    # Data write
    # TODO: write header file

    
    # write timestamps file
    ts_d = [ch_detectors, ts]
    export_data_ts = zip_longest(*ts_d, fillvalue = '')
    with open(a.filename[:-5] + '_timestamps.csv', 'w', newline='') as csvfile:
        wr = csv.writer(csvfile)
        wr.writerow(('Channel', 'Time Tag'))
        wr.writerows(export_data_ts)
    csvfile.close()
            
    # write correlation file
    with open(a.filename[:-5] + '_correlations.csv', 'w', newline='') as csvfile:
        wr = csv.writer(csvfile) 
        
        if len(channels) == 1:
            wr.writerow(('Tau in s', 'G(A,A)'))
            wr.writerows(correlation_d['auto_ch{0}'.format(channels[0])])
        
        else: 
            corr_d = [correlation_d['auto_ch{0}'.format(channels[0])][:,0], correlation_d['auto_ch{0}'.format(channels[0])][:,1], correlation_d['auto_ch{0}'.format(channels[1])][:,1], correlation_d['cross_pos'][:,1], correlation_d['cross_neg'][:,1]]
            export_data_cor = zip_longest(*corr_d, fillvalue = '')
        
            wr.writerow(('Tau in s', 'G(A,A)', 'G(B,B)', 'G(A,B)', 'G(B,A)'))
            wr.writerows(export_data_cor)
    
    csvfile.close()
    

if __name__ == "__main__":
    main()

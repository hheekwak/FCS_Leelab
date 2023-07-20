#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#%%
"""
Created on Wed Jul 12 11:35:53 2023
@author: HYUNHEE KWAK
Purpose: Get correlation curves and data from photon arrival data file *.ptu
        Paul Müller (2012) Python multiple-tau algorithm (Version 0.3.3). Available at https://pypi.python.org/pypi/multipletau/ (Accessed 01 JUN 2023)


Below 5 variables are user input variables 
     filename:      Strings, filepath or filename in current directory
     bin_width_ns:  int, binning width in nano second, desired input value (default = 200 ns) 
     resolution:    int, (even) number of points on each level of lag-time, when to scale up
     time_gate_op:  Strings, user input Y/N
     time_gate_valid: dic ('start': int, 'end': int), desired range of valid time in nano second when time_gate_op is True
"""

import os
import csv
import numpy as np
import math
import matplotlib.pyplot as plt
import multipletau as mt
import user_input
import ptu_reader
from binning import photon_bin
from itertools import zip_longest

#%% Tcspc histogram 

def tcspc_hist(nanotimes, tcspc_num_bins, tcspc_unit, ch_detectors, channels, time_gate):
    '''

    Parameters
    ----------
    nanotimes :         array-like, photon nanotimes in integer number
    tcspc_num_bins :    int, split number between laser pulse
    tcspc_unit :        float, unit time of tcspc in second
    ch_detectors :      array-like, detected channels of each photon
    channels :          array-like, list of existing channels
    time_gate :         dict ('start': int, 'end': int), each channel's desired start and end of valid time in nano second

    Returns
    -------
    nanotimes_d:        dict ('ch_': array-like), each channel's nanotimes split
    plots histogram

    '''
    
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
    '''
    

    Parameters
    ----------
    ts_gated_d :    dict ('ch_': array-like), photon arrival time stamps in integer number of each channel (valid nanotime gated if needed) 
    ts_unit :       float, unit time of timestamp in second, fixed information from fluoroscence spectroscopy
    bin_width_sec : float, binning width in second
    resolution :    int, (even) number of points on each level of lag-time, when to scale up
                        Please read "# Issue_1" in main() for more information
    channels :      array-like, list of existing channels

    Raises
    ------
    ValueError
        maximum count of channels is 2. when more than 2 channels detected, error message raises

    Returns
    -------
    corr_array_d :  dict ('auto_ch_': array-like, 'cross_pos': array-like, 'cross_neg': array-like), 
                        value of autocorrelation (A,A), (B,B) and cross-correlation (A,B), (B,A)

    '''
    
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
    plt.ylim(0 ,1)
    font_y = {'size':15}
    plt.ylabel('g(τ)', fontdict = font_y)
    plt.yticks(fontsize=15)
    
    plt.legend()

    return corr_array_d
    
#%%

def main():
    
    # User input data
    filename = 'no_filename'
    while (not filename[-4:] == '.ptu'):
        filename = input("Enter file name (*.ptu): ") # test.ptu
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
    
    bin_width_sec = a.bin_width_ns * 1E-9  # convert unit from ns to second
    d, meta, tags = ptu_reader.read_ptu(a.filename)
    ts =  d['photon_data']['timestamps']
    ch_detectors = d['photon_data']['detectors']
    ts_unit = d['photon_data']['timestamps_specs']['timestamps_unit']
    nanotimes = d['photon_data']['nanotimes']
    channels = np.unique(ch_detectors)
    
    # split each channel's timestamps
    ts_gated_d = {}
    for ch in channels: 
        ts_gated_d['ch{0}'.format(ch)] = ts[ch_detectors == ch] 
    
    if meta['record_type'][-2:] == 'T2':
        pass
    else:
        tcspc_unit = d['photon_data']['nanotimes_specs']['tcspc_unit']
        tcspc_num_bins = d['photon_data']['nanotimes_specs']['tcspc_num_bins']
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

        # filter gated timestamps
        for ch in channels: 
            ts_gated_d['ch{0}'.format(ch)] = ts_gated_d['ch{0}'.format(ch)][(nanotimes_d['ch{0}'.format(ch)] >= tg['ch{0}_st'.format(ch)]) & (nanotimes_d['ch{0}'.format(ch)] <= tg['ch{0}_end'.format(ch)])]
    
    #-	Issue_1: multiple tau correlation function results in data point interval length of resolution/2 but not resolution from the second level 
    #           e.g.) when m = 8, 8-4-4-4-4-4-, when m = 16, 16-8-8-8-8-8- 
    #           To have the data list in desired data point length for each level, used 2*a.resolution as resolution (m)

    correlation_d = corr(ts_gated_d, ts_unit, bin_width_sec, 2*a.resolution, channels)
        
    # write files
    path_split = os.path.split(a.filename) 
    path_dir = path_split[0]
    name = path_split[1][:-4]
    file_dir = os.path.join(path_dir, name)
    filepath = os.path.join(file_dir, name)
    
    # open a directory to save corresponing files 
    if not os.path.exists(file_dir):
        os.mkdir(file_dir)
        # header file write
        ptu_reader.ptu_header_file(tags, filepath)
        # timestamps file write
        ptu_reader.ptu_ts_file(d, filepath)
    
    else:
        # check if header file and timestamps file exist in the directory, write those files if needed
        headerfile = filepath + '_header.txt'
        if not os.path.exists(headerfile):
            ptu_reader.ptu_header_file(tags, filepath)
            
        ts_file = filepath + '_timestamps.csv'
        if not os.path.exists(ts_file):
            ptu_reader.ptu_ts_file(d, filepath)
   
    # write correlation file
    with open(filepath + '_correlations.csv', 'w', newline='') as csvfile:
        wr = csv.writer(csvfile) 
        
        if len(channels) == 1:
            wr.writerow(('Tau in s', 'G(A,A)'))
            wr.writerows(correlation_d['auto_ch{0}'.format(channels[0])])
        
        else: 
            corr_d = [correlation_d['auto_ch{0}'.format(channels[0])][:,0], correlation_d['auto_ch{0}'.format(channels[0])][:,1], correlation_d['auto_ch{0}'.format(channels[1])][:,1], correlation_d['cross_pos'][:,1], correlation_d['cross_neg'][:,1]]
            export_data_cor = zip_longest(*corr_d, fillvalue = '')
        
            wr.writerow(('Tau in s', 'G(A,A)', 'G(B,B)', 'G(A,B)', 'G(B,A)'))
            wr.writerows(export_data_cor)
    
if __name__ == "__main__":
    main()

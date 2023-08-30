#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#%%
"""
Created on Thu Jul 20 10:46:50 2023

@author: HYUNHEE KWAK
Purpose: Graphical User Interface of ptu_corr_multitau
        Get correlation curves and data from photon arrival data file *.ptu
        For more information about this program, please see ptu_corr_multitau.ptu 

"""
import os
import sys
from PyQt5.QtWidgets import QApplication, QFileDialog, QMainWindow
from PyQt5 import uic
import csv
import numpy as np
import math
import matplotlib.pyplot as plt
import multipletau as mt
import ptu_reader
from itertools import zip_longest
import threading
from queue import Queue
import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas 

#%% photon_bin() function from binning.py

def photon_bin(ts, bin_width, start, end):
    """
    

    Parameters
    ----------
    ts :  array-like,
        Phtone arraival time stamps in integer number
    bin_width : float,
        Desired bin width in timestamp unit.
    start : int
        starting time to bin
    end : int
        ending time of binning

    Returns
    -------
    binned_data : array-like, 
                binned count of photon arrival.

    """

    # Calculate the number of bins
    num_bins = int(np.ceil((end - start + 1) / bin_width))

    # Create an array to store the binned data
    binned_data = np.zeros(num_bins, dtype=int)

    # Iterate through the data points and increment the corresponding bin
    for point in ts:
        bin_index = (point - start) // bin_width
        binned_data[bin_index] += 1

    return binned_data
        


#%% Tcspc histogram without graph

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

    '''
    
    nanotimes_d = {}
    for ch in channels: 
        nanotimes_d['ch{0}'.format(ch)] = nanotimes[ch_detectors == ch]
    
    return nanotimes_d

#%%
# connect ui file 
form_class = uic.loadUiType("ptu_corr_gui_frame.ui")[0]

class WindowClass(QMainWindow, form_class) :
    def __init__(self) :
        super().__init__()
        self.setupUi(self)
        self.filename = 'no_filename'
        self.bin_width_ns = 200
        self.resolution = 8
        self.time_gate_op = False
        self.tg = {}
        # Convert unit from ns to second
        self.bin_width_sec = self.bin_width_ns * 1E-9 
    
        # To save all corr data from files in one file in case of batchRun
        self.batchRun = False
        self.data_dic_by_file = {}
        
        # Graph
        self.fig = plt.Figure()
        self.canvas = FigureCanvas(self.fig)
        self.graph_verticalLayout.addWidget(self.canvas)
        
        # Individual Process Mode
        self.browseButton.clicked.connect(self.browsefiles)
        self.binwidth_line.textChanged.connect(self.setbincheck)
        self.binwidth_line.returnPressed.connect(self.setbinwidth)
        self.resolution_line.textChanged.connect(self.setrescheck)
        self.resolution_line.returnPressed.connect(self.setresolution)
        
        self.tg_groupBox.clicked.connect(self.set_tg_op)
        self.timegate_slider_0.valueChanged.connect(self.set_tg_range_0)
        self.timegate_slider_1.valueChanged.connect(self.set_tg_range_1)
        self.timegate_slider_2.valueChanged.connect(self.set_tg_range_2)
        self.timegate_slider_3.valueChanged.connect(self.set_tg_range_3)
        
        self.widget_list = [self.filename_label, self.filename_line, self.browseButton, self.binwidth_label, self.binwidth_line, self.resolution_label, self.res_expl_label, self.resolution_line, self.tg_groupBox, self.run_pushButton]
        
        self.run_pushButton.clicked.connect(self.onIndividualRun)
        self.cancel_pushButton.clicked.connect(self.onCancelClick)
        
        # Batch Process Mode
        self.browseButton_2.clicked.connect(self.browsedirectory)
        self.run_pushButton_2.clicked.connect(self.batch_process)
        self.cancel_pushButton_2.clicked.connect(self.onCancelClick)
        
    def browsefiles(self):
        fname = QFileDialog.getOpenFileName(self, 'Open file', os.getcwd() , 'PTU files (*.ptu)')
        self.filename_line.setText(fname[0])
        self.filename = fname[0]
        if (self.preparation()):
            self.run_pushButton.setEnabled(True)
        
    def browsedirectory(self):
        self.dname = str(QFileDialog.getExistingDirectory(self, 'Select Directory'))
        self.dirpath_line.setText(self.dname)
        self.run_pushButton_2.setEnabled(True)
         
    def batch_process(self):
        for root, dirs, files in os.walk(self.dname):
            for file in files:
                # Check if the file extension is '.ptu'
                if file.endswith('.ptu'):
                    # Create the filepath by joining the root directory and file name
                    self.filename = (os.path.join(root, file))
                    print(self.filename)
                    if (self.preparation()):
                        self.onBatchRun()                
        print("Batch process completed. ")
                    
    def preparation(self):
        self.d, self.meta, self.tags = ptu_reader.read_ptu(self.filename)
        self.ts =  self.d['photon_data']['timestamps']
        self.ch_detectors = self.d['photon_data']['detectors']
        self.ts_unit = self.d['photon_data']['timestamps_specs']['timestamps_unit']
        self.nanotimes = self.d['photon_data']['nanotimes']
        self.channels = np.unique(self.ch_detectors)
        self.nanotime_max = math.ceil(self.ts_unit *1E9)
        st_nt_max = str(self.nanotime_max)
        
        # If T2 mode, disable TCSPC data process 
        if self.meta['record_type'][-2:] == 'T2':
            self.tg_groupBox.setChecked(False)
        else:
            self.tcspc_unit = self.d['photon_data']['nanotimes_specs']['tcspc_unit']
            self.tcspc_num_bins = self.d['photon_data']['nanotimes_specs']['tcspc_num_bins']
            # TCSPC channel labels and range set
            self.timegate_nt_0.setText(st_nt_max + " ns")
            self.timegate_nt_1.setText(st_nt_max + " ns")
            self.timegate_nt_2.setText(st_nt_max + " ns")
            self.timegate_nt_3.setText(st_nt_max + " ns")         
            self.timegate_slider_0.setRange(0, self.nanotime_max)
            self.timegate_slider_1.setRange(0, self.nanotime_max)
            self.timegate_slider_2.setRange(0, self.nanotime_max)
            self.timegate_slider_3.setRange(0, self.nanotime_max)
            self.timegate_slider_0.setValue(0)
            self.timegate_slider_1.setValue(self.nanotime_max)
            self.tg['ch{0}_st'.format(self.channels[0])] = self.timegate_slider_0.value()*1E-9 / self.tcspc_unit # convert to count number
            self.tg['ch{0}_end'.format(self.channels[0])] = self.timegate_slider_1.value()*1E-9 / self.tcspc_unit # convert to count number
            self.timegate_slider_2.setValue(0)
            self.timegate_slider_3.setValue(self.nanotime_max)
            if len(self.channels) > 1:
                self.tg['ch{0}_st'.format(self.channels[1])] = self.timegate_slider_2.value()*1E-9 / self.tcspc_unit # convert to count number
                self.tg['ch{0}_end'.format(self.channels[1])] = self.timegate_slider_3.value()*1E-9 / self.tcspc_unit # convert to count number    
        
        return True
        
    def setbincheck(self):
        if not self.binwidth_line.text().isnumeric() :
            self.invalid_label.setText("Invalid")
        else:
            self.invalid_label.setText("")
            self.setbinwidth()
        
    def setbinwidth(self):
        if self.binwidth_line.text() == "" or self.binwidth_line.text() == "0" :
            self.binwidth_label.setText("Current Bin width is set as default = 200 ns")
        else:
            self.bin_width_ns = int(self.binwidth_line.text())
            self.bin_width_sec = self.bin_width_ns * 1E-9 
            self.binwidth_label.setText("Current Bin width is set as " + self.binwidth_line.text() + " ns")
        self.resolution_line.setCursorPosition(0)
            
    def setrescheck(self):
        if not self.resolution_line.text().isnumeric() :
            self.invalid_label_2.setText("Invalid")
        else:
            self.invalid_label_2.setText("")
            self.setresolution()
            
    def setresolution(self):
        if self.resolution_line.text() == "" or self.resolution_line.text() == "0" :
            self.resolution_label.setText("Current Resolution is set as default = 8")
            self.res_expl_label.setText("* Data points length at each level will be as 2n-n-n-n-n-n-n-")
        else:
            n_st = self.resolution_line.text()
            n2_st = str(2 * int(n_st))
            self.resolution = int(n_st)
            self.resolution_label.setText("Current Resolution is set as " + n_st)
            self.res_expl_label.setText("* Data points length at each level will be as " + n2_st + "-" + n_st + "-" + n_st + "-" + n_st + "-" + n_st + "-" + n_st + "-" + n_st + "-")

    def set_tg_op(self):
        if self.meta['record_type'][-2:] == 'T2':
            self.tg_groupBox.setChecked(False)
            self.time_gate_op = False
            print("time gating disabled in T2 mode")
        else:
            self.time_gate_op = self.tg_groupBox.isChecked()
            # once activated, range initialized to min, max value
            if (self.time_gate_op):
                self.timegate_slider_0.setValue(0)
                self.tg['ch{0}_st'.format(self.channels[0])] = self.timegate_slider_0.value()*1E-9 / self.tcspc_unit # convert to count number
                self.timegate_range_0s_label.setText('Ch{0} Start: '.format(self.channels[0]) + str(self.timegate_slider_0.value()) + " ns")
                self.timegate_slider_1.setValue(self.nanotime_max)
                self.tg['ch{0}_end'.format(self.channels[0])] = self.timegate_slider_1.value()*1E-9 / self.tcspc_unit # convert to count number
                self.timegate_range_0e_label.setText('Ch{0} End: '.format(self.channels[0]) + str(self.timegate_slider_1.value()) + " ns")
                # If two channels, set second channel range
                if len(self.channels) > 1:
                    self.timegate_slider_2.setValue(0)
                    self.tg['ch{0}_st'.format(self.channels[1])] = self.timegate_slider_2.value()*1E-9 / self.tcspc_unit # convert to count number
                    self.timegate_range_1s_label.setText('Ch{0} Start: '.format(self.channels[1]) + str(self.timegate_slider_2.value()) + " ns")
                    self.timegate_slider_3.setValue(self.nanotime_max)
                    self.tg['ch{0}_end'.format(self.channels[1])] = self.timegate_slider_3.value(*1E-9) / self.tcspc_unit # convert to count number
                    self.timegate_range_1e_label.setText('Ch{0} End: '.format(self.channels[1]) + str(self.timegate_slider_3.value()) + " ns")
        
    def set_tg_range_0(self):
        if self.timegate_slider_0.value() >= self.timegate_slider_1.value():
            self.invalid_label_3.setText("Invalid")
        else:
            self.invalid_label_3.setText("")
            self.invalid_label_4.setText("")
        self.tg['ch{0}_st'.format(self.channels[0])] = self.timegate_slider_0.value()*1E-9 / self.tcspc_unit # convert to count number
        self.timegate_range_0s_label.setText('Ch{0} Start: '.format(self.channels[0]) + str(self.timegate_slider_0.value()) + " ns")
    def set_tg_range_1(self):
        if self.timegate_slider_1.value() <= self.timegate_slider_0.value():
            self.invalid_label_4.setText("Invalid")
        else:
            self.invalid_label_4.setText("")
            self.invalid_label_3.setText("")
        self.tg['ch{0}_end'.format(self.channels[0])] = self.timegate_slider_1.value()*1E-9 / self.tcspc_unit # convert to count number
        self.timegate_range_0e_label.setText('Ch{0} End: '.format(self.channels[0]) + str(self.timegate_slider_1.value()) + " ns")
    def set_tg_range_2(self):
        if self.timegate_slider_2.value() >= self.timegate_slider_3.value():
            self.invalid_label_5.setText("Invalid")
        else:
            self.invalid_label_5.setText("")
            self.invalid_label_6.setText("")
        if len(self.channels) > 1:
            self.tg['ch{0}_st'.format(self.channels[1])] = self.timegate_slider_2.value()*1E-9 / self.tcspc_unit # convert to count number
            self.timegate_range_1s_label.setText('Ch{0} Start: '.format(self.channels[1]) + str(self.timegate_slider_2.value()) + " ns")
    def set_tg_range_3(self):
        if self.timegate_slider_3.value() <= self.timegate_slider_2.value():
            self.invalid_label_6.setText("Invalid")
        else:
            self.invalid_label_6.setText("")
            self.invalid_label_5.setText("")
        if len(self.channels) > 1:
            self.tg['ch{0}_end'.format(self.channels[1])] = self.timegate_slider_3.value()*1E-9 / self.tcspc_unit # convert to count number
            self.timegate_range_1e_label.setText('Ch{0} End: '.format(self.channels[1]) + str(self.timegate_slider_3.value()) + " ns")
    
    def onIndividualRun(self):
        self.batchRun = False
        self.onStart()
        print("Individual process completed. ")
        
    def onBatchRun(self):
        self.batchRun = True
        self.onStart()
    
    def onStart(self):
        # Clear the plot before new work
        self.fig.clear()
        
        
        # Widgets enabled/disabled during the process
        self.cancel_flag = threading.Event()
        self.cancel_pushButton.setEnabled(True)
        self.cancel_pushButton_2.setEnabled(True)
        for widget in self.widget_list:
            widget.setEnabled(False)    
    
        print("Working on " + self.filename)
        # Check if the cancel_flag is set
        if self.cancel_flag.is_set():
            return
        print("step 1")
        self.ch_time_gating()
       
        if self.cancel_flag.is_set():
            return
        print("step 2")
        if self.meta['record_type'][-2:] != 'T2':
            self.tcspc_draw()
       
        if self.cancel_flag.is_set():
            return
        print("step 3")
        self.run_corr()
        
        print("step 4")
        if self.cancel_flag.is_set():
            return
        self.corr_draw()
        
        if self.cancel_flag.is_set():
            return
        print("step 5")
        self.write_files()
        
        self.fig.subplots_adjust(left=0.15, right=0.95, top=0.90, bottom=0.1, hspace=0.5)
        self.canvas.draw()
        
        
    def ch_time_gating(self):
        # split each channel's timestamps
        self.ts_gated_d = {}
        for ch in self.channels: 
            self.ts_gated_d['ch{0}'.format(ch)] = self.ts[self.ch_detectors == ch]
            
        if self.meta['record_type'][-2:] == 'T2':
            pass
        else:
            self.nanotimes_d = tcspc_hist(self.nanotimes, self.tcspc_num_bins, self.tcspc_unit, self.ch_detectors, self.channels, self.tg)
            # filter gated timestamps
            for ch in self.channels: 
                self.ts_gated_d['ch{0}'.format(ch)] = self.ts_gated_d['ch{0}'.format(ch)][(self.nanotimes_d['ch{0}'.format(ch)] >= self.tg['ch{0}_st'.format(ch)]) & (self.nanotimes_d['ch{0}'.format(ch)] <= self.tg['ch{0}_end'.format(ch)])]

    def run_corr(self):
        print("Correlation step takes a while... ")
        #-	Issue_1: multiple tau correlation function results in data point interval length of resolution/2 but not resolution from the second level 
        #           e.g.) when m = 8, 8-4-4-4-4-4-, when m = 16, 16-8-8-8-8-8- 
        #           To have the data list in desired data point length for each level, used 2*a.resolution as resolution (m)
        self.corr_array_d  = {}
        
        binwidth_n = int(self.bin_width_sec / self.ts_unit) # binwidth in numbers : how many cycles per bin
        
        if len(self.channels) == 1: 
            binned_data = photon_bin(self.ts_gated_d['ch{0}'.format(self.channels[0])], binwidth_n, self.ts_gated_d['ch{0}'.format(self.channels[0])][0], self.ts_gated_d['ch{0}'.format(self.channels[0])][-1] ).astype(np.float64)
            self.corr_array_d['auto_ch{0}'.format(self.channels[0])] = mt.autocorrelate(binned_data, m = 2*self.resolution, deltat = self.bin_width_sec, normalize=True, copy=True, dtype=None, compress='average', ret_sum=False)
        
        elif len(self.channels) == 2:
            # set time range (start/end) of binning to match length for all 4 correlations, 
            start = min(self.ts_gated_d['ch{0}'.format(self.channels[0])][0], self.ts_gated_d['ch{0}'.format(self.channels[1])][0])
            end = max(self.ts_gated_d['ch{0}'.format(self.channels[0])][-1], self.ts_gated_d['ch{0}'.format(self.channels[1])][-1])
            
            binneddata_d = {} 
            for ch in self.channels: 
                binneddata_d['ch{0}'.format(ch)] = photon_bin(self.ts_gated_d['ch{0}'.format(ch)], binwidth_n, start, end).astype(np.float64)
                self.corr_array_d['auto_ch{0}'.format(ch)] = mt.autocorrelate(binneddata_d['ch{0}'.format(ch)], m = 2*self.resolution, deltat = self.bin_width_sec, normalize=True, copy=True, dtype=None, compress='average', ret_sum=False)
               
            self.corr_array_d['cross_pos'] = mt.correlate(binneddata_d['ch{0}'.format(self.channels[0])], binneddata_d['ch{0}'.format(self.channels[1])], m = 2*self.resolution, deltat = self.bin_width_sec, normalize=True, copy=True, dtype=None, compress="average", ret_sum=False)
            self.corr_array_d['cross_neg'] = mt.correlate(binneddata_d['ch{0}'.format(self.channels[1])], binneddata_d['ch{0}'.format(self.channels[0])], m = 2*self.resolution, deltat = self.bin_width_sec, normalize=True, copy=True, dtype=None, compress="average", ret_sum=False)
                   
        else:
            raise ValueError('ATTENTION: more than 2 channels detected. Check detectors data.\n')
        
    def write_files(self):
        print("Writing files...")
        # write files
        path_split = os.path.split(self.filename) 
        path_dir = path_split[0]
        name = path_split[1][:-4]
        file_dir = os.path.join(path_dir, name)
        filepath = os.path.join(file_dir, name)
        
        # open a directory to save corresponing files 
        if not os.path.exists(file_dir):
            os.mkdir(file_dir)
            # header file write
            ptu_reader.ptu_header_file(self.tags, filepath)
            # timestamps file write
            ptu_reader.ptu_ts_file(self.d, filepath)
        
        else:
            # check if header file and timestamps file exist in the directory, write those files if needed
            headerfile = filepath + '_header.txt'
            if not os.path.exists(headerfile):
                ptu_reader.ptu_header_file(self.tags, filepath)
                
            ts_file = filepath + '_timestamps.csv'
            if not os.path.exists(ts_file):
                ptu_reader.ptu_ts_file(self.d, filepath)
       
        # write correlation file
        with open(filepath + '_correlations.csv', 'w', newline='') as csvfile:
            wr = csv.writer(csvfile) 
            
            if len(self.channels) == 1:
                wr.writerow(('Tau in s', 'G(A,A)'))
                wr.writerows(self.corr_array_d['auto_ch{0}'.format(self.channels[0])])
            
            else: 
                corr_d = [self.corr_array_d['auto_ch{0}'.format(self.channels[0])][:,0], self.corr_array_d['auto_ch{0}'.format(self.channels[0])][:,1], self.corr_array_d['auto_ch{0}'.format(self.channels[1])][:,1], self.corr_array_d['cross_pos'][:,1], self.corr_array_d['cross_neg'][:,1]]
                export_data_cor = zip_longest(*corr_d, fillvalue = '')
            
                wr.writerow(('Tau in s', 'G(A,A)', 'G(B,B)', 'G(A,B)', 'G(B,A)'))
                wr.writerows(export_data_cor)
        
    def onCancelClick(self):
        # Called when the cancel_pushButton is clicked
        # Set the cancel_flag to signal the thread to stop
        self.cancel_flag.set()
        print("Process stopped.")
        for widget in self.widget_list:
            widget.setEnabled(True)
      
    def tcspc_draw(self):
        ax = self.fig.add_subplot(211)
        font_size = 12
        
        bins = np.arange(self.tcspc_num_bins)
        clr = 'g' # first channel span color green
        for ch in self.channels: 
            ax.hist(self.nanotimes_d['ch{0}'.format(ch)], bins = bins, color=clr)
            # valid range when time gated
            ax.axvspan(self.tg['ch{0}_st'.format(ch)], self.tg['ch{0}_end'.format(ch)], color=clr, alpha=0.1)
            clr = 'r' # second channel span color red
            
        # x-axis (nanotimes in unit count) 
        ax.tick_params(top = True, labeltop = True, bottom = False, labelbottom = False, labelsize=font_size)
        font_x = {'size':font_size}
        ax.set_title('Time', fontdict = font_x) # used as Xaxis(top) label
        
        # second x-axis (nanotimes in nano-second)
        def nt2sec(x):
            return x * self.tcspc_unit * 1E9 # nano second
        def sec2nt(x):
            return x / self.tcspc_unit / 1E9
        
        secax = ax.secondary_xaxis('bottom', functions=(nt2sec, sec2nt))
        secax.set_xlabel('Time (ns)', fontdict = font_x)
        secax.xaxis.set_tick_params(labelsize=font_size)
        
        # y-axis (photone count)
        ax.set_ylabel('Counts in log scale', fontdict = font_x)
        ax.set_yscale("log")
        ax.yaxis.set_tick_params(labelsize=font_size)
        
    def corr_draw(self):
        ax2 = self.fig.add_subplot(212)
        font_size = 12
        
        if len(self.channels) == 1: 
            ax2.plot(self.corr_array_d['auto_ch{0}'.format(self.channels[0])][:, 0], self.corr_array_d['auto_ch{0}'.format(self.channels[0])][:, 1], label = 'autocorrelation')
        
        elif len(self.channels) == 2:
            clr = 'b'   # 1st ch autocorrelation graph color is blue
            for ch in self.channels: 
                lbl = 'Ch{0}'.format(ch)
                ax2.plot(self.corr_array_d['auto_ch{0}'.format(ch)][:,0], self.corr_array_d['auto_ch{0}'.format(ch)][:, 1], color = clr, label = lbl)
                clr = 'g' # 2nd ch autocorrelation graph color is green
                
            ax2.plot(self.corr_array_d['cross_pos'][:,0], self.corr_array_d['cross_pos'][:, 1], color = 'k', label = 'Pos.Cross-corr.') # positive crosscorrelation color is black
            ax2.plot(self.corr_array_d['cross_neg'][:,0], self.corr_array_d['cross_neg'][:, 1], color = 'm', label = 'Neg.Cross-corr.') # negative crosscorrelation color is magenta
            
        else:
            raise ValueError('ATTENTION: more than 2 channels detected. Check detectors data.\n')
        
        # x-axis (lag-time in second) to Logarithmic scale 
        ax2.set_xscale("log")
        font_x = {'size':font_size}
        ax2.set_xlabel('lag time (s)', fontdict = font_x)
        ax2.xaxis.set_tick_params(labelsize=font_size)
        ax2.xaxis.set_major_locator(matplotlib.ticker.LogLocator(numticks=999, subs=(.1,.2,.3,.4,.5,.6,.7,.8,.9)))

        
        # y-axis (extent of correlation)
        ax2.set_ylim(0 ,1)
        font_y = {'size':font_size}
        ax2.set_ylabel('g(τ)', fontdict = font_y)
        ax2.yaxis.set_tick_params(labelsize=font_size)
        ax2.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
        ax2.legend()


        
if __name__ == "__main__" :
    app = QApplication(sys.argv) 
    mainWindow = WindowClass() 
    mainWindow.show()
    sys.exit(app.exec_())

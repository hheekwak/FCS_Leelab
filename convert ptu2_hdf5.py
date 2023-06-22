#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  6 14:51:34 2023

@author: YOUNGKWANG
Purpose: convert PicoQuant photoncounting file to HDF5 format.
Version: 0p2
The code is modified from v0p1 to process all ptu files within a specified top directory and its subfolders. 
The files include those in subdirectories.
"""
# Assign variables for author, author affiliation, description, sample name, dye names, and buffer name.
author = 'Youngkwang Lee'
author_affiliation = 'SDSU'
description = 'time-resolved fluorescence spectroscopy'
Experiment_ID = 'Test data' # Enter experiment code, eg. FCS15
dye_names = 'ATTO488, '
buffer_name = 'PBS pH 7.4'

# Assign the directory path where the data files are stored
directory = "G:\My Drive\Hyunhee work_local\FCS_leelab" # Enter the top directory

#%%
import os

if os.path.exists(directory):
    print('Directory found, it will be proceed.')
else:
    print('ATTENTION: Directory not found, please check the directory name.\n'
          '           (current value "%s")' % directory)
#%%
#%matplotlib inline # This line will work only in a Jupyter Notebook environment
# Import the necessary packages
import numpy as np
import phconvert as phc
import os

# Walk through the directory tree and iterate over each file
for root, dirs, files in os.walk(directory):
    for file in files:
        # Check if the file extension is '.ptu'
        if file.endswith('.ptu'):
            # Create the filepath by joining the root directory and file name
            filepath = os.path.join(root, file)
            
            # Print the version of phconvert package
            print('phconvert version: ' + phc.__version__)
            
            # Load the data from the file using nsalex_pq function from phconvert package
            # Here, we specify various parameters such as donor, acceptor, excitation and detection wavelengths etc.
            d, meta = phc.loader.nsalex_pq(filepath,
                                           donor = 1, #YKL, always donor uses channel 1
                                           acceptor = 0, #YKL, If we have only two channels available (0 and 1), channel 0 is either sync or acceptor.
                                           alex_period_donor = (1500, 4096), #YKL, unit: channel, not sec.
                                           alex_period_acceptor = (0, 0),
                                           excitation_wavelengths = (488e-9, 0), #YKL, nanometers
                                           detection_wavelengths = (520e-9, 0), #YKL, nanometers
                                           )
            
            # Create a copy of metadata and remove the tags field
            m = meta.copy()
            tags = m.pop('tags')
            m
            #YKL, unit of time: sec
            #YKL, timestamps_unit: duration of each pulse. At 20 MHz, it is 50 ns. Each pulse is tagged with an interger.
            #YKL, nanotimes_unit: resolution of photoncounting channels (bins).
            
            # Print the file comment from the tags field
            print(tags['File_Comment']['data'])
            
            # Plot the alternation histogram from the data
            phc.plotter.alternation_hist(d)
            
            # Extract the detectors and their counts from the photon_data dictionary
            detectors = d['photon_data']['detectors']
            print("Detector    Counts")
            print("--------   --------")
            for det, count in zip(*np.unique(detectors, return_counts=True)):
                print("%8d   %8d" % (det, count))
            # YKL: Detertor 1 and 0 are true counts. All other detectors could be overflow.
            
            # Remove non-true detectors (overflow channels) from the photon_data dictionary
            valid = detectors != 15 # enter all nontrue detector numbers. 
            for field in ('detectors', 'timestamps', 'nanotimes'):
                d['photon_data'][field] = d['photon_data'][field][valid]
            
            # Check that the timestamps are monotonically increasing
            ts = d['photon_data']['timestamps']
            assert (np.diff(ts) >= 0).all()
            
            
            # 5. Conversion
            # Update the description, sample, and identity fields in the dictionary
            d['description'] = description
            
            d['sample'] = dict(
                sample_name=Experiment_ID,
                dye_names=dye_names,
                buffer_name=buffer_name,
                num_dyes = len(dye_names.split(',')))
            
            d['identity'] = dict(
                author=author,
                author_affiliation=author_affiliation)
            
            # Remove some empty groups that may cause errors on saving
            _ = meta.pop('dispcurve', None)
            _ = meta.pop('imghdr', None)
            
            # Save the data as an HDF5 file
            d['user'] = {'picoquant': meta}
            phc.hdf5.save_photon_hdf5(d, overwrite=True)


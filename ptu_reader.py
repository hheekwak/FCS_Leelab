# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 10:31:45 2023

@author: Hyunhee Kwak (hkwak9458)

Purpose: Convert PicoQuant photoncounting file to data dictionary and data files (_header.txt, _timestamps.csv)
       Work for PicoQuant .ptu file in both T2 and T3 mode

"""

import os
import csv
import phconvert as phc
import numpy as np
from itertools import zip_longest

#%%

# Assign variables for author, author affiliation, description, sample name, dye names, and buffer name.
author = 'Youngkwang Lee'
author_affiliation = 'SDSU'
description = 'time-resolved fluorescence spectroscopy'
Experiment_ID = 'Test data' # Enter experiment code, eg. FCS15
dye_names = 'ATTO488, '
buffer_name = 'PBS pH 7.4'

# Assign the directory path where the data files are stored
directory = "G:\My Drive\Hyunhee work_local\FCS_Leelab" # Enter the top directory

#%%
# The function "data_load" refers to phconvert.loader:
# https://github.com/Photon-HDF5/phconvert/blob/master/phconvert/loader.py

def data_load(filename):
    """Load PicoQuant file and build a data dictionary"""
    timestamps, detectors, nanotimes, metadata = phc.pqreader.load_ptu(filename)

    
    software = metadata.pop('software')
    software_version = metadata.pop('software_version')
    acquisition_duration = float(metadata.pop('acquisition_duration'))
    timestamps_unit = float(metadata.pop('timestamps_unit'))
    isT2 = metadata.pop('isT2')
    
    if isT2:
        pass
    else:
        laser_repetition_rate = float(metadata.pop('laser_repetition_rate'))
        tcspc_unit = float(metadata.pop('nanotimes_unit'))
        tcspc_num_bins = 4096
        tcspc_range = tcspc_num_bins * tcspc_unit

    # Creation time from the file header
    creation_time = metadata.pop('creation_time')

    provenance = dict(
        filename = filename,
        creation_time = creation_time,
        software = software,
        software_version = software_version,
    )
    if isT2:
        photon_data = dict(
            timestamps = timestamps,
            timestamps_specs = dict(timestamps_unit=timestamps_unit),
            detectors = detectors,
            nanotimes = nanotimes,
            
            measurement_specs = dict(
                measurement_type = 'T2',
            )
        )
    else: 
        photon_data = dict(
            timestamps = timestamps,
            timestamps_specs = dict(timestamps_unit=timestamps_unit),
            detectors = detectors,
            nanotimes = nanotimes,
            
            nanotimes_specs = dict(
                tcspc_unit = tcspc_unit,
                tcspc_range = tcspc_range,
                tcspc_num_bins = tcspc_num_bins),
    
            measurement_specs = dict(
                measurement_type = 'T3',
                laser_repetition_rate = laser_repetition_rate,
            )
        )

    setup = dict(
        num_pixels = 2,
        num_spots = 1,
        num_spectral_ch = 2,
        num_polarization_ch = 1,
        num_split_ch = 1,
        modulated_excitation = True,
        lifetime = True,
    )

    data = dict(
        _filename=filename,
        acquisition_duration = acquisition_duration,
        photon_data = photon_data,
        setup = setup,
        provenance = provenance)

    return data, metadata


#%%

def read_ptu(filepath):
    """Extract valid photon data and check its validation from a PTU file."""        
    d, metadata = data_load(filepath)
    
    # Create a copy of metadata and remove the tags field
    meta = metadata.copy()
    tags = meta.pop('tags')
    
    # Extract the detectors and their counts from the photon_data dictionary
    detectors = d['photon_data']['detectors']
    print("Detector    Counts")
    print("--------   --------")
    for det, count in zip(*np.unique(detectors, return_counts=True)):
        print("%8d   %8d" % (det, count))
    # YKL: Detertor 1 and 0 are true counts. All other detectors could be overflow.
    
    # Remove non-true detectors (overflow channels) from the photon_data dictionary
    valid = detectors != 15 # enter all nontrue detector numbers. 
    
    if meta['record_type'][-2:] == 'T2':
        for field in ('detectors', 'timestamps'):
                d['photon_data'][field] = d['photon_data'][field][valid]
    else:    
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
    
    return d, meta, tags

#%%# The function "ptu_header_file" refers to phconvert.pqreader._ptu_print_tags function :
# https://github.com/Photon-HDF5/phconvert/blob/master/phconvert/pqreader.py

def ptu_header_file(tags, filepath):
    """Write a header file of tags from a PTU file header."""
    headerfile = filepath + '_header.txt'
    
    with open(headerfile, "w+") as header:
        def _byte_to_str(x):
            if isinstance(x, bytes):
                # When loading from HDF5 string are binary
                x = x.decode()
            return x
    
        is_dol = True
        for tagname, tag in tags.items():
            if isinstance(tag, list):
                is_dol = False
                break
        if is_dol:
            tags = phc.pqreader._unconvert_multi_tags(tags)
        for n in tags:
            start = 'D'  # mark for duplicated tags
            tags_n = tags[n]
            if not isinstance(tags[n], list):
                tags_n = [tags_n]
                start = ' '
            for tag in tags_n:
                tag_type = _byte_to_str(tag["type"])
                
                if tag_type == 'tyBool8' or tag_type == 'tyInt8':
                    if tag['idx'] == -1:
                        value = str(tag["value"])
                        line = f'{start} {n:35s} {value}'
                    else:
                        value = "[" + str(tag['idx']) + "]  : " + str(tag["value"])
                        line = f'{start} {n:28s} {value}'    
                else:         
                    if tag_type == 'tyEmpty8':
                        value = "<empty Tag>"
                    
                    elif tag_type == 'tyFloat8':
                        value = f'{tag["value"]:20.4g}'
                    elif tag_type == 'tyAnsiString':
                        value = tag["data"]  
                    else:
                        value = f'{tag["value"]:<20}'
                    line = f'{start} {n:35s} {value}'
                    
                print(line, file= header, end = '\n')


#%%
def ptu_ts_file(d, filepath):
    """ Write a csv data file of channel - timestamps - nanotimes(T3 only) from a PTU file photon data."""
    
    ts_file = filepath + '_timestamps.csv'
    
    channel = d['photon_data']['detectors']
    ts_data = d['photon_data']['timestamps']
    ts_unit = d['photon_data']['timestamps_specs']['timestamps_unit']
    
    if d['photon_data']['measurement_specs']['measurement_type'] == 'T2':
        ts_d = [channel, ts_data, ts_data * ts_unit * 1E12]
        export_data_ts = zip_longest(*ts_d, fillvalue = '')
                    
        with open(ts_file, "w+", newline = '') as ts:
            wr = csv.writer(ts)
            wr.writerow(('channel', 'Time Tag', 'True Time(ps)'))
            wr.writerows(export_data_ts)
    else:
        nanotimes = d['photon_data']['nanotimes']
        ts_d = [channel, ts_data, np.around(ts_data * ts_unit * 1E9, decimals = 2), nanotimes]
        export_data_ts = zip_longest(*ts_d, fillvalue = '')
                    
        with open(ts_file, "w+", newline = '') as ts:
            wr = csv.writer(ts)
            wr.writerow(('channel', 'Time Tag', 'True Time(ns)', 'dTime'))
            wr.writerows(export_data_ts)


#%%

def main():
    if os.path.exists(directory):
        print('Directory found, it will proceed.')
    else:
        print('ATTENTION: Directory not found, please check the directory name.\n'
              '           (current value "%s")' % directory)
    
    # Walk through the directory tree and iterate over each file
    for root, dirs, files in os.walk(directory):
        for file in files:
            # Check if the file extension is '.ptu'
            if file.endswith('.ptu'):
                # Create the filepath by joining the root directory and file name
                filepath = os.path.join(root, file)
                d, meta, tags = read_ptu(filepath)
                name = file[:-4]
                newpath = filepath[:-4]
                
                # make a directory of the corresponing file  
                if not os.path.exists(newpath):
                    os.mkdir(newpath)
                
                filename = os.path.join(newpath, name)
                # header file write
                ptu_header_file(tags, filename)
                # timestamps file write
                ptu_ts_file(d, filename)
        
if __name__ == "__main__":
            main()
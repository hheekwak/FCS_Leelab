# FCS_Leelab project
Calculator for Fluorescence correlation curves from PicoQuant ptu file

## About FCS_Leelab project
This project's purpose is to reduce data process steps to obtain the fluorescence correlation spectroscopy (FCS) data within a valid (customized) range
### Histogram
- Draw TCSPC histrogram with a bin width set by user (defalut 200 ns)
### Correlation
- Time gating option (user defined)
- Resolution, number of points on each level of lag-time when to scale up, can be set by user (default 8)
- Autocorrelation for 1-channel data
- Autocorrelations, positive (a-b) and negative (b-a) cross-correlation values and graphs for 2-channel data
### Batch mode
- Continuously run the process for ptu files in designated directory
### Output files
Files are created in a folder named after .ptu file

- Metadata file: FILENAME_header.txt
  - Basic and environmental information 
  - File would not be overwritten when the file already exists
- Timestamp file: FILENAME_timestamps.csv
  - Record of photon arrival times in a cycle
  - File would not be overwritten once file created
- Correlation file: FILENAME_correlations.csv
  - File overwrites at every run 

## Getting Started

### Development environment
- Python 3
- PyQt5 (GUI)
- IDE: Spyder

### Installing
1. *Paul Müller (2012)* **Python multiple-tau algorithm (Version 0.3.3)** is used for binning and correlation calculation. 

    For more information about multiple-tau algorithm is available at <https://pypi.python.org/pypi/multipletau>

        pip install multipletau

2. "*phconvert* is a python 2 & 3 library that helps writing valid Photon-HDF5 files, a file format for time stamp-based single-molecule spectroscopy."

    For more information about phconvert <https://github.com/Photon-HDF5/phconvert>

        pip install phconvert

3. Replace pqreader.py file in directory “phconvert” for both T3 and T2 mode use

   - find the installed location of “phconvert" 
     - ﻿e.g. C:\ProgramData\anaconda3\Lib\site-packages\phconvert
   - change the name of original 'pqreader.py' file appropriately to keep the original file safely
     - e.g. pqreader_original.py   
   - copy the modified “pqreader.py” file in "pqreader_modified" folder at GitHub into phconvert folder at your local machine

## Running the test
- Run ptu_corr_multitau.py and follow the prompts (e.g. "Enter file name (*.ptu): " -> test.ptu)
- Run ptu_corr_GUI.py and browse test.ptu file and input user setting

## Authors
- Hyunhee Kwak hkwak9458(at)sdsu.edu
- Prof. Youngkwang Lee 

## Resources 
## Acknowledgements

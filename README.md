# FCS_Leelab project
Calculator for Fluorescence correlation curves from PicoQuant ptu file

## About FCS_Leelab project
This project's purpose is to reduce data process steps to obtain the fluorescence correlation spectroscopy (FCS) data within a valid (customized) range

<img width="828" alt="Screenshot 2024-07-23 at 10 46 14 AM" src="https://github.com/user-attachments/assets/78df455c-149f-489d-8a3b-3e226f47fadc">
<img width="1038" alt="Screenshot 2024-07-24 at 12 11 28 PM" src="https://github.com/user-attachments/assets/d4e498d6-0d7c-49f5-bd31-6efc1b4f5613">

### Histogram
- Draw TCSPC histrogram with a bin width set by user (default 200 ns)
  
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
     - e.g. Users/user_name/myenv/lib/python3.11/site-pachages/phconvert
   - change the name of original 'pqreader.py' file appropriately to keep the original file safely
     - e.g. pqreader_original.py   
   - copy the modified “pqreader.py” file in "pqreader_modified" folder at GitHub into phconvert folder at your local machine
  
4. Downgrade numpy to 1.x.x if you need (July 2024)
   - numpy version 2.0.1 is not working with the error message "ValueError: numpy.dtype size changed, may indicate binary incompatibility."

## Running the test
- Run ptu_corr_multitau.py and follow the prompts (e.g. "Enter file name (*.ptu): " -> test.ptu)
- Run ptu_corr_GUI.py and browse test.ptu file and input user setting

## Authors
- Hyunhee Kwak hkwak9458(at)sdsu.edu
- Prof. Youngkwang Lee 

## Resources 
## Acknowledgements

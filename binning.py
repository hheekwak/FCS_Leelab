# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 14:54:08 2023

@author: HYUNHEE KWAK
Purpose: bin the photon arraival counts with desired width
    
    
"""

import numpy as np


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
        
                
    

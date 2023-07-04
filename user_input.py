# -*- coding: utf-8 -*-
"""
Created on Sat Jul  1 03:55:26 2023

@author: hkwak9458
"""
import os

#%%
class UserInput:
    def __init__(self, filename, bin_width_ns = 200, resolution = 8, time_gate_op = False):
        if not (os.path.exists(filename)):
            raise ValueError('ATTENTION: File not found, please check the file name.\n'
                  '           (current value "%s")' % filename)
        
        else:
            self.filename = filename
            self.bin_width_ns = bin_width_ns
            self.resolution = resolution
            self.time_gate_op = time_gate_op
            print('File found, it will proceed.')


    #def set_time_gate(tg1_st = 0, tg1_en = 50, tg2_st = 0, tg2_en = 50): #TODO:ASK DR.LEE if max time (50 ns) changes 
        

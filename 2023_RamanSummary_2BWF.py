import os
import fnmatch
from pandas import *
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


Data = pd.DataFrame()
for file in os.listdir('.'):
    if fnmatch.fnmatch(file, '*_fitfile.txt'):
        fullFilename = file
        filename = file[:-12]
        print(filename)
        flen = len(filename)-9
        
        data = np.genfromtxt(fullFilename, delimiter = '\t', usecols = (1,), skip_header = (2))
        error= np.genfromtxt(fullFilename, delimiter = '\t', usecols = (2,), skip_header = (2))

        position = data[0]        
        exc_laser = data[1]
        NumPks = data[2]
        
        baseline_order = data[3]
        baseline_r2 = data[4]
        baseline_flatness = data[5]
        #model_r2 = data[2]
        qBWF = data[6]
        g_pos = data[7]
        g_width = data[8]
        g_intensity = data[9]
        g_int_error = error[9]
        
        d_pos = data[10]
        d_width = data[11]
        d_intensity = data[12]
        d_ints_error = error[12]
        
        d2_pos = data[13]
        d2_width = data[14]
        d2_intensity = data[15]
        d2_ints_error = error[15]
        
        d3_pos = data[16]
        d3_width = data[17]
        d3_intensity = data[18]
        d3_ints_error = error[18]
        
        d4_pos = data[19]
        d4_width = data[20]
        d4_intensity = data[21]
        d4_ints_error = error[21]
        
        La = data[22]
        La_error = error[22]
        ID_IG = data[24]
        ID_IG_error = error[24]
        #SNRD = data[11]
        #noise = data[12]
        
        ScanInfo = np.genfromtxt(fullFilename, dtype = str, delimiter = '\t', usecols = (1), skip_footer = (len(data)+1))
        
        Data = Data.append({'Position': position, 
            'Scan Info': ScanInfo, 
            'Exc Laser': exc_laser,
            'Num Peaks Fit': NumPks,
            
            'Baseline Flatness': baseline_flatness, 
            
            #'Model R2': model_r2, 
            'qBWF': qBWF, 
            'G Peak Position': g_pos, 
            'G Peak Width':g_width, 
            'G Intensity': g_intensity ,
            #'G Intensity error': g_error  ,
            'D Peak Position': d_pos,
            'D Peak Width': d_width, 
            'D Intensity': d_intensity,
            #'D Intensity error': d_ints_error,
            'D2 Peak Position': d2_pos,
            'D2 Peak Width': d2_width, 
            'D2 Intensity': d2_intensity,
            #'D2 Intensity error': d2_ints_error,
            'D3 Peak Position': d3_pos,
            'D3 Peak Width': d3_width, 
            'D3 Intensity': d3_intensity,
            #'D3 Intensity error': d3_ints_error,
            'D4 Peak Position': d4_pos,
            'D4 Peak Width': d4_width, 
            'D4 Intensity': d4_intensity,
            'Conj Length': La , 
            #'Conj Length error': La_error , 
            'ID IG Ratio': ID_IG ,
            #'ID IG Ratio error': ID_IG_error ,
            'Filename': filename ,
            #'Noise': noise ,
            #'SNR D band': SNRD ,  
            }, ignore_index=True)
Data.to_csv('RamanFit_summary.csv',index=False,header=True)  



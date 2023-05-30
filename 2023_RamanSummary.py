import os
import fnmatch
from pandas import *
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


#Data = pd.DataFrame()
tempData = []

# =============================================================================
# Position=[]
# ScanInfo=[]
# ExcLaser=[]
# NumPeaksFit=[]
# 
# BaselineFlatness=[]
# 
# qBWF=[]
# GPeakPosition=[]
# GPeakWidth=[]
# GIntensity=[]
# DPeakPosition=[]
# DPeakWidth=[]
# DIntensity=[]
# D2PeakPosition=[]
# D2PeakWidth=[]
# D2Intensity=[]
# D3PeakPosition=[]
# D3PeakWidth=[]
# D3Intensity=[]
# D4PeakPosition=[]
# D4PeakWidth=[]
# D4Intensity=[]
# ConjLength=[]
# IDIGRatio=[]
# Filename=[]
# =============================================================================


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
        num_pks = data[2]
        
        baseline_order = data[3]
        baseline_r2 = data[4]
        baseline_flatness = data[5]
        
        pkfit_r2 = data[6]
        pkfit_see = data[7]

        
        q_BWF = data[8]
        g_pos = data[9]
        g_width = data[10]
        g_intensity = data[11]
        g_int_error = error[11]
        
        d_pos = data[12]
        d_width = data[13]
        d_intensity = data[14]
        d_ints_error = error[14]
        
        d2_pos = data[15]
        d2_width = data[16]
        d2_intensity = data[17]
        d2_ints_error = error[17]
        
        d3_pos = data[18]
        d3_width = data[19]
        d3_intensity = data[20]
        d3_ints_error = error[20]
        
        d4_pos = data[21]
        d4_width = data[22]
        d4_intensity = data[23]
        d4_ints_error = error[23]
        
        Lalow_val = data[24]
        Lalow_error = error[24]
        ID_IG = data[26]
        ID_IG_error = error[26]
        
        scan_info = np.genfromtxt(fullFilename, dtype = str, delimiter = '\t', usecols = (1), skip_footer = (len(data)+1))

        #########################  
        # to solve concat/append error.  concat all the bits into individual lists, then after all loops run, concat into a df below

        
        tempData.append({'Position': position, 
            'Scan Info': scan_info, 
            'Exc Laser': exc_laser,
            'Num Peaks Fit': num_pks,
            
            'Baseline Order': baseline_order,
            'Baseline Flatness': baseline_flatness, 
            
            'PkFit R2': pkfit_r2, 
            'PkFit SEE': pkfit_see, 
            
            'qBWF': q_BWF, 
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
            'Conj Length (low)': Lalow_val , 
            #'Conj Length error': La_error , 
            'ID IG Ratio': ID_IG ,
            #'ID IG Ratio error': ID_IG_error ,
            'Filename': filename ,
            #'Noise': noise ,
            #'SNR D band': SNRD ,  
            })
        
# this is where concat the individual lists into dataframes will go--outside the loop  
Data = pd.DataFrame(tempData)          
Data.to_csv('RamanFit_summary.csv',index=False,header=True)  



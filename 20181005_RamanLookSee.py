# Imports .xls file from KOS Raman system, Excel export with both x and y values saved.
# makes spectrum plot for each 

import sys
import numpy as np
from pylab import *
import matplotlib.pyplot as plt
#from scipy import signal
from scipy.interpolate import interp1d
import os
import fnmatch
from pandas import *
import pandas as pd

# File Parameters

Loadfile =  'Kob_3dotJS.xls'
Info = ''


# Data fitting parameters

bkd_bounds = [950, 1100, 1750, 1900] #800-1000, 1800-2000#low wavelength limits (low, high) and high wavelength limits (low, high)
G_bounds = [1600, 100, 60, 30] # G Peak: Center wavelength, wavelength limits, HWHM guess, HWHM limits
D_bounds = [1350, 100, 100, 50]  # D Peak: Center wavelength, wavelength limits, HWHM guess, HWHM limits
Ext_Lambda = 785 #nm
step_size = 1
iterations = 100

lowG = G_bounds[0] - (G_bounds[1]/2)
highG = G_bounds[0] + (G_bounds[1]/2)
lowD = D_bounds[0] - (D_bounds[1]/2)
highD = D_bounds[0] + (D_bounds[1]/2)
    
lowD = 1100

all_data = pd.read_excel(Loadfile, skiprows = (9), header = None)

for n in range(0,len(all_data.columns)/2):
    Position = n
    wave = all_data[2*n].values
    signal = all_data[(2*n)+1].values
    
    new_info = pd.read_excel(Loadfile, dtype = str, usecols = (2*n,), skiprows = (8), skipfooter = (len(all_data)), header = None)
    if len(new_info.index) <> 0:
        Info = str(new_info[0].values)[2:-2]
        print Info
        
    #cut-down region wave and signal to region of interest
    wave_fit = wave[(bkd_bounds[0]-1 <= wave) & (wave <= bkd_bounds[3])]
    signal_fit = signal[(bkd_bounds[0]-1 <= wave) & (wave <= bkd_bounds[3])]
    #wave_base = wave_bkg[(bkd_bounds[0]-1 <= wave_bkg) & (wave_bkg <= bkd_bounds[3])]
    #signal_base = signal_bkg[(bkd_bounds[0]-1 <= wave_bkg) & (wave_bkg <= bkd_bounds[3])]
    
    wavesforfit =  np.arange(bkd_bounds[0], bkd_bounds[3], step_size) # Resizes Wavefile to step size
    tck = interp1d(wave_fit, signal_fit, bounds_error=False, fill_value=0) # Relates wavefile and signal files
    signal_fit = tck(wavesforfit) #Resizes Signalfile to stepsize
    #tck = interp1d(wave_base, signal_base, bounds_error=False, fill_value=0)
    #signal_base = tck(wavesforfit)
    wave_fit = wavesforfit
    #wave_base = wavesforfit
    
    
    # Plotting
    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    ax.plot(wave_fit, signal_fit,'.k')
    #ax.plot(wave_fit, baseline, '-r')
    ax.set_xlabel('Raman Shift (cm-1)')
    ax.set_ylabel('Raman Intensity')
    #ax.set_ylim(0, 1.5*max(signal_fit))
    #ax.legend(loc='best')
    
    plt.savefig(Loadfile + '_POS'+ str(Position).zfill(2)+ '_LookSee.jpg')
    plt.close()
    

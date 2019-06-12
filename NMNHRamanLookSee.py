# Imports .txt file from NMNH Raman spectrometer, single spectrum files!
# Loops to graph all spectra in folder

import sys
import numpy as np
from pylab import *
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from scipy import signal
from scipy.interpolate import interp1d
import os
import fnmatch
from pandas import *
import pandas as pd


for file in os.listdir('.'):
    if fnmatch.fnmatch(file, '*.txt'):
        fullFilename = file
        Loadfile = file[:-4]

        # Data fitting parameters
        bkd_bounds = [950, 1100, 1750, 1900] #800-1000, 1800-2000#low wavelength limits (low, high) and high wavelength limits (low, high)
        G_bounds = [1600, 100, 60, 30] # G Peak: Center wavelength, wavelength limits, HWHM guess, HWHM limits
        D_bounds = [1350, 100, 100, 50]  # D Peak: Center wavelength, wavelength limits, HWHM guess, HWHM limits
        step_size = 1
        iterations = 100
        
        lowG = G_bounds[0] - (G_bounds[1]/2)
        highG = G_bounds[0] + (G_bounds[1]/2)
        lowD = D_bounds[0] - (D_bounds[1]/2)
        highD = D_bounds[0] + (D_bounds[1]/2)
        
        lowD = 1100
        
        
        
        # Loads the Data files
        wave = np.loadtxt(Loadfile + '.txt', usecols = (0,), skiprows = (38))
        signal = np.loadtxt(Loadfile + '.txt', usecols = (1,), skiprows = (38))
        
        #cut-down region wave and signal to region of interest
        wave_fit = wave[(bkd_bounds[0]-1 <= wave) & (wave <= bkd_bounds[3])]
        signal_fit = signal[(bkd_bounds[0]-1 <= wave) & (wave <= bkd_bounds[3])]
        
        # Plotting
        fig = plt.figure(1)
        ax = fig.add_subplot(111)
        ax.plot(wave_fit , signal_fit,'.k')
        ax.set_xlabel('Raman Shift (cm-1)')
        ax.set_ylabel('Raman Intensity')
        #ax.set_ylim(0, 1.5*max(signal_fit))
        ax.legend(loc='best')
        
        plt.savefig(Loadfile + '_LookSee.jpg')
        plt.close()
        
        

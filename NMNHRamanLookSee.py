# Imports .txt file from NMNH Raman spectrometer, single spectrum files!
# Loops to graph all spectra in folder

import matplotlib.pyplot as plt
import os
import fnmatch
import pandas as pd


for file in os.listdir('.'):
    if fnmatch.fnmatch(file, '*.txt'):
        fullFilename = file
        Loadfile = file[:-4]        
        
        # Loads the Data files
        all_data = pd.read_csv(fullFilename, skiprows = (46), sep = None, header = None, engine='python', encoding='latin-1', decimal='.')
        
        wave = all_data[0]
        signal = all_data[1]
                
        # Plotting
        fig = plt.figure(1)
        ax = fig.add_subplot(111)
        ax.plot(wave, signal,'.k')
        ax.set_xlabel('Raman Shift / cm$^{-1}$')
        ax.set_ylabel('Raman Intensity')
        #ax.set_ylim(0, 1.5*max(signal_fit))
        #ax.legend(loc='best')
        
        plt.savefig(Loadfile + '_LookSee.jpg')
        plt.close()
        
        

# Imports .txt file from NMNH Raman spectrometer, single spectrum files!
# Loops to graph all spectra in folder

import matplotlib.pyplot as plt
import os
import fnmatch
import numpy as np
import pandas as pd

def ReadCollDetails(inputFilename):
    global CollDet, Ext_Lambda
    lineList = []
    i=0
    #txtFile = open(inputFilename, "r")
    with open(inputFilename, encoding='latin-1', errors='ignore') as txtFile:
        #print(inputFilename + 'got to subroutine')
        for line in txtFile:
            lineList.append(line.rstrip("\n"))
            if '#Acq. time' in lineList[i]: 
                acqtime = (lineList[i].split('=')[1][1:])
                #print(str(acqtime)+"acq")
            if '#Accumulations=' in lineList[i]: 
                NoScans = (lineList[i].split('=')[1][1:])
                #print(str(NoScans)+"sc")
            if '#Laser=' in lineList[i]:
                Laser = lineList[i].split('=')[1][1:]
                #print(Laser)
                Ext_Lambda = int(Laser[0:3])
                #print(Ext_Lambda)
            if '#ND' in lineList[i]:
                lp = lineList[i].split('=')[1][1:]
            if '#Date' in lineList[i]:
                DateLine = int(lineList[i].split()[1][6:])
                #print(DateLine)
            i = i+1
    CollDet = 'NMNH '+acqtime+'s '+NoScans+ 'sc ' + Laser + ' lp ' +lp
    txtFile.close()
    return CollDet,Ext_Lambda, DateLine

for file in os.listdir('.'):
    if fnmatch.fnmatch(file, '*.txt'):
        Loadfile = file
        print(file)
        filename = file[:-4]
        
        (CollDet,Ext_Lambda,FileDate) = ReadCollDetails(file)
        if FileDate < 2023:
            all_data = pd.read_csv(Loadfile, skiprows = (38), sep = None, header = None, engine='python', encoding='latin-1')
        else:
            all_data = pd.read_csv(Loadfile, skiprows = (46), sep = None, header = None, engine='python', encoding='latin-1')
        locations = all_data[0]
        all_data = all_data.transpose()
        
        xy_data = all_data[1:]
        wave = np.array(xy_data[0])

        for n in range(1,len(xy_data.columns)):
            signal = np.array(xy_data[n])
            # Plotting
            fig = plt.figure(1)
            ax = fig.add_subplot(111)
            ax.plot(wave, signal,'.k')
            ax.set_xlabel('Raman Shift / cm$^{-1}$')
            ax.set_ylabel('Raman Intensity')
            #ax.set_ylim(0, 1.5*max(signal_fit))
            #ax.legend(loc='best')
            plt.savefig(filename + '_POS'+ str(n).zfill(2)+'_LookSee.jpg')
            plt.close()
            
        

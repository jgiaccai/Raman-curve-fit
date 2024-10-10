# Imports .csv file from Wasatch Raman spectrometer
# Loops to graph all spectra in folder  

# Updated 10 Oct 2024 (Wasatch Enlighten v. 3.2.36)

import matplotlib.pyplot as plt
import os
import fnmatch
import numpy.polynomial.polynomial as poly
import linecache
import numpy as np

wavefile = "/Users/jennifergiaccai/OneDrive - Smithsonian Institution/PROJECTS/ChineseInkProject/WasatchWavenumRef-20210624-160600-788505-WP-00665.csv"

base_order = 9 #order of polynomial for bkg fitting, choose 1, 2, or 3
bkd_bounds = [80, 449, 450, 2000] #low wavelength limits (low, high) and high wavelength limits (low, high)


def ReadCollDetails(inputFilename):
    global CollDet, Ext_Lambda

    noscans = int(linecache.getline(inputFilename,9).split(',')[1])
    #print(noscans)
    scantime_ms = int(linecache.getline(inputFilename,28).split(',')[1])  #milliseconds
    Ext_Lambda = float(linecache.getline(inputFilename,40).split(',')[1])  #nm
    laserpower = float(linecache.getline(inputFilename,42).split(',')[1])
    
    scantime = scantime_ms/1000
    
    CollDet = 'Wasach '+str(Ext_Lambda) +'nm '+ str(noscans)+'sc '+ str(scantime) + 's ' + str(laserpower) + 'lp'

    return CollDet,Ext_Lambda



for file in os.listdir('.'):
    if fnmatch.fnmatch(file, '*.csv'):
        Loadfile = file
        print(Loadfile)
        filename = file[:-20]
        
        (CollDet,Ext_Lambda) = ReadCollDetails(file)
        
        wavenum_yes_no = str(linecache.getline(Loadfile,46).split(',')[0])

        if wavenum_yes_no == 'Wavenumber':
            wavenum = np.loadtxt(Loadfile, usecols = 0, skiprows = (47),encoding='latin1', delimiter = ',')
            signal = np.loadtxt(Loadfile, usecols = 1, skiprows = (47),encoding='latin1', delimiter = ',')
        elif wavenum_yes_no == 'Processed\n':
            wavenum = np.loadtxt(wavefile, usecols = 0, skiprows = (47),encoding='latin1', delimiter = ',')
            signal = np.loadtxt(Loadfile, usecols = 0, skiprows = (47),encoding='latin1', delimiter = ',')
        else:
            print('cannot find data')
            continue
                        
        # Plotting
        fig = plt.figure(1)
        ax = fig.add_subplot(111)
        ax.plot(wavenum, signal,'.k')
        ax.set_xlabel('Raman Shift / cm$^{-1}$')
        ax.set_ylabel('Raman Intensity')
        ax.tick_params(direction="in")
        #ax.set_ylim(0, 1.5*max(signal_fit))
        #ax.legend(loc='best')
        
        plt.savefig(filename + '_LookSee.jpg')
        plt.close()
        
        #Fit Baseline and cut-down region wavenum and signal to region of interest
        x_fit = wavenum[(bkd_bounds[0] <= wavenum) & (wavenum <= bkd_bounds[3])] # X data cut down to bkg boundaries
        signal_fit = signal[(bkd_bounds[0] <= wavenum) & (wavenum <= bkd_bounds[3])] # y data cut down to bkg boundaries
        
        lower_bkg = (bkd_bounds[0] <= x_fit) & (x_fit <= bkd_bounds[1]) # True/False selection for low part of bkg
        higher_bkg = (bkd_bounds[2] <= x_fit) & (x_fit <= bkd_bounds[3]) # True/False selection for high part of bkg
        
        bkg_x = x_fit[lower_bkg | higher_bkg] # pipe is essentially acting as numpy.or to select only high and low bkg areas
        bkg_signal = signal_fit[lower_bkg | higher_bkg]
        
        ##  baseline fit
        BasePara, BaseRegrStat = poly.polyfit(bkg_x, bkg_signal, base_order, full = True)
        baseline = poly.polyval(x_fit, BasePara)
        baseline_bkgarea = poly.polyval(bkg_x, BasePara)  
        BaseResiduals = bkg_signal-baseline_bkgarea
        
        fig = plt.figure(1)
        gs1 = fig.add_gridspec(nrows=2, ncols=1,
                  hspace=0.1, wspace=0, height_ratios=[1, 5])
        ax = fig.add_subplot(gs1[1])
        ax.plot(wavenum , signal,'-k')
        ax.plot(x_fit, baseline, '-r')
        ax.plot(bkg_x , bkg_signal,'.b')
        ax.set_xlabel('Raman Shift / cm$^{-1}$')
        ax.set_ylabel('Raman Intensity')
        ax.set_ylim(0, 1.5*max(signal_fit))
        ax.tick_params(direction="in")
        plt.autoscale(enable=True, axis='x', tight=True)
        plt.autoscale(enable=True, axis='y', tight=True)
        
        ax11 = fig.add_subplot(gs1[0])
        ax11.plot(bkg_x,BaseResiduals, '.b')
        ax11.axhline(0, linestyle='--', color = 'gray', linewidth = 1)
        ax11.set_ylabel('Residuals')
        ax11.tick_params(direction='in',labelbottom=True,labelleft=True)
        plt.autoscale(enable=True, axis='x', tight=True) 
        plt.ylim(min(BaseResiduals)*1.15,max(BaseResiduals)*1.15)
        
        gs1.update(left=0.13,right=0.96,top=0.95,bottom=0.12) #as percentages of total figure with 1,1 in upper right
        fig.set_size_inches(6, 5) #width, height
        fname = filename + '_base.jpg'
        plt.savefig(fname, dpi=300)
        plt.close()
        
        signal_base = signal - baseline #corrects for baseline 
        
        fig = plt.figure(2)
        ax = fig.add_subplot(111)
        ax.plot(wavenum, signal_base,'-k')
        ax.set_xlabel('Raman Shift / cm$^{-1}$')
        ax.set_ylabel('Raman Intensity, corrected')
        ax.tick_params(direction="in")
        #ax.set_ylim(0, 1.5*max(signal_fit))
        #ax.legend(loc='best')
        
        plt.savefig(filename + '_BaseCorr.jpg')
        plt.close()
        


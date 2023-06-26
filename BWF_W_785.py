## Uses .csv file from Wasatch instrument
## 20191002 expanding bounds based on scatter plots
## 20201207 getting scan details from metadata
## 20230408 altered to allow fitting of any number of Lorentzian peaks, defined lines 29-34
## 20230522 altered to fit G peak as a BWF with q of -10
## 20230524 choose baseline fitting order at top 

import sys
import numpy as np
from scipy.optimize import curve_fit
from scipy.optimize import minimize
from scipy.optimize import least_squares
import numpy.polynomial.polynomial as poly
import matplotlib.pyplot as plt
#from matplotlib.gridspec import GridSpec
from scipy import signal
import os
import fnmatch
#import random
#import math
#from time import sleep
#import pandas as pd
import linecache
import uncertainties

# File Parameters
# =============================================================================
# # Data fitting parameters
# 
FitGOn = 1 # 1 is yes, 0 is no
FitDOn = 1
FitD2On = 0
FitD3On = 0 
FitD4On = 0
FitU1On = 0 #unidentified peak but need to include in envelope for 405 etc

fitVersion = 2.0 #changing if there is a change to base fitting subtr or peak fitting or stats calc.  Not for making figures or summarizing data.

base_order = 3 #order of polynomial for bkg fitting, choose 1, 2, or 3
bkd_bounds = [520, 950, 1750, 2000] #low wavelength limits (low, high) and high wavelength limits (low, high)

G_bounds = [1590, 50, 50, 40] # Center wavelength, wavelength limits, HWHM guess, HWHM limits (currently unused)
D_bounds = [1350, 60, 100, 60]
D2_bounds = [1620, 10, 20, 10]
D3_bounds = [1500, 10, 45, 40]
D4_bounds = [1225, 10, 60, 40]
U1_bounds = [1725, 20, 10, 8]  #no physical basis, trying because weird peak in some 405 data
IIM = 0.8 #Initial intensity multiplier for G and D peaks 
qBWF = -10
# =============================================================================

#at some point files got switched so not saving the wavenumber column...grr.
#adding a file that has the wavenumbers in a column

wavefile = "/Users/jennifergiaccai/OneDrive - Smithsonian Institution/PROJECTS/ChineseInkProject/WasatchWavenumRef-20210624-160600-788505-WP-00665.csv"

Ext_Lambda = 000 #nm
#iterations = 50

TotalNumPeaks = FitGOn + FitDOn + FitD2On + FitD3On + FitD4On + FitU1On
NumPeaks = FitGOn + FitDOn + FitD2On + FitD3On + FitD4On
NumPksApp = str(NumPeaks)+'BWF'
NumParams    = 3*6     #{Number of parameters to fit}
FitParam =np.zeros(NumParams) 
lobounds = np.zeros(18)
hibounds = np.zeros(18)

# lowG = G_bounds[0] - (G_bounds[1]/2)
# highG = G_bounds[0] + (G_bounds[1]/2)
# lowD = D_bounds[0] - (D_bounds[1]/2)
# highD = D_bounds[0] + (D_bounds[1]/2)
    
# lowD = 1100
CollDet = 'NMNH unk det'  #If need to restart, replace with collection details

# only need to set peak location and width bounds once, at the beginning of the program  
# bounds for optimize.curve_fit: an array of lows, an array of highs

lobounds[0] = (G_bounds[0]-G_bounds[1])
lobounds[1] = (D_bounds[0]-D_bounds[1])
lobounds[2] = (D2_bounds[0]-D2_bounds[1])
lobounds[3] = (D3_bounds[0]-D3_bounds[1])
lobounds[4] = (D4_bounds[0]-D4_bounds[1])
lobounds[5] = (U1_bounds[0]-U1_bounds[1])

lobounds[6] = (G_bounds[2]-G_bounds[3])
lobounds[7] = (D_bounds[2]-D_bounds[3])
lobounds[8] = (D2_bounds[2]-D2_bounds[3])
lobounds[9] = (D3_bounds[2]-D3_bounds[3])
lobounds[10] = (D4_bounds[2]-D4_bounds[3])
lobounds[11] = (U1_bounds[2]-U1_bounds[3])

hibounds[0] = (G_bounds[0]+G_bounds[1])
hibounds[1] = (D_bounds[0]+D_bounds[1])
hibounds[2] = (D2_bounds[0]+D2_bounds[1])
hibounds[3] = (D3_bounds[0]+D3_bounds[1])
hibounds[4] = (D4_bounds[0]+D4_bounds[1])
hibounds[5] = (U1_bounds[0]+U1_bounds[1])

hibounds[6] = (G_bounds[2]+G_bounds[3])
hibounds[7] = (D_bounds[2]+D_bounds[3])
hibounds[8] = (D2_bounds[2]+D2_bounds[3])
hibounds[9] = (D3_bounds[2]+D3_bounds[3])
hibounds[10] = (D4_bounds[2]+D4_bounds[3])
hibounds[11] = (U1_bounds[2]+U1_bounds[3])

hibounds[12] = 500 #made up intensities just to fill the array out will customize per spectrum
hibounds[13] = 500
hibounds[14] = 500
hibounds[15] = 100
hibounds[16] = 100
hibounds[17] = 500

bounds = lobounds,hibounds
  
#print('Bounds are ',bounds)

def ReadCollDetails(inputFilename):
    global CollDet, Ext_Lambda

    noscans = int(linecache.getline(inputFilename,7).split(',')[1])
    scantime_ms = int(linecache.getline(inputFilename,15).split(',')[1])  #milliseconds
    Ext_Lambda = float(linecache.getline(inputFilename,25).split(',')[1])  #nm
    laserpower = float(linecache.getline(inputFilename,27).split(',')[1])
    
    scantime = scantime_ms/1000
    
    CollDet = 'Wasach '+str(Ext_Lambda) +'nm '+ str(noscans)+'sc '+ str(scantime) + 's ' + str(laserpower) + 'lp'

    return CollDet,Ext_Lambda

#physics suggests Raman of solid should be Gaussian--why not here? 
#(gases Lorentzian, liquids Gaussian-Laurentzian or Voigt)

def lorentz(xc, w, I):  
    global x_fit
    s = ((x_fit - xc)/w)
    return (I)*(1/(1+s**2)) #Wikipedia definition.  gives correct value for intensity.  using this definition fit_results you don't need to account for the peak width in the intensity.

def BWF(xc,w,I):
    global qBWF,x_fit
    s = ((x_fit - xc)/w)
    return (I)*((1+s/qBWF)**2/(1+s**2))

def EnterData():
    
    global FitParam, G_ints, D_ints, U1_ints, NumParams, G_bounds, D_bounds, D2_bounds, D3_bounds, D4_bounds, U1_bounds
    global bounds, lobounds, hibounds

    FitParam[0] = G_bounds[0] # G peak position
    FitParam[1]  = D_bounds[0] # D peak position
    FitParam[2]  = D2_bounds[0] # D2 peak position
    FitParam[3]  = D3_bounds[0] # D3 peak position
    FitParam[4]  = D4_bounds[0] # D4 peak position
    FitParam[5]  = U1_bounds[0] # U1 peak position

    FitParam[6]  = G_bounds[2]  # G peak width
    FitParam[7]  = D_bounds[2]  # D peak width
    FitParam[8]  = D2_bounds[2]  # D2 peak width
    FitParam[9]  = D3_bounds[2]  # D3 peak width
    FitParam[10]  = D4_bounds[2]  # D4 peak width
    FitParam[11]  = U1_bounds[2] # U1 peak width

    '''need to work G and D starting peak intensities from initial values'''
    FitParam[12]  = IIM*G_ints  # G peak intensity  
    FitParam[13]  = IIM*D_ints  # D peak intensity
    FitParam[14] = (0.4*FitParam[12])
    FitParam[15] = (0.4*FitParam[13])
    FitParam[16] = (0.4*FitParam[13])
    FitParam[17]  = (0.25*FitParam[12]) 

    Gfit = FitGOn*BWF(FitParam[0],FitParam[6],FitParam[12])
    Dfit = FitDOn*lorentz(FitParam[1],FitParam[7],FitParam[13])
    D2fit = FitD2On*lorentz(FitParam[2],FitParam[8],FitParam[14])
    D3fit = FitD3On*lorentz(FitParam[3],FitParam[9],FitParam[15])
    D4fit = FitD4On*lorentz(FitParam[4],FitParam[10],FitParam[16])
    U1fit = FitU1On*lorentz(FitParam[5],FitParam[11],FitParam[17])

    ModelFit = Gfit + Dfit + D2fit + D3fit + D4fit + U1fit    
    
    # # figure with initial fit of various peaks
    
    # fig = plt.figure(5)
    # ax50 = fig.add_subplot(111)
    # ax50.plot(x_fit, signal_fit,'.k', label = 'Experimental')
    # ax50.plot(x_fit, Gfit,'-g', label = 'G Peak Fit')
    # ax50.plot(x_fit, Dfit,'-b', label = 'D Peak Fit')
    # ax50.plot(x_fit, D2fit,'-y', label = 'D2 Peak Fit')
    # ax50.plot(x_fit, D3fit,'-c', label = 'D3 Peak Fit')
    # ax50.plot(x_fit, D4fit,'-m', label = 'D4 Peak Fit')
    # ax50.plot(x_fit, U1fit, '.g', label = 'U1 Peak Fit')
    # ax50.plot(x_fit, ModelFit,'-r', label = 'Summed Peak Fit')
    # ax50.set_xlabel(r'Raman Shift / cm$^-$$^1$')
    # plt.autoscale(enable=True, axis='x', tight=True)
    # plt.autoscale(enable=True, axis='y', tight=True)
    # ax50.set_ylabel('Raman Intensity')
    # ax50.set_ylim(min(signal_fit), max(signal_fit)*1.2)
    # #plt.text(1075, 14100, 'ink', fontsize=20)
    # plt.tick_params(axis='both', which='major')
    # #plt.tight_layout()
    # plt.savefig(SaveName + '_initialfit.jpg')
    # #plt.show()
    # plt.close()
    
    # need to add intensity bounds for each spectrum      
    # bounds for optimize.curve_fit: an array? of lows, an array? of highs, just need to update highs

    lobounds[12] = 0.5*G_ints  # we might want to let G go to zero depending on D2
    lobounds[13] = 0.5*D_ints
    hibounds[12] = 1.2*G_ints
    hibounds[13] = 1.2*D_ints

    hibounds[14] = 0.5*G_ints
    hibounds[15] = 0.5*D_ints
    hibounds[16] = 0.5*D_ints
    hibounds[17] = 0.5*G_ints
    
    bounds = (lobounds,hibounds)
    
def FitFunc(x_fit, *EvalSimp):   
    #global NumParams, G_bounds, D_bounds, D2_bounds, D3_bounds, D4_bounds, U1_bounds
    #global signal_fit, Residuals
    
    '''
    Need to 
    (1) evaluate Lorentzian for each peak
    (2) Add all peak fits together for total peak fit
    (3) Subtract total peak fit from real data for initial residuals
    '''
    
    Gfit = FitGOn*BWF(EvalSimp[0],EvalSimp[6],EvalSimp[12])
    Dfit = FitDOn*lorentz(EvalSimp[1],EvalSimp[7],EvalSimp[13])
    D2fit = FitD2On*lorentz(EvalSimp[2],EvalSimp[8],EvalSimp[14])
    D3fit = FitD3On*lorentz(EvalSimp[3],EvalSimp[9],EvalSimp[15])
    D4fit = FitD4On*lorentz(EvalSimp[4],EvalSimp[10],EvalSimp[16])
    U1fit = FitU1On*lorentz(EvalSimp[5], EvalSimp[11], EvalSimp[17])

    
    FitY = Gfit + Dfit + D2fit + D3fit + D4fit + U1fit
    
    #Residuals = (signal_fit - EvalFit)
    #ErrorSum = np.sum((Residuals)**2)  #fitting routine minimizes sum of square of residuals
    
    return(FitY)




for file in os.listdir('.'):

    if fnmatch.fnmatch(file, '*.csv'):
        Loadfile = file
        print(Loadfile)
        filename = file[:-20]
        
        (CollDet,Ext_Lambda) = ReadCollDetails(file)
        
        wavenum_yes_no = str(linecache.getline(Loadfile,31).split(',')[0])

        if wavenum_yes_no == 'Wavenumber':
            wavenum = np.loadtxt(Loadfile, usecols = 0, skiprows = (32),encoding='latin1', delimiter = ',')
            signal = np.loadtxt(Loadfile, usecols = 1, skiprows = (32),encoding='latin1', delimiter = ',')
        elif wavenum_yes_no == 'Processed\n':
            wavenum = np.loadtxt(wavefile, usecols = 0, skiprows = (32),encoding='latin1', delimiter = ',')
            signal = np.loadtxt(Loadfile, usecols = 0, skiprows = (32),encoding='latin1', delimiter = ',')
        else:
            print('cannot find data')
            continue
        
        SaveName = 'results/' + filename + '_'+ NumPksApp
        #print(SaveName)
        
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
        BaseTSS = ((bkg_signal - np.mean(bkg_signal))**2).sum()
        BaseR2= 1-(BaseRegrStat[0]/BaseTSS) #should work regardless of order of fit
        # calculate rise/run for fitted baseline immediately before and after the peak region.
        BaseRise = poly.polyval(bkd_bounds[2], BasePara) - poly.polyval(bkd_bounds[1], BasePara)
        BaseRun = bkd_bounds[2]-bkd_bounds[1]
        PseudoSlope = BaseRise/BaseRun
        
        if base_order ==1:
            BaseSlope = BasePara[1] #only works for first order baseline fits
        else:
            BaseSlope = PseudoSlope
        
        BaseResiduals = bkg_signal-baseline_bkgarea
        
        # figure with baseline fit and residuals
        
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
        plt.autoscale(enable=True, axis='x', tight=True)
        plt.autoscale(enable=True, axis='y', tight=True)
        
        ax11 = fig.add_subplot(gs1[0])
        ax11.plot(bkg_x,BaseResiduals, '.b')
        ax11.axhline(0, linestyle='--', color = 'gray', linewidth = 1)
        ax11.set_ylabel('Residuals')
        plt.setp(ax11, xticks=[1000,1200,1400,1600,1800,2000])
        ax11.tick_params(direction='in',labelbottom=True,labelleft=True)
        plt.autoscale(enable=True, axis='x', tight=True) 
        plt.ylim(min(BaseResiduals)*1.15,max(BaseResiduals)*1.15)
        
        #plt.show()
        gs1.update(left=0.13,right=0.96,top=0.95,bottom=0.12) #as percentages of total figure with 1,1 in upper right
        fig.set_size_inches(6, 5) #width, height
        fname = str(SaveName) + '_base.jpg'
        plt.savefig(fname, dpi=300)
        plt.close()
        
        # Baseline Correction
        signal_fit = signal_fit - baseline #corrects for baseline 
        
        #Finds an initial G intensity and D intensity to give a decent start to the peak fitting
        G_ints = max(signal_fit[(G_bounds[0]-G_bounds[1] <= x_fit) & (x_fit <= G_bounds[0]+G_bounds[1])])
        D_ints = max(signal_fit[(D_bounds[0]-D_bounds[1] <= x_fit) & (x_fit <= D_bounds[0]+D_bounds[1])]) 
        U1_ints = max(signal_fit[(U1_bounds[0]-U1_bounds[1] <= x_fit) & (x_fit <= U1_bounds[0]+U1_bounds[1])])
        #print 'initial D/G intensity ratio: ',D_ints/G_ints 
        
        
        EnterData()
        #print 'FitParam ', FitParam
        minres = curve_fit(FitFunc,x_fit,signal_fit,p0=FitParam,method = 'trf', bounds=bounds, full_output=True) 
        if minres[4] > 1 and minres[4] < 4:
            fit_results=minres[0]
            covMatrix = minres[1] 
            dFit = np.sqrt(np.diag(covMatrix))
            #print(fit_results)
            #iterations = minres.nfev # I think? not sure 
        else:
            print(minres[3])
            fit_results = np.zeros(NumParams)
            continue
        Gfit = FitGOn*BWF(fit_results[0],fit_results[6],fit_results[12])
        Dfit = FitDOn*lorentz(fit_results[1],fit_results[7],fit_results[13])
        D2fit = FitD2On*lorentz(fit_results[2],fit_results[8],fit_results[14])
        D3fit = FitD3On*lorentz(fit_results[3],fit_results[9],fit_results[15])
        D4fit = FitD4On*lorentz(fit_results[4],fit_results[10],fit_results[16])
        U1fit = FitU1On*lorentz(fit_results[5],fit_results[11],fit_results[17])

        ModelFit = Gfit + Dfit +D2fit + D3fit + D4fit + U1fit
        
        G_ints = fit_results[12]  #intensities from the fit, not the initial values
        D_ints = fit_results[13]
        D2_ints = FitD2On*fit_results[14]
        D3_ints = FitD3On*fit_results[15]
        D4_ints = FitD4On*fit_results[16]
        U1_ints = fit_results[17]
        TotalIntensity = G_ints + D_ints + D2_ints + D3_ints + D4_ints
       
        Residuals = (signal_fit - ModelFit)
        ss_res = np.sum((Residuals) ** 2)
        ss_tot = np.sum((signal_fit - np.mean(signal_fit)) ** 2)
        R2_fit = 1-ss_res/ss_tot # not meaningful because can't use R2 on Gaussian or Lorentzian fits, so need to use standard error regression
        SEE_fit = np.sqrt(ss_res/(len(signal_fit)-NumPeaks*3))
        
        # Figure 4 Plot of individual peak fits and total peak fit with experimental
        fig = plt.figure(4)
        gs4 = fig.add_gridspec(nrows=2, ncols=1,
                  hspace=0, wspace=0, height_ratios=[1, 5])
        ax40 = fig.add_subplot(gs4[1])
        ax40.plot(x_fit, signal_fit,'.k', label = 'Experimental')
        ax40.plot(x_fit, Gfit,'-g', label = 'G Peak Fit')
        ax40.plot(x_fit, Dfit,'-b', label = 'D Peak Fit')
        if FitD2On == 1:
            ax40.plot(x_fit , D2fit,'-y', label = 'D2 Peak Fit')
        if FitD3On == 1:
            ax40.plot(x_fit, D3fit,'-c', label = 'D3 Peak Fit')
        if FitD4On == 1:
            ax40.plot(x_fit , D4fit,'-m', label = 'D4 Peak Fit')
        if FitU1On == 1:
            ax40.plot(x_fit, U1fit, '-y', label = 'U1 peak fit')
        ax40.plot(x_fit, ModelFit,'-r', label = 'Summed Peak Fit')
        ax40.set_xlabel(r'Raman Shift / cm$^{-1}$', fontsize=16)
        plt.autoscale(enable=True, axis='x', tight=True)
        plt.autoscale(enable=True, axis='y')
        plt.setp(ax40, xticks=[1000,1200,1400,1600,1800,2000])
        ax40.set_ylabel('Raman Intensity', fontsize=16)
        plt.tick_params(axis='both', which='major', labelsize=14)
        
        ax41 = fig.add_subplot(gs4[0])
        ax41.plot(x_fit,Residuals, '.b')
        ax41.axhline(0, linestyle='--', color = 'gray', linewidth = 1)
        ax41.set_ylabel('Residuals')
        plt.setp(ax41, xticks=[800,1000,1200,1400,1600,1800])
        ax41.tick_params(direction='in',labelbottom=False,labelleft=True)
        plt.autoscale(enable=True, axis='x', tight=True) 
        
        plt.ylim(min(Residuals)*1.15,max(Residuals)*1.15)
        
        gs4.update(left=0.13,right=0.94,top=0.95,bottom=0.15) #as percentages of total figure with 1,1 in upper right
        fig.set_size_inches(6, 5) #width, height
        plt.savefig(SaveName + '_fit.jpg', dpi=300)
        plt.close()
        
        # ID/IG Ratio
        Exp_ratio = D_ints/G_ints
        
        # Calculation of uncertainties
        
        Exp_ratio_stdev = Exp_ratio*((dFit[12]/G_ints)**2 + (dFit[13]/D_ints)**2 - 
                                   (2*(covMatrix[13,12])/G_ints/D_ints))**0.5
        TotIntstdev = (dFit[12]**2 + dFit[13]**2 + dFit[14]**2 + dFit[15]**2 + dFit[16]**2 + 
                      2*((covMatrix[13,12]) + (covMatrix[13,14]) + 
                          (covMatrix[13,15]) + (covMatrix[13,16]) + 
                          (covMatrix[14,12]) + (covMatrix[14,15]) +
                          (covMatrix[14,16]) + (covMatrix[15,12]) + 
                          (covMatrix[15,16]) + (covMatrix[16,12])))**0.5
        #IDITstdev =
        #ID2ITstdev = 
        #ID3ITstdev = 
        #ID4ITstdev = 
        #IGITstdev = 
        if FitD2On == 1:
            ID2IGstdev = (D2_ints/G_ints)*((dFit[14]/D2_ints)**2 + (dFit[12]/G_ints)**2 - (2*(covMatrix[14,12])**0.5/G_ints/D2_ints))**0.5
        if FitD3On == 1:
            ID3IDstdev = (D3_ints/D_ints)*((dFit[15]/D3_ints)**2 + (dFit[13]/D_ints)**2 - (2*(covMatrix[13,15])**0.5/D3_ints/D_ints))**0.5
        if FitD4On == 1:
            ID4IDstdev = (D4_ints/D_ints)*((dFit[16]/D4_ints)**2 + (dFit[13]/D_ints)**2 - (2*(covMatrix[13,16])**0.5/D4_ints/D_ints))**0.5
        
        ra = 3.1 # +/- 0.03nm
        rs = 1.0 #
        CA = 160*((1240*Ext_Lambda**-1)**-4) # 160 +/- 48
        CA_min = (160-48)*((1240*Ext_Lambda**-1)**-4)
        CA_max = (160+48)*((1240*Ext_Lambda**-1)**-4)
        
        La_Range = np.arange(0.1, 100, step = 0.001)
        Model_Ratio = [0]*len(La_Range)
        Min_Ratio = [0]*len(La_Range)
        Max_Ratio = [0]*len(La_Range)
        for i in range(0, len(La_Range)):
            La = La_Range[i]
            Model_Ratio[i] =CA*(ra**2 - rs**2)/(ra**2 - 2*(rs**2))*(np.exp((-np.pi*rs**2)/(La**2))-np.exp((-np.pi*(ra**2 - rs**2))/(La**2))) 
            Min_Ratio[i] =CA_min*(ra**2 - rs**2)/(ra**2 - 2*(rs**2))*(np.exp((-np.pi*rs**2)/(La**2))-np.exp((-np.pi*(ra**2 - rs**2))/(La**2))) 
            Max_Ratio[i] =CA_max*(ra**2 - rs**2)/(ra**2 - 2*(rs**2))*(np.exp((-np.pi*rs**2)/(La**2))-np.exp((-np.pi*(ra**2 - rs**2))/(La**2)))
        
        low = np.where(La_Range <=3)[0]
        La_low = [0]*len(low)
        Model_Ratio_low = [0]*len(low)
        Min_Ratio_low = [0]*len(low)
        Max_Ratio_low = [0]*len(low)
        for i in range(0, len(low)):
            La_low[i] = La_Range[low[i]]
            Model_Ratio_low[i] = Model_Ratio[low[i]]
            Min_Ratio_low[i] = Min_Ratio[low[i]]
            Max_Ratio_low[i] = Max_Ratio[low[i]]  
        low_La = La_low[np.where(abs(Model_Ratio_low - Exp_ratio) == min(abs(Model_Ratio_low - Exp_ratio)))[0][0]]
        low_La_min_Ca = La_low[np.where(abs(Min_Ratio_low - Exp_ratio) == min(abs(Min_Ratio_low - Exp_ratio)))[0][0]]
        low_La_max_Ca = La_low[np.where(abs(Max_Ratio_low - Exp_ratio) == min(abs(Max_Ratio_low - Exp_ratio)))[0][0]]
        low_La_stdev_Ca = abs(low_La_min_Ca - low_La_max_Ca)
        
        #low_La_min_ratio = La_low[np.where(abs(Model_Ratio_low - (Exp_ratio-Exp_ratio_stdev) ) == min(abs(Model_Ratio_low - (Exp_ratio-Exp_ratio_stdev))))[0][0]]
        #low_La_max_ratio = La_low[np.where(abs(Model_Ratio_low - (Exp_ratio+Exp_ratio_stdev)) == min(abs(Model_Ratio_low - (Exp_ratio+Exp_ratio_stdev))))[0][0]]
        #low_La_stdev_ratio = abs(low_La_min_ratio - low_La_max_ratio)
        #low_La_stdev = low_La*((low_La_stdev_ratio/low_La)**2 + (low_La_stdev_Ca/low_La)**2)**0.5
        #print 'Conjugation Length (low): ', '%.3f' %low_La,  '+/- ', '%.3f' %low_La_stdev, 'nm'
        
        high = np.where(La_Range >=3)[0]
        La_high = [0]*len(high)
        Model_Ratio_high = [0]*len(high)
        Min_Ratio_high = [0]*len(high)
        Max_Ratio_high = [0]*len(high)
        for i in range(0, len(high)):
            La_high[i] = La_Range[high[i]]
            Model_Ratio_high[i] = Model_Ratio[high[i]]
            Min_Ratio_high[i] = Min_Ratio[high[i]]
            Max_Ratio_high[i] = Max_Ratio[high[i]]  
        high_La = La_high[np.where(abs(Model_Ratio_high - Exp_ratio) == min(abs(Model_Ratio_high - Exp_ratio)))[0][0]]
        high_La_min_Ca = La_high[np.where(abs(Min_Ratio_high - Exp_ratio) == min(abs(Min_Ratio_high - Exp_ratio)))[0][0]]
        high_La_max_Ca = La_high[np.where(abs(Max_Ratio_high - Exp_ratio) == min(abs(Max_Ratio_high - Exp_ratio)))[0][0]]
        high_La_stdev_Ca = abs(high_La_min_Ca - high_La_max_Ca)
        
        #high_La_min_ratio = La_high[np.where(abs(Model_Ratio_high - (Exp_ratio-Exp_ratio_stdev) ) == min(abs(Model_Ratio_high - (Exp_ratio-Exp_ratio_stdev))))[0][0]]
        #high_La_max_ratio = La_high[np.where(abs(Model_Ratio_high - (Exp_ratio+Exp_ratio_stdev)) == min(abs(Model_Ratio_high - (Exp_ratio+Exp_ratio_stdev))))[0][0]]
        #high_La_stdev_ratio = abs(high_La_min_ratio - high_La_max_ratio)
        #high_La_stdev = high_La*((high_La_stdev_ratio/high_La)**2 + (high_La_stdev_Ca/high_La)**2)**0.5
        #print 'Conjugation Length (high): ', '%.3f' %high_La,  '+/- ', '%.3f' %high_La_stdev, 'nm'
        
        Ratio = [Exp_ratio, Exp_ratio]
        La = [low_La, high_La]
        
        fig = plt.figure(3)
        ax = fig.add_subplot(111)
        ax.plot(La_Range , Model_Ratio,'-k')
        ax.fill_between(La_Range, Min_Ratio, Max_Ratio,  facecolor='gray', alpha = 0.25)
        #ax.plot(La , Ratio,'or')
        ax.loglog(La , Ratio,'or')
        high_label = '%.3f' %high_La #+ '+/- ' + '%.3f' %high_La_stdev + 'nm'
        low_label = '%.3f' %low_La #+ '+/- ' + '%.3f' %low_La_stdev + 'nm'
        ax.text(high_La, Exp_ratio, high_label, va="top", ha = "center", size = 10)
        ax.text(low_La, Exp_ratio, low_label, va="top", ha = "center", size = 10)
        ax.set_xlabel('La')
        ax.set_ylabel('ID/IG')
        #ax.set_xlim(0, 10)
        #ax.set_ylim(0, 2)
        ax.set_xlim(0.1, 100)
        ax.set_ylim(0.01, 100)
        plt.savefig(SaveName + '_Ratio.jpg')
        
        #plt.show()
        plt.close()
        
        
        f = open(SaveName+'_fitfile.txt',"w")
        f.write("{}\t{}\t{}\n".format('Collection Details',CollDet,''))
        f.write("{}\t{}\t{}\n".format('Original File',Loadfile,''))
        f.write("{}\t{}\t{}\n".format('position','0','0'))
        f.write("{}\t{}\t{}\n".format('Laser Wavelength', Ext_Lambda,'0'))
        f.write("{}\t{}\t{}\n".format('Num Peaks Fit', NumPeaks,'0'))

        f.write("{}\t{}\t{}\n".format('Baseline Order', base_order, 0) )            
        f.write("{}\t{}\t{}\n".format('Baseline R2', BaseR2[0], 0) )
        f.write("{}\t{}\t{}\n".format('Baseline Slope', BaseSlope, 0) )
        
        f.write("{}\t{}\t{}\n".format('Peak Fit R2?', R2_fit, 0) )
        f.write("{}\t{}\t{}\n".format('Peak Fit SEE', SEE_fit, 0) )
        
        f.write("{}\t{}\t{}\n".format('qBWF', qBWF, 0) )
        
        f.write("{}\t{}\t{}\n".format('G Band Peak Position', fit_results[0], dFit[0]) )
        f.write("{}\t{}\t{}\n".format('G Band Peak Width', fit_results[6], dFit[6] ))
        f.write("{}\t{}\t{}\n".format('G Band Peak Intensity', G_ints,dFit[12]))
        
        f.write("{}\t{}\t{}\n".format('D Band Peak Position',fit_results[1],dFit[1]))
        f.write("{}\t{}\t{}\n".format('D Band Peak Width',fit_results[7],dFit[7]))
        f.write("{}\t{}\t{}\n".format('D Band Peak Intensity', D_ints,dFit[13]))
        
        f.write("{}\t{}\t{}\n".format('D2 Band Peak Position',fit_results[2],dFit[2]))
        f.write("{}\t{}\t{}\n".format('D2 Band Peak Width',fit_results[8],dFit[8]))
        f.write("{}\t{}\t{}\n".format('D2 Band Peak Intensity', D2_ints,dFit[14]))
        
        f.write("{}\t{}\t{}\n".format('D3 Band Peak Position',fit_results[3],dFit[3]))
        f.write("{}\t{}\t{}\n".format('D3 Band Peak Width',fit_results[9],dFit[9]))
        f.write("{}\t{}\t{}\n".format('D3 Band Peak Intensity', D3_ints,dFit[15]))
        
        f.write("{}\t{}\t{}\n".format('D4 Band Peak Position',fit_results[4],dFit[4]))
        f.write("{}\t{}\t{}\n".format('D4 Band Peak Width',fit_results[10],dFit[10]))
        f.write("{}\t{}\t{}\n".format('D4 Band Peak Intensity', D4_ints,dFit[16]))
        
        f.write("{}\t{}\t{}\n".format('Conjugation Length (low)',low_La,0))
        f.write("{}\t{}\t{}\n".format('Conjugation Length (high)',high_La,0))
        f.write("{}\t{}\t{}\n".format('ID/IG',Exp_ratio,Exp_ratio_stdev))

        f.write("{}\t{}\t{}\n".format('ID/total ratio', D_ints/TotalIntensity,0))
        f.write("{}\t{}\t{}\n".format('ID2/total ratio', D2_ints/TotalIntensity,0))
        f.write("{}\t{}\t{}\n".format('ID3/total ratio', D3_ints/TotalIntensity,0))
        f.write("{}\t{}\t{}\n".format('ID4/total ratio', D4_ints/TotalIntensity,0))
        f.write("{}\t{}\t{}\n".format('IG/total ratio',G_ints/TotalIntensity,0))
        f.write("{}\t{}\t{}\n".format('ID2/IG ratio', D2_ints/G_ints,ID2IGstdev))
        f.write("{}\t{}\t{}\n".format('ID3/ID ratio', D3_ints/D_ints,ID3IDstdev))
        f.write("{}\t{}\t{}\n".format('ID4/ID ratio', D4_ints/D_ints,ID4IDstdev))
        f.write("{}\t{}\t{}\n".format('bkd_low',bkd_bounds[0],bkd_bounds[1]))
        f.write("{}\t{}\t{}\n".format('bkd_hi',bkd_bounds[2],bkd_bounds[3])) 
        f.close()

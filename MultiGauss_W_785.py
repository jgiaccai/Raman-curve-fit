## Uses .csv file from Wasatch instrument
## 20191002 expanding bounds based on scatter plots
## 20201207 getting scan details from metadata
## 20230408 altered to allow fitting of any number of Lorentzian peaks, defined lines 29-34
## 20230524 choose baseline fitting order at top 
## position and location info to maintain parity with other mapping/queueing systems
## 20230628 v.3 now with uncertainties!
## 20230919 using adjusted R2 instead of R2 for baseline and signal fits
## 20231210 add in correlation matrix for individual fits and avgCorr calc


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
from uncertainties import ufloat
from uncertainties import unumpy
from uncertainties.umath import sqrt as usqrt
from uncertainties.umath import exp as uexp
from uncertainties.umath import log as ulog

# File Parameters
# =============================================================================
# # Data fitting parameters
# 
FitGOn = 1 # 1 is yes, 0 is no
FitDOn = 1
FitD2On = 1
FitD3On = 1
FitD4On = 1
FitU1On = 0 #unidentified peak but need to include in envelope for 405 etc

fitVersion = 3.01 #changing if there is a change to base fitting subtr or peak fitting or stats calc.  Not for making figures or summarizing data.

base_order = 3 #order of polynomial for bkg fitting, choose 1, 2, or 3
bkd_bounds = [520, 950, 1750, 2000] #low wavelength limits (low, high) and high wavelength limits (low, high)

G_bounds = [1590, 50, 40, 35] # Center wavelength, wavelength limits, HWHM guess, HWHM limits (currently unused)
D_bounds = [1310, 60, 100, 60]
D2_bounds = [1630, 50, 30, 25]
D3_bounds = [1470, 80, 85, 80]
D4_bounds = [1180, 90, 80, 75]
U1_bounds = [1725, 20, 10, 8]  #no physical basis, trying because weird peak in some 405 data
IIM = 0.8 #Initial intensity multiplier for G and D peaks 
qBWF = -10
# =============================================================================

#at some point files got switched so not saving the wavenumber column...grr.
#adding a file that has the wavenumbers in a column

wavefile = "/Users/jennifergiaccai/OneDrive - Smithsonian Institution/PROJECTS/ChineseInkProject/WasatchWavenumRef-20210624-160600-788505-WP-00665.csv"

Ext_Lambda = 000 #nm

TotalNumPeaks = FitGOn + FitDOn + FitD2On + FitD3On + FitD4On + FitU1On
NumPeaks = FitGOn + FitDOn + FitD2On + FitD3On + FitD4On
NumPksApp = str(NumPeaks)+'Gau'
NumParams    = 3*6     #{Number of parameters to fit}
FitParam =np.zeros(NumParams) 
lobounds = np.zeros(18)
hibounds = np.zeros(18)

CollDet = 'Wasatch unk det'  #If need to restart, replace with collection details

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

def Gaussian(xc,w,I):
    global x_fit
    s = ((x_fit - xc)/w)
    return (I)*np.exp((-1*np.log(2)*s**2))

def GaussianWithUnc(xc,w,I):
    global x_fit
    s = ((x_fit - xc)/w)
    return (I)*unumpy.exp((-1*np.log(2)*s**2))

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
    FitParam[14] = (0.4*FitParam[12])*FitD2On
    FitParam[15] = (0.4*FitParam[13])*FitD3On
    FitParam[16] = (0.4*FitParam[13])*FitD4On
    FitParam[17]  = (0.25*FitParam[12])*FitU1On

    Gfit = FitGOn*Gaussian(FitParam[0],FitParam[6],FitParam[12])
    Dfit = FitDOn*Gaussian(FitParam[1],FitParam[7],FitParam[13])
    D2fit = FitD2On*Gaussian(FitParam[2],FitParam[8],FitParam[14])
    D3fit = FitD3On*Gaussian(FitParam[3],FitParam[9],FitParam[15])
    D4fit = FitD4On*Gaussian(FitParam[4],FitParam[10],FitParam[16])
    U1fit = FitU1On*Gaussian(FitParam[5],FitParam[11],FitParam[17])

    ModelFit = Gfit + Dfit + D2fit + D3fit + D4fit + U1fit    
    
    # # figure with initial fit of various peaks
    
    fig = plt.figure(5)
    ax50 = fig.add_subplot(111)
    ax50.plot(x_fit, signal_fit,'.k', label = 'Experimental')
    ax50.plot(x_fit, Gfit,'-g', label = 'G Peak Fit')
    ax50.plot(x_fit, Dfit,'-b', label = 'D Peak Fit')
    ax50.plot(x_fit, D2fit,'-y', label = 'D2 Peak Fit')
    ax50.plot(x_fit, D3fit,'-c', label = 'D3 Peak Fit')
    ax50.plot(x_fit, D4fit,'-m', label = 'D4 Peak Fit')
    ax50.plot(x_fit, U1fit, '.g', label = 'U1 Peak Fit')
    ax50.plot(x_fit, ModelFit,'-r', label = 'Summed Peak Fit')
    ax50.set_xlabel(r'Raman Shift / cm$^-$$^1$')
    plt.autoscale(enable=True, axis='x', tight=True)
    plt.autoscale(enable=True, axis='y', tight=True)
    ax50.set_ylabel('Raman Intensity')
    ax50.set_ylim(min(signal_fit), max(signal_fit)*1.2)
    #plt.text(1075, 14100, 'ink', fontsize=20)
    plt.tick_params(axis='both', which='major')
    #plt.tight_layout()
    plt.savefig(SaveName + '_initialfit.jpg')
    #plt.show()
    plt.close()
    
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
      
    Gfit = FitGOn*Gaussian(EvalSimp[0],EvalSimp[6],EvalSimp[12])
    Dfit = FitDOn*Gaussian(EvalSimp[1],EvalSimp[7],EvalSimp[13])
    D2fit = FitD2On*Gaussian(EvalSimp[2],EvalSimp[8],EvalSimp[14])
    D3fit = FitD3On*Gaussian(EvalSimp[3],EvalSimp[9],EvalSimp[15])
    D4fit = FitD4On*Gaussian(EvalSimp[4],EvalSimp[10],EvalSimp[16])
    U1fit = FitU1On*Gaussian(EvalSimp[5], EvalSimp[11], EvalSimp[17])
    
    FitY = Gfit + Dfit + D2fit + D3fit + D4fit + U1fit
        
    return(FitY)

def FitFuncWithUnc(x_fit, *EvalSimp):      
    
    Gfit = FitGOn*GaussianWithUnc(EvalSimp[0],EvalSimp[6],EvalSimp[12])
    Dfit = FitDOn*GaussianWithUnc(EvalSimp[1],EvalSimp[7],EvalSimp[13])
    D2fit = FitD2On*GaussianWithUnc(EvalSimp[2],EvalSimp[8],EvalSimp[14])
    D3fit = FitD3On*GaussianWithUnc(EvalSimp[3],EvalSimp[9],EvalSimp[15])
    D4fit = FitD4On*GaussianWithUnc(EvalSimp[4],EvalSimp[10],EvalSimp[16])
    U1fit = FitU1On*GaussianWithUnc(EvalSimp[5], EvalSimp[11], EvalSimp[17])
   
    FitY = Gfit + Dfit + D2fit + D3fit + D4fit + U1fit
        
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
        AdjBaseR2 = 1-(1-BaseR2)*(len(baseline_bkgarea)-1)/((len(baseline_bkgarea))-(base_order+1)-1) #wiki definition independent of baseline order fit

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
        plt.setp(ax11, xticks=[600,800,1000,1200,1400,1600,1800,2000])
        ax11.tick_params(direction='in',labelbottom=True,labelleft=True)
        plt.autoscale(enable=True, axis='x', tight=True) 
        plt.ylim(min(BaseResiduals)*1.15,max(BaseResiduals)*1.15)
        
        gs1.update(left=0.13,right=0.96,top=0.95,bottom=0.12) #as percentages of total figure with 1,1 in upper right
        fig.set_size_inches(6, 5) #width, height
        fname = str(SaveName) + '_base.jpg'
        plt.savefig(fname, dpi = 600, bbox_inches='tight')
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
        
        try:
            minres = curve_fit(FitFunc,x_fit,signal_fit,p0=FitParam,method = 'trf', bounds=bounds, full_output=True) 
        except RuntimeError as error:  #if fit is unsuccessful
                print("No fit found")
                print(error)
                fit_results = np.zeros(NumParams)  #shouldn't need this, but just in case, not passing bad fit data
                continue
        else:  #if fit is successful
            fit_results=minres[0]
            covMatrix = minres[1] 
            dFit = np.sqrt(np.diag(covMatrix))
        
        #setting fit results to zero if peak is off
        
        if FitGOn == 0:
            fit_results[0],fit_results[6],fit_results[12] = 0,1,0
        if FitDOn == 0:
            fit_results[1],fit_results[7],fit_results[13] = 0,1,0
        if FitD2On == 0:
            fit_results[2],fit_results[8],fit_results[14] = 0,1,0
        if FitD3On == 0:
            fit_results[3],fit_results[9],fit_results[15] = 0,1,0
        if FitD4On == 0:
            fit_results[4],fit_results[10],fit_results[16] = 0,1,0
        if FitU1On == 0:
            fit_results[5],fit_results[11],fit_results[17] = 0,1,0
        
        # setting intensities using uncertainties library
        (Gloc, Dloc, D2loc, D3loc, D4loc, U1loc, Gwid, Dwid, D2wid, D3wid, D4wid, U1wid, 
             G_ints, D_ints, D2_ints, D3_ints, D4_ints, U1_ints) = uncertainties.correlated_values(fit_results, covMatrix)
        
            
        Gfit_nom = FitGOn*Gaussian(fit_results[0],fit_results[6],fit_results[12])
        Dfit_nom = FitDOn*Gaussian(fit_results[1],fit_results[7],fit_results[13])
        D2fit_nom = FitD2On*Gaussian(fit_results[2],fit_results[8],fit_results[14])
        D3fit_nom = FitD3On*Gaussian(fit_results[3],fit_results[9],fit_results[15])
        D4fit_nom = FitD4On*Gaussian(fit_results[4],fit_results[10],fit_results[16])
        U1fit_nom = FitU1On*Gaussian(fit_results[5],fit_results[11],fit_results[17])
        
        Gfit = FitGOn*GaussianWithUnc(Gloc,Gwid,G_ints)
        Dfit = FitDOn*GaussianWithUnc(Dloc,Dwid,D_ints)
        D2fit = FitD2On*GaussianWithUnc(D2loc,D2wid,D2_ints)
        D3fit = FitD3On*GaussianWithUnc(D3loc,D3wid,D3_ints)
        D4fit = FitD4On*GaussianWithUnc(D4loc,D4wid,D4_ints)
        U1fit = FitU1On*GaussianWithUnc(U1loc,U1wid,U1_ints)

        ModelFit = Gfit + Dfit +D2fit + D3fit + D4fit + U1fit
        
        
        D2_ints = FitD2On*D2_ints
        D3_ints = FitD3On*D3_ints
        D4_ints = FitD4On*D4_ints
        
        TotalIntensity = G_ints + D_ints + D2_ints + D3_ints + D4_ints
       
        Residuals = (signal_fit - ModelFit)
        ss_res = np.sum((Residuals) ** 2)
        ss_tot = np.sum((signal_fit - np.mean(signal_fit)) ** 2)
        R2_fit = 1-ss_res/ss_tot # not meaningful because can't use R2 on Gaussian or Lorentzian fits, so need to use standard error regression
        AdjR2_fit = 1-(1-R2_fit)*(len(ModelFit)-1)/((len(ModelFit))-(NumParams)-1) # wiki definition, adjusted to account for no. variables but still non-linear so ish?
        SEE_fit = usqrt(ss_res/(len(signal_fit)-NumPeaks*3))
        
        # Figure 4 Plot of individual peak fits and total peak fit with experimental
        
        ModelFit_nom = np.fromiter((ModelFit[i].n for i in range(len(ModelFit))),float, count = len(ModelFit))
        ModelFit_unc = np.fromiter((ModelFit[i].s for i in range(len(ModelFit))),float, count = len(ModelFit))
        Residuals_nom = np.fromiter((Residuals[i].n for i in range(len(ModelFit))),float, count = len(ModelFit))

        
        fig = plt.figure(4)
        gs4 = fig.add_gridspec(nrows=2, ncols=1,
                  hspace=0, wspace=0, height_ratios=[1, 5])
        ax40 = fig.add_subplot(gs4[1])
        ax40.plot(x_fit, signal_fit,'.k', label = 'Experimental')
        ax40.plot(x_fit, Gfit_nom,'-g', label = 'G Peak Fit')
        ax40.plot(x_fit, Dfit_nom,'-b', label = 'D Peak Fit')
        if FitD2On == 1:
            ax40.plot(x_fit , D2fit_nom,'-y', label = 'D2 Peak Fit')
        if FitD3On == 1:
            ax40.plot(x_fit, D3fit_nom,'-c', label = 'D3 Peak Fit')
        if FitD4On == 1:
            ax40.plot(x_fit , D4fit_nom,'-m', label = 'D4 Peak Fit')
        if FitU1On == 1:
            ax40.plot(x_fit, U1fit_nom, '-y', label = 'U1 peak fit')
        ax40.plot(x_fit, ModelFit_nom,'-r', label = 'Summed Peak Fit')
        ax40.fill_between(x_fit, ModelFit_nom - ModelFit_unc, ModelFit_nom + ModelFit_unc,  facecolor='red', alpha = 0.25)
        # will need to double the ModelFit_unc in fill line for 95% confidence only one stdev now
        ax40.set_xlabel(r'Raman Shift / cm$^{-1}$', fontsize=16)
        plt.xlim(750,1900)
        plt.autoscale(enable=True, axis='y')
        plt.setp(ax40, xticks=[800,1000,1200,1400,1600,1800])
        ax40.set_ylabel('Raman Intensity', fontsize=16)
        plt.tick_params(axis='both', which='major', labelsize=14)
        
        ax41 = fig.add_subplot(gs4[0])
        ax41.plot(x_fit,Residuals_nom, '.b')
        ax41.axhline(0, linestyle='--', color = 'gray', linewidth = 1)
        ax41.set_ylabel('Residuals')
        plt.setp(ax41, xticks=[800,1000,1200,1400,1600])
        ax41.tick_params(direction='in',labelbottom=False,labelleft=True)
        plt.xlim(750,1900)        
        plt.ylim(min(Residuals_nom)*1.15,max(Residuals_nom)*1.15)
        
        gs4.update(left=0.16,right=0.94,top=0.95,bottom=0.15) #as percentages of total figure with 1,1 in upper right
        fig.set_size_inches(6, 5) #width, height
        plt.savefig(SaveName + '_fit.jpg', dpi = 600, bbox_inches='tight')
        plt.close()
        
        # ID/IG Ratio
        Exp_ratio = D_ints/G_ints #with uncertainties
        Exp_ratio = fit_results[13]/fit_results[12] # w/o uncertainties calc
        
        # Calculation of uncertainties with covariances
        
        Exp_ratio_stdev = Exp_ratio*((dFit[12]/fit_results[12])**2 + (dFit[13]/fit_results[13])**2 - 
                                   (2*(covMatrix[13,12])/fit_results[12]/fit_results[13]))**0.5
# =============================================================================
#         TotIntstdev = (dFit[12]**2 + dFit[13]**2 + dFit[14]**2 + dFit[15]**2 + dFit[16]**2 + 
#                       2*((covMatrix[13,12]) + (covMatrix[13,14]) + 
#                           (covMatrix[13,15]) + (covMatrix[13,16]) + 
#                           (covMatrix[14,12]) + (covMatrix[14,15]) +
#                           (covMatrix[14,16]) + (covMatrix[15,12]) + 
#                           (covMatrix[15,16]) + (covMatrix[16,12])))**0.5
# =============================================================================
        IDIG = D_ints/G_ints  #with uncertainties calc
        ID2IG = D2_ints/G_ints
        ID3ID = D3_ints/D_ints
        ID4ID = D4_ints/D_ints
        IDIT = D_ints/TotalIntensity
        IGIT = G_ints/TotalIntensity
        ID2IT = D2_ints/TotalIntensity
        ID3IT = D3_ints/TotalIntensity
        ID4IT = D4_ints/TotalIntensity
        
        
# =============================================================================
#         Calculating the conjugation length, La
# =============================================================================
        
        # from Herdman and Miller 2011, based on Cancado et al. 2011 and Luchhese et al. 2010

        ra = 3.1 #Cancado has 3.1 nm,  Luchhese has 3.00 +/- 0.03 nm
        rs = 1.0 # Luchhese has 1.00 +/- 0.04 nm, Cancado also uses 1.0
        CA_mult = ufloat(160,48)# 160 +/- 48 from Cancado et al.
        CA = CA_mult*((1240*Ext_Lambda**-1)**-4)  # unitless
        La_model = np.arange(0.1, 100, step = 0.001)
        Ratio_model = [0]*len(La_model)
        Ratio_model =CA*(ra**2 - rs**2)/(ra**2 - 2*(rs**2))*(np.exp((-np.pi*rs**2)/(La_model**2))-np.exp((-np.pi*(ra**2 - rs**2))/(La_model**2))) 
        
        # find the La from ID/IG values for the lower part of the curve and the higher part of the curve
        # no equations from ID/IG for high defect frequency so have to do it the hard way
                
        Ratio_model_low = np.where(La_model <=3, Ratio_model, np.nan)
        La_low = np.where(La_model <=3, La_model, np.nan)
        low_La = La_low[np.where(abs(Ratio_model_low - Exp_ratio) == min(abs(Ratio_model_low - Exp_ratio)))[0][0]]

        Ratio_model_high = np.where(La_model >3, Ratio_model, 100)
        La_high = np.where(La_model >3, La_model, np.nan)
        high_La = La_high[np.where(abs(Ratio_model_high - Exp_ratio) == min(abs(Ratio_model_high - Exp_ratio)))[0][0]]
        high_La = La_high[np.where(abs(Ratio_model_high - Exp_ratio) == min(abs(Ratio_model_high - Exp_ratio)))[0][0]]
        
        # when resolving eq from Herdman and Miller for LsubA, if LsubA approaches 1, 
        # we can neglect second exponential term as it gets very small compared to first term (by e-10)
        # however for LsubA close to 20, the two terms are close in value, so uncertainty only for low_La
        
        low_La_calc = usqrt((-1*np.pi)/(ulog(((ra**2-2)/(ra**2-1))*(IDIG/CA))))

        Ratio = [Exp_ratio, Exp_ratio]
        La_calc = [low_La, high_La]
        
        if high_La > 8:  # should be >=10 via Cancado et al. Eq 2 but when taking in mind uncertainty...okay
            high_La_calc = usqrt(1/((IDIG)/CA/(np.pi*(3.1**2-1))))
            La_calc2 = [low_La_calc.n, high_La_calc.n]
            high_label = u'{:.2fP}'.format(high_La_calc)
        else:
            La_calc2 = [low_La_calc.n, high_La]
            high_label = '{:.2f}'.format(high_La)
            high_La_calc = ufloat(high_La, np.nan)


        
# =============================================================================
        Ratio = [Exp_ratio, Exp_ratio]
        La_calc = [low_La, high_La]
        La_calc2 = [low_La_calc.n, high_La_calc.n]
        
        Ratio_model_plot = np.fromiter((Ratio_model[i].n for i in range(len(Ratio_model))),float, count = len(Ratio_model))
        dRatio_model = np.fromiter((Ratio_model[i].s for i in range(len(Ratio_model))),float, count = len(Ratio_model))
        
        
        # fig = plt.figure(3)
        # ax = fig.add_subplot(111)
        # ax.plot(La_model , Ratio_model_plot ,'-k')
        # ax.fill_between(La_model, Ratio_model_plot - dRatio_model, Ratio_model_plot + dRatio_model,  facecolor='gray', alpha = 0.25)
        # # will need to double the dRatio_model in fill line for 95% confidence only one stdev now
        # ax.loglog(La_calc2 , Ratio,'or')
        # high_label = u'{:.2fP}'.format(high_La_calc)
        # low_label = u'{:.2fP}'.format(low_La_calc)
        # ax.text(high_La, Exp_ratio-0.5, high_label, va="top", ha = "center", size = 10) #offset for legibility
        # ax.text(low_La, Exp_ratio-0.5, low_label, va="top", ha = "center", size = 10) #offset for legibility
        # ax.set_xlabel('$L_a$', fontsize=16)
        # ax.set_ylabel('$I_D$ / $I_G$', fontsize=16)
        # ax.set_xlim(0.1, 100)
        # ax.set_ylim(0.01, 100)
        # fig.set_size_inches(6, 5) #width, height
        # plt.savefig(SaveName + '_Ratio.jpg',dpi=300, bbox_inches='tight')
        
        # plt.close()
        
        
        f = open(SaveName+'_fitfile.txt',"w")
        f.write("{}\t{}\t{}\n".format('Collection Details',CollDet,''))
        f.write("{}\t{}\t{}\n".format('Original File',Loadfile,''))
        f.write("{}\t{}\t{}\n".format('Location','0','0'))
        f.write("{}\t{}\t{}\n".format('position','0','0'))
        f.write("{}\t{}\t{}\n".format('Laser Wavelength', Ext_Lambda,'0'))
        f.write("{}\t{}\t{}\n".format('Num Peaks Fit and fit version', NumPeaks,fitVersion))

        f.write("{}\t{}\t{}\n".format('Baseline Order', base_order, 0) )            
        f.write("{}\t{}\t{}\n".format('Baseline R2', AdjBaseR2[0], 0) )
        f.write("{}\t{}\t{}\n".format('Baseline Slope', BaseSlope, 0) )
        
        f.write("{}\t{}\t{}\n".format('Peak Fit R2ish', float(AdjR2_fit.n), float(AdjR2_fit.s)) )
        f.write("{}\t{}\t{}\n".format('Peak Fit SEE', float(SEE_fit.n), float(SEE_fit.s)) )
        
        f.write("{}\t{}\t{}\n".format('qBWF', 0, 0) )
        
        f.write("{}\t{}\t{}\n".format('G Band Peak Position', fit_results[0], dFit[0]) )
        f.write("{}\t{}\t{}\n".format('G Band Peak Width', fit_results[6], dFit[6] ))
        f.write("{}\t{}\t{}\n".format('G Band Peak Intensity', fit_results[12],dFit[12]))
        
        f.write("{}\t{}\t{}\n".format('D Band Peak Position',fit_results[1],dFit[1]))
        f.write("{}\t{}\t{}\n".format('D Band Peak Width',fit_results[7],dFit[7]))
        f.write("{}\t{}\t{}\n".format('D Band Peak Intensity', fit_results[13],dFit[13]))
        
        f.write("{}\t{}\t{}\n".format('D2 Band Peak Position',fit_results[2],dFit[2]))
        f.write("{}\t{}\t{}\n".format('D2 Band Peak Width',fit_results[8],dFit[8]))
        f.write("{}\t{}\t{}\n".format('D2 Band Peak Intensity', fit_results[14],dFit[14]))
        
        f.write("{}\t{}\t{}\n".format('D3 Band Peak Position',fit_results[3],dFit[3]))
        f.write("{}\t{}\t{}\n".format('D3 Band Peak Width',fit_results[9],dFit[9]))
        f.write("{}\t{}\t{}\n".format('D3 Band Peak Intensity', fit_results[15],dFit[15]))
        
        f.write("{}\t{}\t{}\n".format('D4 Band Peak Position',fit_results[4],dFit[4]))
        f.write("{}\t{}\t{}\n".format('D4 Band Peak Width',fit_results[10],dFit[10]))
        f.write("{}\t{}\t{}\n".format('D4 Band Peak Intensity', fit_results[16],dFit[16]))
        
        f.write("{}\t{}\t{}\n".format('Conjugation Length (low)',low_La_calc.n,low_La_calc.s))
        f.write("{}\t{}\t{}\n".format('Conjugation Length (high)',high_La,high_La_calc.s))
        f.write("{}\t{}\t{}\n".format('ID/IG',Exp_ratio,Exp_ratio_stdev))

        f.write("{}\t{}\t{}\n".format('ID/total ratio', IDIT.n,IDIT.s))
        f.write("{}\t{}\t{}\n".format('ID2/total ratio', ID2IT.n,ID2IT.s))
        f.write("{}\t{}\t{}\n".format('ID3/total ratio', ID3IT.n,ID3IT.s))
        f.write("{}\t{}\t{}\n".format('ID4/total ratio', ID4IT.n,ID4IT.s))
        f.write("{}\t{}\t{}\n".format('IG/total ratio',IGIT.n,IGIT.s))
        f.write("{}\t{}\t{}\n".format('ID2/IG ratio', ID2IG.n,ID2IG.s))
        f.write("{}\t{}\t{}\n".format('ID3/ID ratio', ID3ID.n,ID3ID.s))
        f.write("{}\t{}\t{}\n".format('ID4/ID ratio', ID4ID.n,ID4ID.s))
        f.write("{}\t{}\t{}\n".format('bkd_low',bkd_bounds[0],bkd_bounds[1]))
        f.write("{}\t{}\t{}\n".format('bkd_hi',bkd_bounds[2],bkd_bounds[3])) 
        f.close()

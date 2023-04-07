## Uses .txt file from NMNH mapping results for 405 nm
## need to add section determining mapping versus individual data text files
## need to add section writing original name of data file to the summary text file and add into summary prog
## 20191002 expanding bounds based on scatter plots
## 20201207 getting scan details from metadata

import sys
import numpy as np
#from pylab import *
#from scipy.optimize import curve_fit, fmin
from scipy.optimize import minimize
from scipy.optimize import least_squares
import matplotlib.pyplot as plt
from scipy import signal
#from scipy.signal import medfilt
from scipy.interpolate import interp1d
import os
import fnmatch
import random
import math
from time import sleep
#from pandas import *
import pandas as pd

# File Parameters
# Data fitting parameters

FitGOn = 1 # 1 is yes, 0 is no
FitDOn = 1
FitD2On = 0
FitD3On = 0
FitD4On = 0
FitU1On = 1 #unidentified peak but need to include in envelope

fitVersion = 2.0

NumPeaks = FitGOn + FitDOn + FitD2On + FitD3On + FitD4On + FitU1On
NumParams    = 3*NumPeaks      #{Number of parameters to fit}
FitParam =np.zeros(NumParams) 
bounds = np.zeros((NumParams,2))

bkd_bounds = [950, 1020, 1800, 2050] #low wavelength limits (low, high) and high wavelength limits (low, high)
G_bounds = [1590, 35, 50, 40] # Center wavelength, wavelength limits, HWHM guess, HWHM limits (currently unused)
D_bounds = [1350, 60, 100, 40]
D2_bounds = [1620, 10, 20, 10]
D3_bounds = [1500, 10, 45, 40]
D4_bounds = [1200, 10, 60, 40]
U1_bounds = [1725, 20, 10, 8]  #no physical basis, trying because weird peak in 405 data
IIM = 0.8 #Initial intensity multiplier for G and D peaks 
qBWF = -10
Ext_Lambda = 000 #nm
step_size = 1 #nm x-wave spacing for final data
iterations = 50

lowG = G_bounds[0] - (G_bounds[1]/2)
highG = G_bounds[0] + (G_bounds[1]/2)
lowD = D_bounds[0] - (D_bounds[1]/2)
highD = D_bounds[0] + (D_bounds[1]/2)
    
lowD = 1100
CollDet = 'NMNH unk det'  #If need to restart, replace with collection details

def ReadCollDetails(inputFilename):
    global CollDet, Ext_Lambda
    lineList = []
    i=0
    #txtFile = open(inputFilename, "r")
    with open(inputFilename, encoding="utf8", errors='ignore') as txtFile:
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
            i = i+1
    CollDet = 'NMNH '+acqtime+'s '+NoScans+ 'sc ' + Laser + 'lp ' +lp
    txtFile.close()
    return CollDet,Ext_Lambda

def lorentz(xc, w, I):  
    global x_fit
    s = ((x_fit - xc)/w)
    return (I)*(1/(1+s**2)) #Wikipedia definition.  gives correct value for intensity.  using this definition fit_results you don't need to account for the peak width in the intensity.

def BWF(xc,w,I):
    global qBWF,x_fit
    s = ((x_fit - xc)/w)
    return (I)*((1+s/qBWF)**2/(1+s**2))

def EnterData():
    
    global FitParam ,G_ints, D_ints, D5_ints, NumParams, G_bounds, D_bounds, D5_bounds, bounds

    FitParam[0] = G_bounds[0] # G peak position
    FitParam[1]  = D_bounds[0] # D peak position
    FitParam[2] = D5_bounds[0] # D5 peak position

    FitParam[3]  = G_bounds[2]  # G peak width
    FitParam[4]  = D_bounds[2]  # D peak width
    FitParam[5] = D5_bounds[2]  # D5 peak width

    '''need to work G and D starting peak intensities from initial values'''
    FitParam[6]  = IIM*G_ints  # G peak intensity  
    FitParam[7]  = IIM*D_ints  # D peak intensity
    FitParam[8] = 0.25 * FitParam[6]  # D5 peak intensity relative to G peak intensity
    
    Gfit = BWF(FitParam[0],FitParam[3],FitParam[6])
    Dfit = lorentz(FitParam[1],FitParam[4],FitParam[7])
    D5fit = lorentz(FitParam[2],FitParam[5],FitParam[8])
    
    ModelFit = Gfit + Dfit +D5fit
    
    
    # fig = plt.figure(5)
    # ax50 = fig.add_subplot(111)
    # ax50.plot(x_fit, signal_fit,'.k', label = 'Experimental')
    # ax50.plot(x_fit, Gfit,'-g', label = 'G Peak Fit')
    # ax50.plot(x_fit, Dfit,'-b', label = 'D Peak Fit')
    # ax50.plot(x_fit, D5fit, '-y', label = 'D5 peak fit')
    # ax50.plot(x_fit, ModelFit,'-r', label = 'Summed Peak Fit')
    # ax50.set_xlabel(r'Raman Shift / cm$^-$$^1$', fontsize=18)
    # plt.autoscale(enable=True, axis='x', tight=True)
    # plt.autoscale(enable=True, axis='y', tight=True)
    # ax50.set_ylabel('Raman Intensity', fontsize=18)
    # ax50.set_ylim(0, max(signal_fit)*1.2)
    # plt.text(1075, 14100, 'ink', fontsize=20)
    # plt.tick_params(axis='both', which='major', labelsize=18)
    # plt.tight_layout()
    # #ax40.legend(bbox_to_anchor=(0.00, 0.70, .3, .152), loc=3,
    # #        ncol=1, mode='expand', frameon=True, borderaxespad=0.)
    # #
    # 
    # 
    # plt.savefig(SaveName + 'initialfit_3BWFLS.jpg')
    # #plt.show()
    # plt.close()
    
        
    
    bounds[0] = (G_bounds[0]-G_bounds[1],G_bounds[0]+G_bounds[1])
    bounds[1] = (D_bounds[0]-D_bounds[1],D_bounds[0]+D_bounds[1])
    bounds[2] = (D5_bounds[0]-D5_bounds[1],D5_bounds[0]+D5_bounds[1])
    
    bounds[3] = (G_bounds[2]-G_bounds[3],G_bounds[2]+G_bounds[3])
    bounds[4] = (D_bounds[2]-D_bounds[3],D_bounds[2]+D_bounds[3])
    bounds[5] = (D5_bounds[2]-D5_bounds[3],D5_bounds[2]+D5_bounds[3])
    
    bounds[6] = (0.5*G_ints,1.2*G_ints)
    bounds[7] = (0.5*D_ints,1.2*D_ints)
    bounds[8] = (1,0.5*G_ints)
    
    #print('Bounds are ',bounds)
    
def Evaluate(EvalSimp):   
    global NumParams, G_bounds, D_bounds, D5_bounds
    global x_fit, signal_fit, Residuals
    
    '''
    Need to 
    (1) evaluate Lorentzian for each peak
    (2) Add all peak fits together for total peak fit
    (3) Subtract total peak fit from real data for initial residuals
    (4) punish residuals if various parameters are too far from real
    '''
    
    Gfit = BWF(EvalSimp[0],EvalSimp[3],EvalSimp[6])
    Dfit = lorentz(EvalSimp[1],EvalSimp[4],EvalSimp[7])
    D5fit = lorentz(EvalSimp[2],EvalSimp[5],EvalSimp[8])
    
    EvalFit = Gfit + Dfit + D5fit
    
    Residuals = (EvalFit - signal_fit)
    ErrorSum = np.linalg.norm(Residuals)
    
    '''Residual punishment--similar to JDH
    (1) D2/D3/D4 intensity is higher than G or D intensity
    (2) peak locations vary too much from my ideal
    (3) peak intensity is negative
    (4) peak width is negative
    (5) between 1200 and 1700 signal fit is more than 1% different than true signal
    '''
    ErrorWeight = 1
    # CheckThese = [1290,1300,1310,1320,1330,1340,1350,1360,1470,1490,1520,1550,1560,1570,1580,1590,1600]
    # for n in CheckThese:
    #     if Residuals[n-925] > 0.1*signal_fit[n-925]: #xvalue of resds is wavenumber - 925
    #         ErrorWeight = ErrorWeight*(100*Residuals[n-925])/signal_fit[n-925]
    #         #print n, 'is ',Residuals[n-925]/signal_fit[n-925]
    # 
    # for n in range(0,len(EvalSimp)):
    #     if EvalSimp[n] <0:
    #         ErrorWeight = ErrorWeight *10
    #         #print 'cond 1'
    #     if n == 0:
    #         if abs(1580 - EvalSimp[0])>G_bounds[1]:
    #             ErrorWeight = ErrorWeight *10
    #             #print 'cond 3'
    #     if n == 1:
    #         if abs(1350 - EvalSimp[1])>D_bounds[1]:
    #             ErrorWeight = ErrorWeight *10
    #             #print 'cond 4'
    # 
    # #print 'ErrorWeight for array ', iterations,' is ',ErrorWeight
    WeightResiduals = ErrorSum*ErrorWeight
    
    return(WeightResiduals)




for file in os.listdir('.'):
    if fnmatch.fnmatch(file, '*.txt'):
        print(file)
        Loadfile = file
        filename = Loadfile[:-4]
        RefMat = Loadfile[:-4]
        
        (CollDet,Ext_Lambda) = ReadCollDetails(file)
        Ext_Lambda = int(Ext_Lambda)
        print(Ext_Lambda)
        
        all_data = pd.read_csv(Loadfile, skiprows = (38), sep = None, header = None, engine='python')
        locations = all_data[0]
        all_data = all_data.transpose()
        
        xy_data = all_data[1:]
        
        wavenum = xy_data[0]
        for n in range(1,len(xy_data.columns)):
            signal = xy_data[n]
            SpotNo = n
            print(n)
            SaveName = 'results/' + RefMat +'_POS'+ str(n).zfill(2)
            
            #Fit Background and cut-down region wavenum and signal to region of interest
            x_fit = wavenum[(bkd_bounds[0] <= wavenum) & (wavenum <= bkd_bounds[3])] # X data cut down to bkg boundaries
            signal_fit = signal[(bkd_bounds[0] <= wavenum) & (wavenum <= bkd_bounds[3])] # y data cut down to bkg boundaries
            
            
            wavesforfit =  np.arange(bkd_bounds[0], bkd_bounds[3], step_size) # x data between bounds at step size
            tck = interp1d(x_fit, signal_fit, bounds_error=False, fill_value=0) # 
            signal_fit = tck(wavesforfit) # y data between bounds at step size
            x_fit = wavesforfit 
            
            lower_bkg = np.where((bkd_bounds[0] <= x_fit) & (x_fit <= bkd_bounds[1]))[0] # x values low part of bkg
            higher_bkg = np.where((bkd_bounds[2] <= x_fit) & (x_fit <= bkd_bounds[3]))[0] # x values high part of bkg
            
            bkg_x = [0]*(len(lower_bkg)+len(higher_bkg)) # wave length of bkg only for x bkg values in loops below
            bkg_signal= [0]*(len(lower_bkg)+len(higher_bkg)) # wave length of bkg only for y bkg values
            count = 0
            for i in range(0, len(lower_bkg)):
                bkg_x[count] = x_fit[lower_bkg[i]]
                bkg_signal[count] = signal_fit[lower_bkg[i]]   
                count += 1
            for j in range(0, len(higher_bkg)):
                bkg_x[count] = x_fit[higher_bkg[j]]
                bkg_signal[count] = signal_fit[higher_bkg[j]]   
                count += 1
                
            ## Option for Linear fit for baseline fit
            
            # BasePara, residuals, rank, singular_values, rcond = np.polyfit(bkg_x,bkg_signal,1, full = True)
            # m = BasePara[0]
            # b = BasePara[1]
            # baseline = np.polyval([m,b],x_fit)
            #print 'Baseline: ', '%.5f' %m, 'x + ', '%.5f' %b
            #print 'R2 Baseline: ', 1-(residuals[0]/sum(bkg_signal))
            
            ## Option for 2d Order polynomial fit for baseline fit
            BasePara, residuals, rank, singular_values, rcond  =  np.polyfit(bkg_x,bkg_signal,2, full = True)
            a = BasePara[0]
            b = BasePara[1]
            c = BasePara[2]
            baseline = np.polyval([a,b,c],x_fit)
            
            
            ## Option for 3rd Order polynomial fit for baseline fit
            # BasePara, residuals, rank, singular_values, rcond  =  np.polyfit(bkg_x,bkg_signal,3, full = True)
            # a = BasePara[0]
            # b = BasePara[1]
            # c = BasePara[2]
            # d = BasePara[3]
            # baseline = np.polyval([a,b,c,d],x_fit)
            ##print 'Baseline: ', a, 'x**3 + ', b, 'x**2 + ', c , 'x + ', d
            ##print 'R2 Baseline: ', 1-(residuals[0]/sum(bkg_signal))
            
            fig = plt.figure(1)
            ax = fig.add_subplot(111)
            ax.plot(wavenum , signal,'-k')
            ax.plot(x_fit, baseline, '-r')
            ax.plot(bkg_x , bkg_signal,'.b')
            ax.set_xlabel('Raman Shift (cm-1)')
            ax.set_ylabel('Raman Intensity')
            ax.set_ylim(0, 1.5*max(signal_fit))
            #plt.show()
            fname = str(SaveName) + '_base_3LS.jpg'
            plt.savefig(fname)
            plt.close()
            
            # Baseline Correction
            signal_fit = signal_fit - baseline #corrects for baseline 
            Signal_hi = signal[np.where(abs(wavenum - 1910) ==min(abs(wavenum - 1910)))[0][0]]
            #Signal_2500 = signal[np.where(abs(wavenum - 2500) ==min(abs(wavenum - 2500)))[0][0]]
            #Signal_800 = signal[np.where(abs(wavenum - 800) ==min(abs(wavenum - 800)))[0][0]]
            Signal_lo = signal[np.where(abs(wavenum - 925) ==min(abs(wavenum - 925)))[0][0]]
            #print 'Baseline Flatness:', Signal_hi/Signal_lo #NOTE if there's noise right here this value will be not so hot
            
            
            #Finds an initial G intensity and D intensity to give a decent start to the peak fitting
            G_ints = max(signal_fit[(G_bounds[0]-G_bounds[1] <= x_fit) & (x_fit <= G_bounds[0]+G_bounds[1])])
            D_ints = max(signal_fit[(D_bounds[0]-D_bounds[1] <= x_fit) & (x_fit <= D_bounds[0]+D_bounds[1])]) 
            D5_ints = max(signal_fit[(D5_bounds[0]-D5_bounds[1] <= x_fit) & (x_fit <= D5_bounds[0]+D5_bounds[1])])
            #print('D5 raw intensity ',D5_ints)
            #print 'initial D/G intensity ratio: ',D_ints/G_ints 
            
            
            EnterData()
            #print 'FitParam ', FitParam
            minres = minimize(Evaluate,FitParam,method='SLSQP',tol=1e-8,bounds=bounds,options={'maxiter':1e5}) 
            if minres.success == 1:
                fit_results=minres.x
                #print(fit_results)
                iterations = minres.nit
            else:
                print(minres.message)
                fit_results = np.zeros(NumParams)
                continue
            Gfit = BWF(fit_results[0],fit_results[3],fit_results[6])
            Dfit = lorentz(fit_results[1],fit_results[4],fit_results[7])
            D5fit = lorentz(fit_results[2],fit_results[5],fit_results[8])
                
            ModelFit = Gfit + Dfit +D5fit
            
            G_ints = fit_results[6]  #intensities from the fit, not the initial values
            D_ints = fit_results[7]
            D5_ints = fit_results[8]
                                    
            #fit_error = np.corrcoef(ModelFit, signal_fit)[0][1] # Calculates R2 from Erin June 2016
            ss_res = np.sum((Residuals) ** 2)
            ss_tot = np.sum((signal_fit - np.mean(signal_fit)) ** 2)
            r2 = 1 - (ss_res / ss_tot)
            #print(r2)
            fit_error = r2  #Coefficient of Determination to determine r-squared
            
            # Figure 4 Plot of individual peak fits and total peak fit with experimental
            fig = plt.figure(4)
            ax40 = fig.add_subplot(111)
            ax40.plot(x_fit, signal_fit,'.k', label = 'Experimental')
            ax40.plot(x_fit, Gfit,'-g', label = 'G Peak Fit')
            ax40.plot(x_fit, Dfit,'-b', label = 'D Peak Fit')
            ax40.plot(x_fit, D5fit, '-y', label = 'D5 peak fit')
            ax40.plot(x_fit, ModelFit,'-r', label = 'Summed Peak Fit')
            ax40.set_xlabel(r'Raman Shift / cm$^{-1}$', fontsize=18)
            plt.autoscale(enable=True, axis='x', tight=True)
            plt.autoscale(enable=True, axis='y', tight=True)
            ax40.set_ylabel('Raman Intensity', fontsize=18)
            ax40.set_ylim(0, max(signal_fit)*1.2)
            plt.text(1075, 14100, 'ink', fontsize=20)
            plt.tick_params(axis='both', which='major', labelsize=18)
            plt.tight_layout()
            #ax40.legend(bbox_to_anchor=(0.00, 0.70, .3, .152), loc=3,
            #        ncol=1, mode='expand', frameon=True, borderaxespad=0.)
            #
            # ax40 = fig.add_subplot(111)
            ax41 = plt.axes([0.185, .88, .787, .1])
            ax41.plot(x_fit,Residuals, '.b')
            #ax41.set_ylabel('Residuals')
            plt.setp(ax41, xticks=[1000,1200,1400,1600,1800], yticks=[-100,0,100])
            ax41.tick_params(direction='in',labelbottom=False,labelleft=False)
            plt.autoscale(enable=True, axis='x', tight=True) 
            
            #plt.ylim(min(Residuals)*1.15,max(Residuals)*1.15)
            plt.ylim(min(Residuals)*1.15,130)
            
            plt.savefig(SaveName + 'fit_3BWFLS.jpg')
            #plt.show()
            plt.close()
                    
            # Uncertainty analysis
            G_error = [0]*iterations
            D_error = [0]*iterations
            for i in range(0, iterations):
                model_error = ModelFit + fit_error*np.random.randn(1)/4*ModelFit
                G_error[i] =model_error[np.where(x_fit == np.around(fit_results[0],0) )[0][0]]
                D_error[i] = model_error[np.where(x_fit == np.around(fit_results[1],0) )[0][0]]
            
            G_ints_error =  np.std(G_error) 
            D_ints_error =  np.std(D_error)
                
            #print 'G Band Peak Position: ', fit_results[0]
            #print 'G Band Peak Width: ', fit_results[2]
            #print 'G Band Peak Intensity: ', '%.3f' %G_ints#, '+/- ',  '%.3f' %G_ints_error
            #print 'D Band Peak Position: ', fit_results[1]
            #print 'D Band Peak Width: ', fit_results[3]
            #print 'D Band Peak Intensity: ', '%.3f' %D_ints  # , '+/- ',  '%.3f' %D_ints_error 
            
            # ID/IG Ratio
            Exp_ratio = D_ints/G_ints
            Exp_ratio_err = Exp_ratio*((G_ints_error/G_ints)**2 + (D_ints_error/D_ints)**2)**0.5
            
            #print 'ID/IG Ratio: ', '%.3f' %Exp_ratio #,'+/- ',  '%.3f' % Exp_ratio_err 
            #print 'q = ',qBWF
            #print 'Sum Sq Residuals = ',np.sum(Residuals**2)
            
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
            low_La_error_Ca = abs(low_La_min_Ca - low_La_max_Ca)
            
            low_La_min_ratio = La_low[np.where(abs(Model_Ratio_low - (Exp_ratio-Exp_ratio_err) ) == min(abs(Model_Ratio_low - (Exp_ratio-Exp_ratio_err))))[0][0]]
            low_La_max_ratio = La_low[np.where(abs(Model_Ratio_low - (Exp_ratio+Exp_ratio_err)) == min(abs(Model_Ratio_low - (Exp_ratio+Exp_ratio_err))))[0][0]]
            low_La_error_ratio = abs(low_La_min_ratio - low_La_max_ratio)
            low_La_error = low_La*((low_La_error_ratio/low_La)**2 + (low_La_error_Ca/low_La)**2)**0.5
            #print 'Conjugation Length (low): ', '%.3f' %low_La,  '+/- ', '%.3f' %low_La_error, 'nm'
            
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
            high_La_error_Ca = abs(high_La_min_Ca - high_La_max_Ca)
            
            high_La_min_ratio = La_high[np.where(abs(Model_Ratio_high - (Exp_ratio-Exp_ratio_err) ) == min(abs(Model_Ratio_high - (Exp_ratio-Exp_ratio_err))))[0][0]]
            high_La_max_ratio = La_high[np.where(abs(Model_Ratio_high - (Exp_ratio+Exp_ratio_err)) == min(abs(Model_Ratio_high - (Exp_ratio+Exp_ratio_err))))[0][0]]
            high_La_error_ratio = abs(high_La_min_ratio - high_La_max_ratio)
            high_La_error = high_La*((high_La_error_ratio/high_La)**2 + (high_La_error_Ca/high_La)**2)**0.5
            #print 'Conjugation Length (high): ', '%.3f' %high_La,  '+/- ', '%.3f' %high_La_error, 'nm'
            
            Ratio = [Exp_ratio, Exp_ratio]
            La = [low_La, high_La]
            
            fig = plt.figure(3)
            ax = fig.add_subplot(111)
            ax.plot(La_Range , Model_Ratio,'-k')
            ax.fill_between(La_Range, Min_Ratio, Max_Ratio,  facecolor='gray', alpha = 0.25)
            #ax.plot(La , Ratio,'or')
            ax.loglog(La , Ratio,'or')
            high_label = '%.3f' %high_La #+ '+/- ' + '%.3f' %high_La_error + 'nm'
            low_label = '%.3f' %low_La #+ '+/- ' + '%.3f' %low_La_error + 'nm'
            ax.text(high_La, Exp_ratio, high_label, va="top", ha = "center", size = 10)
            ax.text(low_La, Exp_ratio, low_label, va="top", ha = "center", size = 10)
            ax.set_xlabel('La')
            ax.set_ylabel('ID/IG')
            #ax.set_xlim(0, 10)
            #ax.set_ylim(0, 2)
            ax.set_xlim(0.1, 100)
            ax.set_ylim(0.01, 100)
            plt.savefig(SaveName + '_Ratio_3BWFLS.jpg')
            
            #plt.show()
            plt.close()
            
            
            f = open(SaveName+'_3BWFLS_fitfile.txt',"w")
            f.write("{}\t{}\t{}\n".format('Collection Details',CollDet,''))
            f.write("{}\t{}\t{}\n".format('Original File',Loadfile,''))
            f.write("{}\t{}\t{}\n".format('Baseline R2??', 1-(residuals[0]/sum(bkg_signal)), 0) )
            f.write("{}\t{}\t{}\n".format('Baseline Flatness??', Signal_hi/Signal_lo, 0) )
            f.write("{}\t{}\t{}\n".format('Model Fit R2', r2, 0) )
            f.write("{}\t{}\t{}\n".format('Q factor', qBWF, 0) )
            
            f.write("{}\t{}\t{}\n".format('G Band Peak Position', fit_results[0], 0) )
            f.write("{}\t{}\t{}\n".format('G Band Peak Width', fit_results[3], 0 ))
            f.write("{}\t{}\t{}\n".format('G Band Peak Intensity', fit_results[6],0))
            
            f.write("{}\t{}\t{}\n".format('D Band Peak Position',fit_results[1],0))
            f.write("{}\t{}\t{}\n".format('D Band Peak Width',fit_results[4],0))
            f.write("{}\t{}\t{}\n".format('D Band Peak Intensity',fit_results[7],0))
            
            f.write("{}\t{}\t{}\n".format('D5 Band Peak Position',fit_results[2],0))
            f.write("{}\t{}\t{}\n".format('D5 Band Peak Width',fit_results[5],0))
            f.write("{}\t{}\t{}\n".format('D5 Band Peak Intensity',fit_results[8],0))
            
            f.write("{}\t{}\t{}\n".format('Conjugation Length (low)',low_La,0))
            f.write("{}\t{}\t{}\n".format('Conjugation Length (high)',high_La,0))
            f.write("{}\t{}\t{}\n".format('ID/IG',Exp_ratio,0))
            f.write("{}\t{}\t{}\n".format('Location', locations[n],''))
            f.write("{}\t{}\t{}\n".format('Laser Wavelength', Ext_Lambda,''))
            f.close()

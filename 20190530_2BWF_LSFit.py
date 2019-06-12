## Uses .xls file from KOS Raman system, auto Export to Excel 
## Collection details go in X-column, Reference Material description from filename
## Auto increments the spot location, all with same distance between Raman spots (TravDist)

import sys
import numpy as np
from pylab import *
#from scipy.optimize import curve_fit, fmin
from scipy.optimize import minimize
from scipy.optimize import least_squares
import matplotlib.pyplot as plt
from scipy import signal
#from scipy.signal import medfilt
from scipy.interpolate import interp1d
import random
import math
from time import sleep
from pandas import *
import pandas as pd

# File Parameters

Loadfile =  'Kob_3dotJS.xls'
RefMat = Loadfile[:-4]
TravDist = 0.20 #mm distance stage is moving between points
CollDet = 'not supplied'  #If need to restart, replace with collection details
Location = 0

# Data fitting parameters

NumPeaks = 2
NumParams    = 3*NumPeaks      #{Number of parameters to fit}
FitParam =np.zeros(NumParams) 
bounds = np.zeros((NumParams,2))

bkd_bounds = [925, 1000, 1800, 1870] #low wavelength limits (low, high) and high wavelength limits (low, high)
G_bounds = [1580, 15, 50, 40] # Center wavelength, wavelength limits, HWHM guess, HWHM limits (currently unused)
D_bounds = [1350, 20, 120, 40]
D2_bounds = [1620, 10, 20, 10]
D3_bounds = [1500, 10, 45, 40]
D4_bounds = [1200, 10, 60, 40]
IIM = 0.6 #Initial intensity multiplier for G and D peaks 
qBWF = -10
Ext_Lambda = 785 #nm
step_size = 1 #nm x-wave spacing for final data
iterations = 500

lowG = G_bounds[0] - (G_bounds[1]/2)
highG = G_bounds[0] + (G_bounds[1]/2)
lowD = D_bounds[0] - (D_bounds[1]/2)
highD = D_bounds[0] + (D_bounds[1]/2)
    
lowD = 1100

all_data = pd.read_excel(Loadfile, skiprows = (9), header = None)

for n in range(2,len(all_data.columns)/2):
    SpotNo = n
    print n
    SaveName = RefMat +'_POS'+ str(n).zfill(2)
    wave = all_data[2*n].values
    signal = all_data[(2*n)+1].values
    
    new_info = pd.read_excel(Loadfile, dtype = str, usecols = (2*n,), skiprows = (8), skipfooter = (len(all_data)), header = None)
    if len(new_info.index) <> 0:
        CollDet = str(new_info[0].values)[2:-2]
        print CollDet
    Location = Location + TravDist
    
    #Fit Background and cut-down region wave and signal to region of interest
    x_fit = wave[(bkd_bounds[0] <= wave) & (wave <= bkd_bounds[3])] # X data cut down to bkg boundaries
    signal_fit = signal[(bkd_bounds[0] <= wave) & (wave <= bkd_bounds[3])] # y data cut down to bkg boundaries
    
    
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
    
    BasePara, residuals, rank, singular_values, rcond = polyfit(bkg_x,bkg_signal,1, full = True)
    m = BasePara[0]
    b = BasePara[1]
    baseline = polyval([m,b],x_fit)
    #print 'Baseline: ', '%.5f' %m, 'x + ', '%.5f' %b
    #print 'R2 Baseline: ', 1-(residuals[0]/sum(bkg_signal))
    
    ### Option for 3rd Order polynomial fit for baseline fit
    #BasePara, residuals, rank, singular_values, rcond  =  polyfit(bkg_wave,bkg_signal,3, full = True)
    #a = BasePara[0]
    #b = BasePara[1]
    #c = BasePara[2]
    #d = BasePara[3]
    #baseline = polyval([a,b,c,d],x_fit)
    ##print 'Baseline: ', a, 'x**3 + ', b, 'x**2 + ', c , 'x + ', d
    ##print 'R2 Baseline: ', 1-(residuals[0]/sum(bkg_signal))
    
    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    ax.plot(x_fit , signal_fit,'.k')
    #ax.plot(wave , signal,'--k')
    ax.plot(x_fit, baseline, '-r')
    ax.set_xlabel('Raman Shift (cm-1)')
    ax.set_ylabel('Raman Intensity')
    ax.set_ylim(0, 1.5*max(signal_fit))
    #plt.show()
    fname = str(SaveName) + '_base_2LS.jpg'
    #plt.savefig(fname)
    plt.close()
      
    # Baseline Correction
    signal_fit = signal_fit - baseline #corrects for baseline 
    Signal_hi = signal[np.where(abs(wave - 1910) ==min(abs(wave - 1910)))[0][0]]
    #Signal_2500 = signal[np.where(abs(wave - 2500) ==min(abs(wave - 2500)))[0][0]]
    #Signal_800 = signal[np.where(abs(wave - 800) ==min(abs(wave - 800)))[0][0]]
    Signal_lo = signal[np.where(abs(wave - 925) ==min(abs(wave - 925)))[0][0]]
    #print 'Baseline Flatness:', Signal_hi/Signal_lo #NOTE if there's noise right here this value will be not so hot
    
    
    def lorentz(xc, w, I):  
        global x_fit
        s = ((x_fit - xc)/w)
        return (I)*(1/(1+s**2)) #Wikipedia definition.  gives correct value for intensity.  using this definition means you don't need to account for the peak width in the intensity.
    
    def BWF(xc,w,I):
        global qBWF,x_fit
        q = qBWF
        s = ((x_fit - xc)/w)
        return (I)*((1+s/qBWF)**2/(1+s**2))

    def EnterData():
        
        global FitParam ,G_ints, D_ints, NumParams, G_bounds, D_bounds, bounds
       
        FitParam[0] = G_bounds[0] # G peak position
        FitParam[1]  = D_bounds[0] # D peak position
        FitParam[2]  = G_bounds[2]  # G peak width
        FitParam[3]  = D_bounds[2]  # D peak width
        
        '''need to work G and D starting peak intensities from initial values'''
        FitParam[4]  = IIM*G_ints  # G peak intensity  
        FitParam[5]  = IIM*D_ints  # D peak intensity
        
        bounds[0] = (G_bounds[0]-G_bounds[1],G_bounds[0]+G_bounds[1])
        bounds[1] = (D_bounds[0]-D_bounds[1],D_bounds[0]+D_bounds[1])
        
        bounds[2] = (G_bounds[2]-G_bounds[3],G_bounds[2]+G_bounds[3])
        bounds[3] = (D_bounds[2]-D_bounds[3],D_bounds[2]+D_bounds[3])
        
        bounds[4] = (0.5*G_ints,1.2*G_ints)
        bounds[5] = (0.5*D_ints,1.2*D_ints)
        #print 'Bounds are ',bounds
        
    def Evaluate(EvalSimp):   
        global NumParams, G_bounds, D_bounds
        global x_fit, signal_fit, Residuals
         
        #This is the function that gives a goodness of fit
        
        '''
        Need to 
        (1) evaluate Lorentzian for each peak
        (2) Add all peak fits together for total peak fit
        (3) Subtract total peak fit from real data for initial residuals
        (4) punish residuals if various parameters are too far from real
        '''
        
        Gfit = BWF(EvalSimp[0],EvalSimp[2],EvalSimp[4])
        Dfit = lorentz(EvalSimp[1],EvalSimp[3],EvalSimp[5])
         
        EvalFit = Gfit + Dfit
         
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
        #CheckThese = [1290,1300,1310,1320,1330,1340,1350,1360,1470,1490,1520,1550,1560,1570,1580,1590,1600]
        #for n in CheckThese:
        #    if Residuals[n-925] > 0.1*signal_fit[n-925]: #xvalue of resds is wavenumber - 925
        #        ErrorWeight = ErrorWeight*(100*Residuals[n-925])/signal_fit[n-925]
        #        print n, 'is ',Residuals[n-925]/signal_fit[n-925]
        
        for n in range(0,len(EvalSimp)):
            if EvalSimp[n] <0:
                ErrorWeight = ErrorWeight *10
                #print 'cond 1'
            if n == 0:
                if abs(1580 - EvalSimp[0])>G_bounds[1]:
                    ErrorWeight = ErrorWeight *10
                    #print 'cond 3'
            if n == 1:
                if abs(1350 - EvalSimp[1])>D_bounds[1]:
                    ErrorWeight = ErrorWeight *10
                    #print 'cond 4'
        
        #print 'ErrorWeight for array ', NumIter,' is ',ErrorWeight
        WeightResiduals = ErrorSum*ErrorWeight
        
        return(WeightResiduals)
    
    
    G_ints = max(signal_fit[(G_bounds[0]-G_bounds[1] <= x_fit) & (x_fit <= G_bounds[0]+G_bounds[1])])
    D_ints = max(signal_fit[(D_bounds[0]-D_bounds[1] <= x_fit) & (x_fit <= D_bounds[0]+D_bounds[1])]) 
    #print 'initial D/G intensity ratio: ',D_ints/G_ints
 
      
    
    EnterData()
    #print 'FitParam ', FitParam
    minres = minimize(Evaluate,FitParam,method='SLSQP',tol=1e-8,bounds=bounds,options={'maxiter':1e5}) 
    if minres.success == 1:
        Mean=minres.x
        NumIter = minres.nit
    else:
        print minres.message
        Mean = np.zeros(NumParams)
        continue
    Gfit = lorentz(Mean[0],Mean[2],Mean[4])
    Dfit = lorentz(Mean[1],Mean[3],Mean[5])
        
    ModelFit = Gfit + Dfit 
      
    '''
    fit_error = abs(Residuals)/np.sum(Model_fit)
    '''
    G_ints = Mean[4]
    D_ints = Mean[5]
    
    # Figure 4 Plot of individual peak fits and total peak fit with experimental
    fig = plt.figure(4)
    ax40 = fig.add_subplot(111)
    ax40.plot(x_fit, signal_fit,'.k', label = 'Experimental')
    ax40.plot(x_fit, Gfit,'-g', label = 'G Peak Fit')
    ax40.plot(x_fit, Dfit,'-b', label = 'D Peak Fit')
    ax40.plot(x_fit, ModelFit,'-r', label = 'Summed Peak Fit')
    ax40.set_xlabel(r'Raman Shift / cm$^-$$^1$', fontsize=18)
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
    ax40 = fig.add_subplot(111)
    ax41 = plt.axes([0.185, .88, .787, .1])
    ax41.plot(x_fit,Residuals, '.b')
    #ax41.set_ylabel('Residuals')
    plt.setp(ax41, xticks=[1000,1200,1400,1600,1800], yticks=[-100,0,100])
    ax41.tick_params(direction='in',labelbottom='off',labelleft='off')
    plt.autoscale(enable=True, axis='x', tight=True) 
    
    #plt.ylim(min(Residuals)*1.15,max(Residuals)*1.15)
    plt.ylim(min(Residuals)*1.15,130)
    
    plt.savefig(SaveName + 'fit_2BWFLS.jpg')
    #plt.show()
    plt.close()
    '''        
    # Uncertainty analysis
    G_error = [0]*iterations
    D_error = [0]*iterations
    for i in range(0, iterations):
        model_error = Model_fit + fit_error*randn(1)/4*Model_fit
        G_error[i] =model_error[np.where(x_fit == around(fit_params[1],0) )[0][0]]
        D_error[i] = model_error[np.where(x_fit == around(fit_params[4],0) )[0][0]]
    
    G_ints_error =  np.std(G_error) 
    D_ints_error =  np.std(D_error)
        
    #print "R2 Model Fit: ", 1-fit_error
    '''
    #print 'G Band Peak Position: ', Mean[0]
    #print 'G Band Peak Width: ', Mean[2]
    #print 'G Band Peak Intensity: ', '%.3f' %G_ints#, '+/- ',  '%.3f' %G_ints_error
    #print 'D Band Peak Position: ', Mean[1]
    #print 'D Band Peak Width: ', Mean[3]
    #print 'D Band Peak Intensity: ', '%.3f' %D_ints  # , '+/- ',  '%.3f' %D_ints_error 
    
    # ID/IG Ratio
    Exp_ratio = D_ints/G_ints
    '''
    Exp_ratio_err = Exp_ratio*((G_ints_error/G_ints)**2 + (D_ints_error/D_ints)**2)**0.5
    '''
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
    '''
    low_La_min_ratio = La_low[np.where(abs(Model_Ratio_low - (Exp_ratio-Exp_ratio_err) ) == min(abs(Model_Ratio_low - (Exp_ratio-Exp_ratio_err))))[0][0]]
    low_La_max_ratio = La_low[np.where(abs(Model_Ratio_low - (Exp_ratio+Exp_ratio_err)) == min(abs(Model_Ratio_low - (Exp_ratio+Exp_ratio_err))))[0][0]]
    low_La_error_ratio = abs(low_La_min_ratio - low_La_max_ratio)
    low_La_error = low_La*((low_La_error_ratio/low_La)**2 + (low_La_error_Ca/low_La)**2)**0.5
    #print 'Conjugation Length (low): ', '%.3f' %low_La,  '+/- ', '%.3f' %low_La_error, 'nm'
    '''
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
    '''
    high_La_min_ratio = La_high[np.where(abs(Model_Ratio_high - (Exp_ratio-Exp_ratio_err) ) == min(abs(Model_Ratio_high - (Exp_ratio-Exp_ratio_err))))[0][0]]
    high_La_max_ratio = La_high[np.where(abs(Model_Ratio_high - (Exp_ratio+Exp_ratio_err)) == min(abs(Model_Ratio_high - (Exp_ratio+Exp_ratio_err))))[0][0]]
    high_La_error_ratio = abs(high_La_min_ratio - high_La_max_ratio)
    high_La_error = high_La*((high_La_error_ratio/high_La)**2 + (high_La_error_Ca/high_La)**2)**0.5
    #print 'Conjugation Length (high): ', '%.3f' %high_La,  '+/- ', '%.3f' %high_La_error, 'nm'
    '''
    Ratio = [Exp_ratio, Exp_ratio]
    La = [low_La, high_La]
    '''   '''
    fig = plt.figure(3)
    ax = fig.add_subplot(111)
    ax.loglog(La_Range , Model_Ratio,'-k')
    ax.fill_between(La_Range, Min_Ratio, Max_Ratio,  facecolor='gray', alpha = 0.25)
    ax.loglog(La , Ratio,'or')
    high_label = '%.3f' %high_La #+ '+/- ' + '%.3f' %high_La_error + 'nm'
    low_label = '%.3f' %low_La #+ '+/- ' + '%.3f' %low_La_error + 'nm'
    ax.text(high_La, Exp_ratio, high_label, va="top", ha = "center", size = 10)
    ax.text(low_La, Exp_ratio, low_label, va="top", ha = "center", size = 10)
    ax.set_xlabel('La')
    ax.set_ylabel('ID/IG')
    ax.set_xlim(0.1, 100)
    ax.set_ylim(0.01, 100)
    plt.savefig(SaveName + '_Ratio_2BWFLS.jpg')
    
    #plt.show()
    plt.close()
    '''   '''
    
    f = open(SaveName+'_2BWFLS_fitfile.txt',"w")
    f.write("{}\t{}\t{}\n".format('Collection Details',CollDet,''))
    
    f.write("{}\t{}\t{}\n".format('Baseline R2', 1-(residuals[0]/sum(bkg_signal)), 0) )
    f.write("{}\t{}\t{}\n".format('Baseline Flatness', Signal_hi/Signal_lo, 0) )
    #f.write("{}\t{}\t{}\n".format('Model Fit R2', 1-fit_error, 0) )
    f.write("{}\t{}\t{}\n".format('Q factor', qBWF, 0) )
    f.write("{}\t{}\t{}\n".format('G Band Peak Position', Mean[0], 0) )
    f.write("{}\t{}\t{}\n".format('G Band Peak Width', Mean[2], 0 ))
    f.write("{}\t{}\t{}\n".format('G Band Peak Intensity', G_ints,0))
    f.write("{}\t{}\t{}\n".format('D Band Peak Position',Mean[1],0))
    f.write("{}\t{}\t{}\n".format('D Band Peak Width',Mean[3],0))
    f.write("{}\t{}\t{}\n".format('D Band Peak Intensity',D_ints,0))
    f.write("{}\t{}\t{}\n".format('Conjugation Length (low)',low_La,0))
    f.write("{}\t{}\t{}\n".format('Conjugation Length (high)',high_La,0))
    f.write("{}\t{}\t{}\n".format('ID/IG',Exp_ratio,0))
    f.write("{}\t{}\t{}\n".format('Location', Location,''))
    f.close()

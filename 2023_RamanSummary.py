# 202309 automatically produces a correlation plot colored according to average laser color
#        automatically drops any fit with a peak fit R2 less than 0.7 (can be set at top)

import os
import fnmatch
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
plt.style.use('default')
plt.rcParams.update({'font.family': 'Times New Roman'})
plt.rcParams.update({'font.size': 20})

DropR2 = 0.7
tempData = []


for file in os.listdir('.'):
    if fnmatch.fnmatch(file, '*_fitfile.txt'):
        fullFilename = file
        filename = file[:-12]
        print(filename)
        flen = len(filename)-9
        
        data = np.genfromtxt(fullFilename, delimiter = '\t', usecols = (1,), skip_header = (2))
        unc= np.genfromtxt(fullFilename, delimiter = '\t', usecols = (2,), skip_header = (2))

        position = data[0]
        location = data[1]        
        exc_laser = data[2]
        num_pks = data[3]
        fit_version = unc[3]
        
        baseline_order = data[4]
        baseline_r2 = data[5]
        baseline_flatness = data[6]
        
        pkfit_r2 = data[7]
        pkfit_see = data[8]
        
        if pkfit_r2 < DropR2:
            continue

        
        q_BWF = data[9]
        avgCorr = unc[9]
        g_pos = data[10]
        g_pos_unc = unc[10]
        g_width = data[11]
        g_width_unc = unc[11]
        g_intensity = data[12]
        g_ints_unc = unc[12]
        
        d_pos = data[13]
        d_pos_unc = unc[13]
        d_width = data[14]
        d_width_unc = unc[14]
        d_intensity = data[15]
        d_ints_unc = unc[15]
        
        d2_pos = data[16]
        d2_pos_unc = unc[16]
        d2_width = data[17]
        d2_width_unc = unc[17]
        d2_intensity = data[18]
        d2_ints_unc = unc[18]
        
        d3_pos = data[19]
        d3_pos_unc = unc[19]
        d3_width = data[20]
        d3_width_unc = unc[20]
        d3_intensity = data[21]
        d3_ints_unc = unc[21]
        
        d4_pos = data[22]
        d4_pos_unc = unc[22]
        d4_width = data[23]
        d4_width_unc = unc[23]
        d4_intensity = data[24]
        d4_ints_unc = unc[24]

        tpa_pos = data[25]
        tpa_pos_unc = unc[25]
        tpa_width = data[26]
        tpa_width_unc = unc[26]
        tpa_intensity = data[27]
        tpa_ints_unc = unc[27]

        Lalow_val = data[28]
        Lalow_unc = unc[28]
        ID_IG = data[30]
        ID_IG_unc = unc[30]       
        
        ID_IT = data[31]
        ID_IT_unc = unc[31] 
        ID2_IT = data[32]
        ID2_IT_unc = unc[32]        
        ID3_IT = data[33]
        ID3_IT_unc = unc[33]        
        ID4_IT = data[34]
        ID4_IT_unc = unc[34]        ### this is where TPA stuff will go
        IG_IT = data[35]
        IG_IT_unc = unc[35]        
        
        ID2_IG = data[36]
        ID2_IG_unc = unc[36]        
        ID3_ID = data[37]
        ID3_ID_unc = unc[37]        
        ID4_ID = data[38]
        ID4_ID_unc = unc[38]      ### more TPA stuff will go here
        
        lowbkg_st = data[39]
        lowbkg_end = unc[39]        
        hibkg_st = data[38]
        hibkg_end = unc[38]       
           
        scan_info = np.genfromtxt(fullFilename, dtype = str, delimiter = '\t', usecols = (1), skip_footer = (len(data)+1))
        fitting_info = filename.split(sep='_')[-1]
        
        #########################  
        # to solve concat/append unc. here we concat all the bits into individual lists, then after all loops run, concat into a df below

        
        tempData.append({'Position': position, 
            'Location': location,
            'Scan Info': scan_info, 
            'Exc Laser': exc_laser, #in nm
            'Exc Energy': 299792458 * 6.6261e-34 /exc_laser/1e-9/1.602176565E-19, # in eV
            'Num Peaks Fit': num_pks,
            'Fitting routine': fitting_info + 'v.' +str(fit_version),
            
            'Baseline Order': baseline_order,
            'Baseline Adj R2': baseline_r2,
            'Baseline Flatness': baseline_flatness, 
            'low bkg': str(lowbkg_st) + '-' + str(lowbkg_end),
            'hi_bkg': str(hibkg_st) + '-' + str(hibkg_end),
            
            'PkFit Adj R2': pkfit_r2, 
            'PkFit SEE %D': pkfit_see/d_intensity, 
                        
            'qBWF': q_BWF,
            'average Corr': avgCorr,
            
            'G Peak Position': g_pos, 
            'G PeakPos Unc':g_pos_unc,
            'G Peak Width':g_width, 
            'G PeakWid Unc':g_width_unc,
            'G Intensity': g_intensity ,
            'G Intensity unc': g_ints_unc,
            
            'D Peak Position': d_pos,
            'D PeakPos Unc': d_pos_unc,
            'D Peak Width': d_width, 
            'D PeakWid Unc': d_width_unc,
            'D Intensity': d_intensity,
            'D Intensity unc': d_ints_unc,
            
            'D2 Peak Position': d2_pos,
            'D2 PeakPos Unc': d2_pos_unc,
            'D2 Peak Width': d2_width, 
            'D2 PeakWid Unc': d2_width_unc,
            'D2 Intensity': d2_intensity,
            'D2 Intensity unc': d2_ints_unc,
            
            'D3 Peak Position': d3_pos,
            'D3 PeakPos Unc': d3_pos_unc,
            'D3 Peak Width': d3_width, 
            'D3 PeakWid Unc': d3_width_unc,
            'D3 Intensity': d3_intensity,
            'D3 Intensity unc': d3_ints_unc,
            
            'D4 Peak Position': d4_pos,
            'D4 PeakPos Unc': d4_pos_unc,
            'D4 Peak Width': d4_width, 
            'D4 PeakWid Unc': d4_width_unc,
            'D4 Intensity': d4_intensity,
            'D4 Intensity unc': d4_ints_unc,
            
            'TPA Peak Position': tpa_pos,
            'TPA PeakPos Unc': tpa_pos_unc,
            'TPA Peak Width': tpa_width, 
            'TPA PeakWid Unc': tpa_width_unc,
            'TPA Intensity': tpa_intensity,
            'TPA Intensity unc': tpa_ints_unc,
            
            'Conj Length (low)': Lalow_val , 
            'Conj Length unc': Lalow_unc , 
            'ID IG Ratio': ID_IG ,
            'ID IG Ratio unc': ID_IG_unc ,
            
            'ID IT Ratio': ID_IT ,
            'ID IT Ratio unc': ID_IT_unc ,
            'ID2 IT Ratio': ID2_IT ,
            'ID2 IT Ratio unc': ID2_IT_unc ,
            'ID3 IT Ratio': ID3_IT ,
            'ID3 IT Ratio unc': ID3_IT_unc ,
            'ID4 IT Ratio': ID4_IT ,
            'ID4 IT Ratio unc': ID4_IT_unc ,
            #'ITPA IT Ratio': ITPA_IT ,
            #'ITPA IT Ratio unc': ITPA_IT_unc ,      
            'IG IT Ratio': IG_IT ,
            'IG IT Ratio unc': IG_IT_unc ,            
            
            'ID2 IG Ratio': ID2_IG ,
            'ID2 IG Ratio unc': ID2_IG_unc ,
            'ID3 ID Ratio': ID3_ID ,
            'ID3 ID Ratio unc': ID3_ID_unc ,
            'ID4 ID Ratio': ID4_ID ,
            'ID4 ID Ratio unc': ID4_ID_unc ,
            #'ITPA ID Ratio': ITPA_ID ,
            #'ITPA ID Ratio unc': ITPA_ID_unc ,
            
            'Filename': filename ,
            #'Noise': noise ,
            #'SNR D band': SNRD ,  
            })
        
# this is where concat the individual lists into dataframes is--outside the loop  
Data = pd.DataFrame(tempData)          
Data.to_csv((fitting_info + 'RamanFit_summary.csv'),index=False,header=True)  

# Let's create a correlation plot
if num_pks.max() == 2:
    IndVar = Data[['D Peak Position','D Peak Width','G Peak Position','G Peak Width','ID IG Ratio']]
elif num_pks.max() == 3:
    IndVar = Data[['D Peak Position','D Peak Width','G Peak Position','G Peak Width',
                   'ID IG Ratio', 'D3 Peak Position','D3 Peak Width', 'ID3 ID Ratio']]
elif num_pks.max() == 4:
    #print('found 4 loop')
    IndVar = Data[['D Peak Position','D Peak Width','G Peak Position','G Peak Width',
                   'ID IG Ratio', 'D3 Peak Position','D3 Peak Width',  
                   'D4 Peak Position','D4 Peak Width']]
elif num_pks.max() == 5:
    IndVar = Data[['D Peak Position','D Peak Width','G Peak Position','G Peak Width',
                   'ID IG Ratio', 'D3 Peak Position','D3 Peak Width', 
                   'D4 Peak Position','D4 Peak Width', 'D2 Peak Position',
                   'D2 Peak Width']]
else:
    IndVar = Data[['D Peak Position','D Peak Width','G Peak Position','G Peak Width',
                   'ID IG Ratio', 'D3 Peak Position','D3 Peak Width', 
                   'D4 Peak Position','D4 Peak Width', 'D2 Peak Position',
                   'D2 Peak Width', 'TPA Peak Position', 'TPA Peak Width']]
if exc_laser.mean() < 500:
    lasercolor = 'blue'
elif exc_laser.mean() <600:
    lasercolor = 'green'
elif exc_laser.mean() < 650:
    lasercolor = 'orange'
elif exc_laser.mean() < 800:
    lasercolor = 'red'
else:
    lasercolor = 'black'

plottitle = fitting_info + ' ' + str(round(exc_laser.mean())) + 'nm'
savename = 'PRUNEDcorr_' + fitting_info + '_' + str(round(exc_laser.mean())) + '.jpg'

g = sns.PairGrid(IndVar, diag_sharey=False)
print('step 1')
g.map_lower(sns.scatterplot, alpha=0.5, color = lasercolor)
print('step 2')
g.map_upper(sns.kdeplot, alpha=0.5, color = lasercolor)
print('step 3')
g.map_diag(sns.kdeplot, color = lasercolor)
print('step 4')
g.fig.suptitle(plottitle)
g.fig.subplots_adjust(top=.95)
g.savefig(savename)

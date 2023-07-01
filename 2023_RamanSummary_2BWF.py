import os
import fnmatch
import pandas as pd
import numpy as np


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
        
        q_BWF = data[9]
        
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
        
        Lalow_val = data[25]
        Lalow_unc = unc[25]
        ID_IG = data[27]
        ID_IG_unc = unc[27]
        
        scan_info = np.genfromtxt(fullFilename, dtype = str, delimiter = '\t', usecols = (1), skip_footer = (len(data)+1))
        fitting_info = filename.split(sep='_')[-1]
        
        #########################  
        # to solve concat/append unc.  concat all the bits into individual lists, then after all loops run, concat into a df below

        
        tempData.append({'Position': position, 
            'Location': location,
            'Scan Info': scan_info, 
            'Exc Laser': exc_laser,
            'Num Peaks Fit': num_pks,
            'Fitting routine': fitting_info + 'v.' +str(fit_version),
            
            'Baseline Order': baseline_order,
            'Baseline Flatness': baseline_flatness, 
            
            'PkFit R2': pkfit_r2, 
            'PkFit SEE': pkfit_see, 
            
            'qBWF': q_BWF, 
            
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
            
            'Conj Length (low)': Lalow_val , 
            #'Conj Length unc': La_unc , 
            'ID IG Ratio': ID_IG ,
            'ID IG Ratio unc': ID_IG_unc ,
            'Filename': filename ,
            #'Noise': noise ,
            #'SNR D band': SNRD ,  
            })
        
# this is where concat the individual lists into dataframes will go--outside the loop  
Data = pd.DataFrame(tempData)          
Data.to_csv((fitting_info + 'RamanFit_summary.csv'),index=False,header=True)  



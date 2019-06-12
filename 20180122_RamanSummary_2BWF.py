import os
import fnmatch
from pandas import *
import pandas as pd
import matplotlib.pyplot as plt


Data = pd.DataFrame()
for file in os.listdir('.'):
    if fnmatch.fnmatch(file, '*_fitfile.txt'):
        fullFilename = file
        filename = file[:-12]
        print filename
        flen = len(filename)-9
        position = int(filename[flen:flen+2])
        
        data = np.genfromtxt(fullFilename, delimiter = '\t', usecols = (1,), skip_header = (1))
        error= np.genfromtxt(fullFilename, delimiter = '\t', usecols = (2,), skip_header = (1))
        baseline_r2 = data[0]
        baseline_flatness = data[1]
        #model_r2 = data[2]
        qBWF = data[2]
        g_pos = data[3]
        g_width = data[4]
        g_intensity = data[5]
        g_error = error[5]
        d_pos = data[6]
        d_width = data[7]
        d_intensity = data[8]
        d_ints_error = error[8]
        La = data[9]
        La_error = error[9]
        ID_IG = data[11]
        ID_IG_error = error[11]
        Location = data[12]
        #SNRD = data[11]
        #noise = data[12]
        
        ScanInfo = np.genfromtxt(fullFilename, dtype = str, delimiter = '\t', usecols = (1), skip_footer = (len(data)))
        
        Data = Data.append({'Position': position, 
            'Location': Location,
            'Scan Info': ScanInfo, 
            'Baseline Flatness': baseline_flatness, 
            #'Model R2': model_r2, 
            'qBWF': qBWF, 
            'G Peak Position': g_pos, 
            'G Peak Width':g_width, 
            'G Intensity': g_intensity ,
            'G Intensity error': g_error  ,
            'D Peak Position': d_pos,
            'D Peak Width': d_width, 
            'D Intensity': d_intensity,
            'D Intensity error': d_ints_error,
            'Conj Length': La , 
            'Conj Length error': La_error , 
            'ID IG Ratio': ID_IG ,
            'ID IG Ratio error': ID_IG_error ,
            'Filename': filename ,
            #'Noise': noise ,
            #'SNR D band': SNRD ,  
            }, ignore_index=True)
Data.to_csv(filename + '_summary.csv',index=True,header=True)  


        


fig = plt.figure(1)
ax = fig.add_subplot(111)
#ax.plot(Data['Position'].values , Data['SNR D band'].values ,'ob')
ax.plot(Data['Position'].values , Data['Baseline Flatness'].values ,'ob')
#ax.set_ylim(0, 1.2*max(abs_coeff))
#ax.set_xlim(0, max(hv))
ax.set_xlabel('Position')
ax.set_ylabel('Baseline Flatness')
plt.savefig(filename + '_summary.jpg')
plt.show()
plt.close()

fig = plt.figure(2)
ax = fig.add_subplot(111)
ax.plot(Data['Position'].values , Data['ID IG Ratio'].values ,'og')
#ax.plot(Data['Position'].values , Data['SNR D band'].values ,'ob')
#ax.plot(Data['Position'].values , Data['Baseline Flatness'].values ,'ob')
#ax.set_ylim(0, 1.2*max(abs_coeff))
#ax.set_xlim(0, max(hv))
ax.set_xlabel('Position')
ax.set_ylabel('ID/IG ratio')
plt.savefig(filename + '_summary2.jpg')
plt.show()
plt.close()
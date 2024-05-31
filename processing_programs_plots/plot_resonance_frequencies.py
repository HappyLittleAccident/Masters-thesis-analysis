# -*- coding: utf-8 -*-
"""
Created on Wed May 22 14:01:57 2024

@author: Marek
"""

import numpy as np
import matplotlib.pyplot as plt
from glob import glob

plt.close('all')

temperatures_ampsweeps = [
 'T1325',
 'T1350_4',
 'T1375',
 'T1400',
 'T1425',
 'T1450_2',
 'T1475',
 'T1500',
 'T1550',
 'T1650_2',
 'T1700',
 'T1750',
 'T1850',
 'T1950']
temperatures = [
 'T1325',
 'T1350',
 'T1375',
 'T1400',
 'T1425',
 'T1450',
 'T1475',
 'T1500',
 'T1550',
 'T1650',
 'T1700',
 'T1750',
 'T1850',
 'T1950']
resonators = ['A','B','C','D']


folder_params = r'D:\OneDrive_copy\OneDrive - Univerzita Karlova\DATA\2023_4_Helmholtz_resonators_recal_2\resonance_freq'

folder_ampsweeps_raw = r'D:\Github\Masters-thesis-analysis\filip_calibration\Helmholtz_res_drive_dep'

folder_vld = r'D:\OneDrive_copy\OneDrive - Univerzita Karlova\DATA\2023_4_Helmholtz_resonators_recal_2\Helmholtz_res_VLD'

folder_fsweeps = r'D:\OneDrive_copy\OneDrive - Univerzita Karlova\DATA\2023_4_Helmholtz_resonators_recal_2\Helmholtz_fsweeps'


for letter in resonators:
    fig,ax = plt.subplots(2,1,sharex=True)
    d_par = np.load( folder_params + fr'\500{letter}.npy',allow_pickle=True).item()
    
    f0_a = []
    f0_f = []
    f0s_pars = d_par['frequency (Hz)']
    f0_pars = [np.mean(x) for x in f0s_pars]
    w_vld = []

    
    ws_pars = d_par['width (Hz)']
    w_pars = [np.mean(x)/2 for x in ws_pars]
    
    temps = d_par['temperature (K)']
    print('letter: ',letter)
    
    for temp_a,temp in zip(temperatures_ampsweeps,temperatures):
        files_ampsweeps = glob(folder_ampsweeps_raw+rf'\{temp_a}\ampsweeps_fast\resonator_500{letter}_updowndown\*.npy')
        files_vld = glob(folder_vld +rf'\{temp}\500{letter}\*.npy')
        files_fsweeps = glob(folder_fsweeps +rf'\{temp}\500{letter}\*.npy')
        
        d_a = np.load(files_ampsweeps[0],allow_pickle = True).item()
        d_vld = np.load(files_vld[0],allow_pickle = True).item()
        d_f = np.load(files_fsweeps[20],allow_pickle = True).item()
        
        f0_a.append(d_a['freq (Hz)'])
        
        w_vld.append(d_vld['Width (Hz)'])
        
        r = d_f['Velocity (m/s)']
        f = d_f['Freq (Hz)']
        f0_f.append(f[np.argmax(r)])
        if temp == 'T1850':
            plt.figure()
            plt.plot(f,r)
        
        
    ax[0].plot(temps,f0_a,'ks-',label = 'ampsweeps')
    ax[0].plot(temps,f0_f,'bo-',label = 'fsweeps')
    ax[0].plot(temps,f0_pars,'ro-',label = 'pars file')
        
    ax[1].plot(temps,w_vld,'b',label = 'vld')
    ax[1].plot(temps,w_pars,'r',label = 'pars file')    
        
        
    ax[0].legend()
    ax[1].legend()
        
    ax[1].set_xlabel('Temp (K)')
    ax[0].set_ylabel('Res freq (Hz)')
    ax[1].set_ylabel('Width (Hz)')
    while not plt.waitforbuttonpress():
        pass



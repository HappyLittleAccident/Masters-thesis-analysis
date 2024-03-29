# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 16:55:10 2023

@author: Marek
"""

import numpy as np
import matplotlib.pyplot as plt
import os
from glob import glob

plt.close('all')
folder = r'D:\Github\Masters-thesis-analysis\filip_calibration\2023_4_Helmholtz_resonators_recal_2\Helmholtz_fsweeps'
folder_ampsweeps = r'D:\OneDrive_copy\OneDrive - Univerzita Karlova\DATA\2023_4_Helmholtz_resonators_recal_2\Helmholtz_res_drive_dep'
folder_fsweep_amps = r'D:\OneDrive_copy\OneDrive - Univerzita Karlova\DATA\2023_4_Helmholtz_resonators_recal_2\Helmholtz_res_drive_dep_fsweeps'

temps = ['T1325', 'T1350', 'T1375', 'T1400', 'T1425', 'T1450','T1475', 'T1500', 'T1550', 'T1650', 'T1700', 'T1750', 'T1850', 'T1950']
resonators = ['500A','500B','500C','500D']

cmap = plt.get_cmap('viridis')

fig,ax = plt.subplots(1,1,sharex=True)
for resonator in resonators[0:3]:
    for temp in temps:
        ax.clear()
        filenames = glob(os.path.join(folder,temp,f'{resonator}','*.npy'))
        filenames_ampsweeps = glob(os.path.join(folder_ampsweeps,temp,fr'ampsweeps_fast\resonator_{resonator}_updowndown','*.npy'))
        filename_famps = os.path.join(folder_fsweep_amps,fr'{resonator}',f'{temp}.npy')
        
        # A_max = np.max([np.load(file_temp,allow_pickle=True).item()['App (V)'] for file_temp in filenames])
        for file in filenames:
            data = np.load(file,allow_pickle=True).item()
            # f = data['Freq (Hz)']
            r = np.max(data['Velocity (m/s)'])
            A = data['App_gen (V)']
            p = data['Pressure (Pa/m)']
            ax.plot(p,r,'.',c = 'k',alpha=0.5)
        
        # file_ampsweep  = filenames_ampsweeps[0]  
        # data_a = np.load(file_ampsweep,allow_pickle=True).item()
        # # f = data['Freq (Hz)']
        # r_a = data_a['upVelocity (m/s)']
        # p_a = data_a['upPressure (Pa/m)']
        
        # axx.plot(p_a,r_a,'.',c = 'b',alpha=0.5)
          
        data_f = np.load(filename_famps,allow_pickle=True).item()
        # f = data['Freq (Hz)']
        r_f = data_f['Velocity (m/s)']
        p_f = data_f['Pressure grad (Pa/m)']
        
        ax.plot(p_f,r_f,'.',c = 'r',alpha=0.5)
        ax.set_title(f'{resonator} {temp}')
        fig.canvas.draw()
        while not plt.waitforbuttonpress():
            pass
        
            
            
    
    

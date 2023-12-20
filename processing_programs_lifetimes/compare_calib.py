# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 16:09:48 2023

@author: Marek
"""

import matplotlib.pyplot as plt
from glob import glob
import numpy as np

folder_old = r"D:\OneDrive_copy\OneDrive - Univerzita Karlova\DATA\2023_4_Helmholtz_resonators_recal\Helmholtz_buffer_measurements"
folder_new = r"D:\OneDrive_copy\OneDrive - Univerzita Karlova\DATA\2023_4_Helmholtz_resonators_recal_2\Helmholtz_buffer_measurements"


temps = ['T1350','T1450','T1650','T1850']
plt.close('all')

for temp in temps:
    files_new = glob(folder_new + fr'\{temp}\buffer\resonator_500B\*.npy')
    files_old = glob(folder_old + fr'\{temp}\buffer\resonator_500B\*.npy')
    
    fig,ax = plt.subplots()
    for fn,fo in zip(files_new,files_old):
        dn = np.load(fn,allow_pickle=True).item()
        do = np.load(fo,allow_pickle=True).item()
        
        vn = dn['Velocity (m/s)']
        vo = do['Velocity (m/s)']
        
        pn = dn['Pressure (Pa/m)']
        po = do['Pressure (Pa/m)']
        
        if (temp == 'T1850' and po < 6000) or np.mean(vo)>500:
            continue
        
        ax.plot(pn,np.mean(vn),'r.')
        ax.plot(po,np.mean(vo),'k.')
        
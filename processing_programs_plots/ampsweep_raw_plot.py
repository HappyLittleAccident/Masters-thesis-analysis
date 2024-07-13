# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 16:55:10 2023

@author: Marek
"""

import numpy as np
import matplotlib.pyplot as plt
import os
from glob import glob
from format_figure import use_latex,polish
import matplotlib as mpl

plt.close('all')
folder = r'D:\Github\Masters-thesis-analysis\filip_calibration\Helmholtz_res_drive_dep'
temps = ['T1325', 'T1350_4', 'T1375', 'T1400', 'T1425', 'T1450_2', 'T1475',
         'T1500', 'T1550', 'T1650_2', 'T1700', 'T1750', 'T1850', 'T1950']
resonators = ['500A', '500B', '500C']
names = ['C', 'J', 'R']
folder_parameters = r'D:\OneDrive_copy\OneDrive - Univerzita Karlova\DATA\2023_4_Helmholtz_resonators_recal_2\resonance_freq'


# use_latex()
plt.rcdefaults()

cmap = plt.get_cmap('winter')
viridis = plt.get_cmap('viridis')

# plt.ioff()
plt.ion()

for i,resonator in enumerate(resonators[0:3]):
    # fig,ax = plt.subplots(3,1,sharex=False,dpi=250,gridspec_kw={'height_ratios': [1, 1,0.05]})
    d_pars = np.load(folder_parameters+'\\'+resonator+'.npy',allow_pickle=True).item()
    gamma = d_pars['width (Hz)']
    print(f'\nresonator {resonator}\n')
    for j,temp in enumerate(temps[0::]):
        gamma_temp = np.mean(gamma[j])/2
        filenames = glob(os.path.join(folder,f'{temp}','ampsweeps_fast',f'resonator_{resonator}_updowndown','*.npy'))
        # A_max = np.max([np.load(file_temp,allow_pickle=True).item()['App (V)'] for file_temp in filenames])
        for file in filenames[5:6]:
            data = np.load(file,allow_pickle=True).item()
            
            
            x = data['x (A)']
            y = data['y (A)']
            print(gamma_temp)
            print(data['time constant (s)'],data['wait time (s)'],f'{2/gamma_temp:.2f}')
            
                
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 16:55:10 2023

@author: Marek
"""

import numpy as np
import matplotlib.pyplot as plt
import os
from glob import glob
from format_figure import polish,use_latex

plt.close('all')
folder = r'..\filip_calibration\2023_4_Helmholtz_resonators_recal_2\Helmholtz_fsweeps'

temp = 'T1350'
resonators = ['500A','500B','500C']

cmap = plt.get_cmap('viridis')
use_latex()

for resonator in resonators[:]:
    filenames = glob(os.path.join(folder,temp,f'{resonator}','*.npy'))
    fig,ax = plt.subplots(1,1,sharex=True)
    A_max = np.max([np.load(file_temp,allow_pickle=True).item()['App_gen (V)'] for file_temp in filenames])
    for file in filenames[0:]:
        data = np.load(file,allow_pickle=True).item()
        f = data['Freq (Hz)']
        r = data['Velocity (m/s)']
        A = data['App_gen (V)']
        ax.plot(f,r,'.-',ms=5,c = cmap(A/A_max))
    polish(fig,1,name=f'{resonator}',extension='.png')
        
        



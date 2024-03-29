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
folder = r'D:\OneDrive_copy\OneDrive - Univerzita Karlova\DATA\2023-03-Helmholtz_resonators\Helmholtz_res_drive_dep'

temp = 'T1350'
resonators = ['500A','500B','500C','500D']

cmap = plt.get_cmap('viridis')

use_latex()
# plt.rcdefaults()


for resonator in resonators:
    filenames = glob(os.path.join(folder,f'{temp}','fsweeps',f'resonator_{resonator}','*.npy'))
    fig,(axx,axy) = plt.subplots(2,1,sharex=True)
    A_max = np.max([np.load(file_temp,allow_pickle=True).item()['App (V)'] for file_temp in filenames])
    for file in filenames[:10]:
        data = np.load(file,allow_pickle=True).item()
        f = data['freq (Hz)']
        x = data['x (A)']
        y = data['y (A)']
        A = data['App (V)']
        axx.plot(f,x,'.',c = cmap(A/A_max))
        axy.plot(f,y,'.',c = cmap(A/A_max))
        axx.set_xlabel(r"sometext with some equation $\sum_{i=0}^{\inf} \Omega (i)$")
    polish(fig,1,extension='.png')



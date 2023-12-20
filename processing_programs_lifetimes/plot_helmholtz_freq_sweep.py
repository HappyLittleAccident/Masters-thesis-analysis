# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 10:28:31 2023

@author: sflab
"""
import numpy as np
import matplotlib.pyplot as plt
import os.path as path
from glob import glob
from matplotlib import cm

folder = r'T1450_2/fsweeps'
basedir = fr'D:\OneDrive\OneDrive - Univerzita Karlova\DATA\2023-03-Helmholtz_resonators\Helmholtz_res_drive_dep\{folder}'

date = '*'
letter = 'D'
distance = 500

resdir = path.join(basedir, 'resonator_{}{}'.format(distance, letter))
files = glob(path.join(resdir, fr'*{date}*.npy'))
# files.sort(key=lambda f: np.load(f, allow_pickle=True).item()['App (V)'])

plt.close('all')
# fig, ax = plt.subplots(3, 1,sharex = True)
fig, ax = plt.subplots()

cmap = plt.get_cmap('copper')

As = []
Rs = []
i=0

for file in files[:]:
    i+=1
    print(file)
    d = np.load(file, allow_pickle=True).item()
    if not isinstance(d['x (A)'][0],list): 
        f = np.array(d['freq (Hz)'])
        x = np.array(d['x (A)'])
        y = np.array(d['y (A)'])
        R = np.abs(x+1j*y)/3e-8
        A = d['App (V)']
        As.append(A)
        Rs.append(R.max() - R[-1:].mean())
        #print(d['bias (V)'])
       
        if A > 0 and i % 3==1:
            # ax[0].plot(f, x, '-o', ms=2, color=cmap(A/3.5))
            # ax[1].plot(f, y, '-o', ms=2, color=cmap(A/3.5))
            # ax[2].plot(f, R, '-o', ms=2, color=cmap(A/3.5))
            im = ax.plot(f,R,'-', ms=2, color=cmap(A/3.5))    
        
    print(A)


fig.colorbar(cm.ScalarMappable(norm=None, cmap=cmap))
fig.suptitle(f'{letter}{distance} nm')
ax.set_xlabel('Frekvence / Hz')
ax.set_ylabel('Rychlost / a.u.')
# f = fig.ginput(1)
#print(f)
fig1,ax1 = plt.subplots()
ax1.plot(As, Rs, 'o')

# -*- coding: utf-8 -*-
"""
Created on Fri Apr 21 20:36:56 2023

@author: Marek
"""

import numpy as np
import matplotlib.pyplot as plt
import os.path as path
from glob import glob
from scipy import interpolate
import os
from scipy.ndimage import gaussian_filter
from format_figure import use_latex,polish


plt.close('all')
basefolder = r'D:\OneDrive_copy\OneDrive - Univerzita Karlova\DATA\2023_4_Helmholtz_resonators_recal_2'
use_latex()

switching_sweeps_names = {'1.35':'T1350','1.45':'T1450','1.65':'T1650','1.85':'T1850'}
temps = ['1.35','1.45','1.65','1.85']
letters = ['A','B','C','D']
fig,ax = plt.subplots(3,2,sharey=True,dpi=250,height_ratios = [0.05,1,1])
gs = ax[0,0].get_gridspec()
for a in ax[0,0:]:
    a.remove()
axcbar = fig.add_subplot(gs[0, 0:])    
# dummy = fig.add_subplot(2,1,2,frameon=True)
# dummy.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
# dummy.grid(visible=None)
fig.subplots_adjust(hspace=0.2, wspace=0.1)
for letter in letters:
    print(f'Resonator: {letter}')
    if letter in ['B']:
        for i,temp in enumerate(temps[:]):

            print(f'Temperature: {temp}')
            date = '*'
            distance = 500
               
            folder = switching_sweeps_names[temp]
            basedir = fr'{basefolder}\Helmholtz_buffer_measurements\{folder}\buffer'
        
            resdir = path.join(basedir,r'resonator_{}{}'.format(distance,letter))
            files = glob(path.join(resdir, r'*.npy'))            
            files.sort(key=lambda f: np.load(f, allow_pickle=True).item()['App_gen (V)'])
            
            amps = []

            amp_last = 0
            vs=[]
            ps=[]
            
            for f in files:
                data = np.load(f, allow_pickle=True).item()
                
                r = np.array(data['Velocity (m/s)'])*100
                P = np.array(data['Pressure (Pa/m)'])/1000
                amp = np.array(data['App_gen (V)'])            
                
                p = np.full_like(r,P)
                amps.append(amp)     
            
                ps.append(p)
                vs.append(r)
                amp_last = amp
                                    
            ps = np.concatenate(ps)
            vs = np.concatenate(vs)
            cmap = plt.get_cmap('magma_r')
                    
            vs_hists = []
            max_vs = np.max(vs)
            min_vs = np.min(vs)
            
            min_ps = np.min(ps)
            max_ps = np.max(ps)
            # min_ps = 2
            # max_ps = 17
            
            max_vs = 30.9 #m/s
            min_vs = 13.1 #m/s
            

            hist,xedges,yedges = np.histogram2d(ps,vs,bins=[(-(i<2)*50 +100),100],range = [[min_ps,max_ps],[min_vs,max_vs]])
            vs_values=0.5*(yedges[1:]+yedges[:-1])
            vs_edges=yedges
            ps_values=0.5*(xedges[1:]+xedges[:-1])
            
            Ps,Vs = np.meshgrid(ps_values,vs_values)
            hist = gaussian_filter(hist,0.1)
            hist = hist.T
            # ax.pcolormesh(Ps,Vs,hist.T,cmap = cmap,shading = 'nearest',norm='linear')
            first = i % 2
            second =1+ i // 2
            

            h = ax[second][first].imshow(hist, extent=[xedges.min(),xedges.max(), yedges.min(), yedges.max()], origin='lower',
                aspect='auto', norm='linear', cmap='gist_heat_r', vmin=0, vmax=3e4)
fig.colorbar(h, ax=None,cax = axcbar, label='Counts',orientation='horizontal',location='top')
for a in ax[:,:].flatten():
    a.tick_params(top=False, labeltop=False, bottom=True, labelbottom=True)
for a in ax[:,1]:
    a.tick_params(left = False, right = False,labelleft=False, labelright=False)

# fig.tight_layout()    
# axcbar.set_aspect(1)
fig.supxlabel('Pressure gradient (kPa/mm)')
fig.supylabel('Superfluid velocity (cm/s)')
# ax.set_xlabel('Pressure (Pa/mm)')
# ax.set_ylabel('Velocity (cm/s)')
# ax.set_title(f'Temperature {temp}')

polish(fig, 1, name='images//switching', extension='.png', grid=True,width_to_height = 0.8,tight_layout=False)        


#%% plot single measurements
plt.close('all')
basefolder = r'D:\OneDrive_copy\OneDrive - Univerzita Karlova\DATA\2023_4_Helmholtz_resonators_recal_2'
switching_sweeps_names = {'1.35':'T1350','1.45':'T1450','1.65':'T1650','1.85':'T1850'}
temps = ['1.35','1.45','1.65','1.85']

plt.rcdefaults()
indices = [45,95,100,200,220,300,200,250]

for i,temp in enumerate(temps):
    folder = basefolder+f'\\Helmholtz_buffer_measurements\\{switching_sweeps_names[temp]}\\buffer\\resonator_500B\\*npy'
    files = glob(folder)
    A_last = 0
    fig,ax= plt.subplots(1,2,sharey= True)
    len_files = len(files)
    print('len: ',len_files)
    files.sort(key = lambda x: np.load(x,allow_pickle=True).item()['App_gen (V)'])


    d = np.load(files[indices[2*i]],allow_pickle=True).item()
    # print(d.keys())
    t = d['time (s)']
    vs = d['Velocity (m/s)']
    P = d['Pressure (Pa/m)']

    print('P1',P)
    ax[0].plot(t,vs,'ko',ms=0.5)
    
    d = np.load(files[indices[2*i+1]],allow_pickle=True).item()
    t = d['time (s)']
    vs = d['Velocity (m/s)']    
    P = d['Pressure (Pa/m)']

    print('P2',P)
    ax[1].plot(t,vs,'ko-',ms=0.5,lw=0.1)
    fig.supxlabel('Time (s)')
    fig.supylabel('Superfluid velocity')
    fig.suptitle(f'Temperature = {temp}K')





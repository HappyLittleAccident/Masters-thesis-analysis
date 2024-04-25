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
folder = r'Helmholtz_res_drive_dep'
temps = ['T1325', 'T1350', 'T1375', 'T1400', 'T1425', 'T1450', 'T1475',
         'T1500', 'T1550', 'T1650', 'T1700', 'T1750', 'T1850', 'T1950']
resonators = ['500A', '500B', '500C']
names = ['C', 'J', 'R']

# use_latex()
plt.rcdefaults()

cmap = plt.get_cmap('winter')
viridis = plt.get_cmap('viridis')



for i,resonator in enumerate(resonators[0:]):
    fig,ax = plt.subplots(1,1,sharex=True,dpi=200)
    for temp in temps[0::]:
        if temp in ['T1325','T1650','T1850']:
            scale = (float(temp[1:])*1e-3 - 1.325)/(1.95 - 1.325)
    
            filenames = glob(os.path.join(folder,f'{temp}','ampsweeps_fast',f'resonator_{resonator}_updowndown','*.npy'))
            # A_max = np.max([np.load(file_temp,allow_pickle=True).item()['App (V)'] for file_temp in filenames])
            for file in filenames[5:6]:
                data = np.load(file,allow_pickle=True).item()
                   
                upP = data['upPressure (Pa/m)']
                downP = data['downPressure (Pa/m)']
                jumpP = data['jumpPressure (Pa/m)']
                up = data['upVelocity (m/s)']*100
                down =  data['downVelocity (m/s)']*100
                jump = data['jumpVelocity (m/s)']*100
                if resonator == '500C':
                    upP/=25000
                    downP/=25000
                    up/=11
                    down/=11
                
                ax.plot(upP,up,'o-',ms=0.5,lw=0.5,c='tab:blue',alpha=1)
                ax.plot(downP,down,'o-',ms=0.5,lw=0.2,c='tab:orange',alpha=0.3)
                
                maxx = ax.get_xlim()[1]-ax.get_xlim()[0]
                maxy = ax.get_ylim()[1]-ax.get_ylim()[0]
                z = (upP[-1] - upP[upP > upP[-1]*0.8][0])/maxx + 1j*(up[-1] - up[upP > upP[-1]*0.8][0])/maxy               
                ax.text(upP[-1]*0.8,up[upP > upP[-1]*0.8][0]+maxy*0.025,temp[1:2]+'.'+temp[2:]+'K',rotation=np.angle(z,deg=True))
                
                
                
    ax.set_xlabel('Pressure gradient (Pa/m)')
    ax.set_ylabel('Superfluid velocity (cm/s)')
    fig.tight_layout()
    # fig.colorbar(plt.cm.ScalarMappable(
    # norm=mpl.colors.Normalize(vmin=1.325, vmax=1.85, clip=False), cmap=cmap),
    # ax=ax, label='Temperature (K)')


    polish(fig, 1, name=f'images//ampsweeps_single{names[i]}', extension='.png', grid=True,width_to_height = 0.8)        
#%%
for i,resonator in enumerate(resonators[0:]):
    fig,ax = plt.subplots(1,1,sharex=True,dpi=200)
    for temp in temps[0::]:
        scale = (float(temp[1:])*1e-3 - 1.325)/(1.95 - 1.325)

        filenames = glob(os.path.join(folder,f'{temp}','ampsweeps_fast',f'resonator_{resonator}_updowndown','*.npy'))
        # A_max = np.max([np.load(file_temp,allow_pickle=True).item()['App (V)'] for file_temp in filenames])
        for file in filenames[:]:
            data = np.load(file,allow_pickle=True).item()
            upP = data['upPressure (Pa/m)']
            downP = data['downPressure (Pa/m)']
            jumpP = data['jumpPressure (Pa/m)']
            up = data['upVelocity (m/s)']*100
            down =  data['downVelocity (m/s)']*100
            jumpt = data['jumpVelocity (m/s)']*100
            ax.plot(upP,up,'o-',ms=0.5,lw=0.5,c=cmap(scale),alpha=1)
            # ax.plot(downP,down,'o-',ms=0.5,lw=0.2,c=,alpha=0.3)
            
            
    ax.set_xlabel('Pressure gradient (Pa/m)')
    ax.set_ylabel('Superfluid velocity (cm/s)')
    fig.tight_layout()
    fig.colorbar(plt.cm.ScalarMappable(norm=mpl.colors.Normalize(
    vmin=1.325, vmax=1.85, clip=False), cmap=cmap), ax=ax, label='Temperature (K)')


    polish(fig, 1, name=f'images//ampsweeps{names[i]}', extension='.png', grid=True,width_to_height = 1.3)        




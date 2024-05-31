# -*- coding: utf-8 -*-
"""
Created on Wed May 22 16:27:30 2024

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

# plt.ioff()
plt.ion()

for i,resonator in enumerate(resonators[1:2]):
    fig,ax = plt.subplots(3,1,sharex=False,dpi=250,gridspec_kw={'height_ratios': [1,1, 0.05]})


    for j,temp in enumerate(temps[0::]):
        Vs = []
        scale = (float(temp[1:])*1e-3 - 1.325)/(1.95 - 1.325)

        filenames = glob(os.path.join(folder,f'{temp}','ampsweeps_fast',f'resonator_{resonator}_updowndown','*.npy'))
        # A_max = np.max([np.load(file_temp,allow_pickle=True).item()['App (V)'] for file_temp in filenames])
        vs = []
        ps = []    
        for file in filenames[:]:
            data = np.load(file,allow_pickle=True).item()
            if len(data['upPressure (Pa/m)']) != 2000:
                continue
            
            indent= j*(700+ (i==1)*5000 + (i==2)*1500)*0
            upP = data['upPressure (Pa/m)'] + indent
            downP = data['downPressure (Pa/m)'] + indent
            jumpP = data['jumpPressure (Pa/m)'] + indent
            upV = data['upVelocity (m/s)']*100 
            downV =  data['downVelocity (m/s)']*100
            jumpV = data['jumpVelocity (m/s)']*100  
            
            if len(upV)==2000:
                Vs.append(upV)
                # Vs.append(downV[::-1])
            else:
                continue
            
            if resonator == '500C':
                upP/=25000
                downP/=25000
                jumpP/=25000
                upV/=11
                downV/=11
                jumpV/=11
                if j>len(temps)-3:
                    continue
            
            

            vs.append(upV)
            vs.append(downV)
            # print(temp)

            ps.append(upP)
            ps.append(downP)
            ax[1].plot(upP,upV,'.',markeredgewidth=0.0,ms=0.5,c=cmap(scale),alpha=0.2)
            ax[1].plot(downP,downV,'.',markeredgewidth=0.0,ms=0.5,c=cmap(scale),alpha=0.2)
            
        Vs = np.column_stack(Vs)
        
        std_V = np.std(Vs,axis = 1)
        ax[0].plot(upP,std_V/np.mean(Vs,axis=1),c=cmap(scale))
        ax[0].set_xlabel('Pressure grad (Pa/m)')
        ax[0].set_ylabel(r'Standard deviation of vs ($\mathrm{cm s^{-1}}$)')
        fig.canvas.draw()
        fig.tight_layout()
        while not plt.waitforbuttonpress():
            pass
        ax[0].clear()
        ax[1].clear()
    if i==2:
        ax[1].set_xlabel('Normalised pressure gradient (Arb.)')
        ax[1].set_ylabel('Normalised superfluid velocity (Arb.)')
    else:
        ax[1].set_xlabel('Pressure gradient (Pa/m)')
        ax[1].set_ylabel('Superfluid velocity (cm/s)')
    
    
    fig.colorbar(plt.cm.ScalarMappable(norm=mpl.colors.Normalize(
    vmin=1.325, vmax=1.95, clip=False), cmap=cmap),orientation='horizontal',location="bottom", ax=None,cax=ax[2], label='Temperature (K)')
    # 

    # polish(fig, 1, name=f'images//ampsweeps{names[i]}', extension='.png', grid=True,width_to_height = 1.3)        
    


    polish(fig, 1, name=f'images//ampsweeps{names[i]}', extension='.png', grid=True,width_to_height = 0.5,tight_layout = True)        



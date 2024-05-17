# -*- coding: utf-8 -*-
"""
Created on Thu Jul 27 11:22:38 2023

@author: Marek
"""
import numpy as np
import matplotlib.pyplot as plt
import os
from glob import glob
import matplotlib as mpl
from format_figure import use_latex,polish

plt.close('all')

folder = r'D:\OneDrive_copy\OneDrive - Univerzita Karlova\DATA\2023_4_Helmholtz_resonators_recal_2\Helmholtz_res_VLD'


temperature_folders = glob(folder+'\\*')
temps_strings = [name[-4:] for name in temperature_folders]

resonators = ['500A', '500B', '500C']
names = ['C', 'J', 'R']

use_latex()
# plt.rcdefaults()

cmap = plt.get_cmap('winter')
viridis = plt.get_cmap('viridis')



for i,resonator in enumerate(resonators[0:]):
    fig_hist,ax_hist = plt.subplots(dpi=200)
    for j,temp in enumerate(temps_strings[:]):
        filenames = glob(folder + f'\\T{temp}\\{resonator}\\*.npy')
    
        scale = (float(temp[:])*1e-3 - 1.325)/(1.95 - 1.325)
        for filename in filenames[:]:
            
            data = np.load(filename,allow_pickle= True).item()
            
            indent= j*(10 - ((i==1)*0.06+0.32)*j + (i==1)*20 + (i==2)*1)*5e8
            
            
            upV = data['upV (m/s)']*100 
            downV = data['downV (m/s)']*100
            jumpV = data['jumpV (m/s)']*100
            
            
            upL = data['upL (m^-2)']     + indent
            downL = data['downL (m^-2)'] + indent
            jumpL = data['jumpL (m^-2)'] + indent

            
            if resonator == '500C':
                upV/=11
                downV/=11
                jumpV/=11
                
                upL/=5e10
                downL/=5e10
                if j>len(temps_strings)-4:
                    continue
                
            ax_hist.plot(upV,upL,'.',markeredgewidth=0.0,ms=0.5,c=cmap(scale),alpha=0.2)
            ax_hist.plot(downV,downL,'.',markeredgewidth=0.0,ms=0.5,c=cmap(scale),alpha=0.2)
            # ax_hist.semilogy(upV,upL,'.',markeredgewidth=0.0,ms=0.5,c=cmap(scale),alpha=0.2)
            
    if i==2:
        ax_hist.set_ylim(-0.1,1)
        ax_hist.set_xlabel('Normalised superfluid velocity (Arb.)')
        ax_hist.set_ylabel('Normalised VLD (Arb.)')
    else:
        ax_hist.set_xlabel('Superfluid velocity (cm/s)')
        ax_hist.set_ylabel('VLD ($m^{-2}$)')
        
    fig_hist.colorbar(plt.cm.ScalarMappable(norm=mpl.colors.Normalize(
    vmin=1.325, vmax=1.95, clip=False), cmap=cmap), ax=ax_hist, label='Temperature (K)')
    # fig_hist.tight_layout()
    

    polish(fig_hist, 1, name=f'images//VLD{names[i]}', extension='.png', grid=True,width_to_height = 2.5,tight_layout=True)        


    
            
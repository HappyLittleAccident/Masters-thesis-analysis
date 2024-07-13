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
# 
use_latex()
# plt.rcdefaults()

cmap = plt.get_cmap('winter')
viridis = plt.get_cmap('viridis')


fig_hist,ax_hist = plt.subplots(1,3,dpi=200)
fig_comp,ax_comp = plt.subplots(1,2,dpi=300)

for i,resonator in enumerate(resonators[0:3]):
    
    
    for j,temp in enumerate(temps_strings[:]):

        filenames = glob(folder + f'\\T{temp}\\{resonator}\\*.npy')
    
        scale = (float(temp[:])*1e-3 - 1.325)/(1.95 - 1.325)
        for filename in filenames[9:10]:
            
            
            data = np.load(filename,allow_pickle= True).item()
            
            indent= j*(10 - ((i==1)*0.06+0.32)*j + (i==1)*20 + (i==2)*1)*5e8
            
            
            upV = data['upV (m/s)']*100 
            downV = data['downV (m/s)']*100
            jumpV = data['jumpV (m/s)']*100
            
            
            upL = data['upL (m^-2)']     
            downL = data['downL (m^-2)'] 
            jumpL = data['jumpL (m^-2)']
            
            
            
            if i<2 and temp == '1500':
                ax_comp[1].plot(upV+7*i,upL,label='Resonator ' + ('C' if i==0 else 'J'))


            upL+=indent
            downL+=indent
            
            if resonator == '500C':
                upV/=11
                downV/=11
                jumpV/=11
                
                upL/=5e10
                downL/=5e10
                if j>len(temps_strings)-4:
                    continue
                
            ax_hist[i].plot(upV,upL,'.',markeredgewidth=0.0,ms=0.5,c=cmap(scale),alpha=0.2)
            ax_hist[i].plot(downV,downL,'.',markeredgewidth=0.0,ms=0.5,c=cmap(scale),alpha=0.2)
            # ax_hist.semilogy(upV,upL,'.',markeredgewidth=0.0,ms=0.5,c=cmap(scale),alpha=0.2)

        
    if i==2:
        ax_hist[i].set_ylim(-0.1,1)
        # ax_hist[i].set_xlabel('Normalised superfluid velocity (Arb.)')
        # ax_hist[i].set_ylabel('Normalised VLD (Arb.)')
    else:

        
        # ax_hist[i].set_xlabel('Superfluid velocity (cm/s)')
        ax_hist[0].set_ylabel('VLD ($\mathrm{m^{-2}}$)')
        

fig_hist.supxlabel('Superfluid velocity (cm/s)')
fig_hist.colorbar(plt.cm.ScalarMappable(norm=mpl.colors.Normalize(
vmin=1.325, vmax=1.95, clip=False), cmap=cmap), ax=ax_hist, label='Temperature (K)',location = 'bottom')
    # fig_hist.tight_layout()
#%%

ax_comp[1].axvline(31.8,c='k',lw=0.5)
ax_comp[1].axvline(34.8,c='k',lw=0.5)
gs = np.linspace(0,6,1000)
Ls = np.linspace(0,6,1000)

GS,LS = np.meshgrid(gs,Ls)

def Ldot(L,d,g):
    return -L -d*(L/(L**2 + 1))**2 + g

colors = ['tab:blue','tab:orange']

ds = [6,15]
for i,d in enumerate(ds):
    condition = np.logical_or(Ldot(LS-0.1,d,GS) < 0,Ldot(LS+0.1,d,GS) > 0)
    
    result = np.ma.MaskedArray(Ldot(LS,d,GS),condition)
    result_inverse=  np.ma.MaskedArray(Ldot(LS,d,GS),np.logical_not(condition))
    
    ax_comp[0].contour(GS,LS,result,levels=[0],colors=[colors[i]],linewidths=[1])
    ax_comp[0].contour(GS,LS,result_inverse,levels=[0],colors=[colors[i]],linestyles=['--'],linewidths=[1])
    ax_comp[0].plot([],[],'-',c=colors[i],label=f'$d = {d}$')

ax_comp[0].axvline(4.27,lw=0.5,c='k')
ax_comp[0].axvline(4.85,lw=0.5,c='k')



ax_comp[0].set_xlabel('Vortex generation rate $g$')
ax_comp[0].set_ylabel('Dimensionless vortex line density $L$')
ax_comp[0].legend()

ax_comp[0].annotate('bistable regime',[4.32,4.5],[0.1,4.25],xycoords='data',arrowprops=dict(arrowstyle="simple",color='k'))

ax_comp[1].set_xlabel('Superfluid velocity (cm/s)')
ax_comp[1].set_ylabel('Vortex line density ($\mathrm{m^{-2}}$)')
ax_comp[1].legend()


polish(fig_hist, 1, name=f'images//VLD{names[i]}', extension='.png', grid=True,width_to_height = 6,tight_layout=False)        
polish(fig_comp, 1, name='images//compare_VLD', extension='.png', grid=False,width_to_height = 1.4,tight_layout=True)

    
            
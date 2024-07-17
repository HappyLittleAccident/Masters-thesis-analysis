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
from scipy.optimize import curve_fit


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
            ax[second][first].annotate(f'T = {temp}K',(0.1,0.9),xycoords='axes fraction')
            
fig.colorbar(h, ax=None,cax = axcbar, label='Counts',orientation='horizontal',location='top')
for a in ax[:,:].flatten():
    a.tick_params(top=False, labeltop=False, bottom=True, labelbottom=True)
for a in ax[:,1]:
    a.tick_params(left = False, right = False,labelleft=False, labelright=False)

# fig.tight_layout()    
# axcbar.set_aspect(1)
fig.supxlabel('Pressure gradient (Pa/mm)')
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

use_latex()
# plt.rcdefaults()
fig_dist,ax_dist = plt.subplots(dpi=300)
for i,temp in enumerate(temps[0:3]):

    folder =  basefolder+f'\\Helmholtz_buffer_measurements\\{switching_sweeps_names[temp]}\\buffer\\resonator_500B\\*npy'
    folder_populations = rf'D:\OneDrive_copy\OneDrive - Univerzita Karlova\DATA\switching_analysis_code\files_populations\{switching_sweeps_names[temp]}\buffer\resonator_500B'
    files = glob(folder)
    A_last = 0

    len_files = len(files)
    print('len: ',len_files)
    files.sort(key = lambda x: np.load(x,allow_pickle=True).item()['App_gen (V)'])

    
    
    if i == 2:
        index = 250
        fig,ax= plt.subplots(1,1,dpi=300)
        name = files[index].split('\\')[-1]
        d = np.load(files[index],allow_pickle=True).item()
        d_pops = np.load(folder_populations+'\\'+name,allow_pickle=True).item()
    
        
    
        # print(d.keys())
        t = d['time (s)']
        vs = d['Velocity (m/s)']*100
        P = d['Pressure (Pa/m)']
        populations = d_pops['populations']*(np.max(vs)-np.min(vs))+np.min(vs)
        
        ax.axhline(23.2,c='k',ls='--')
        ax.plot(t,vs,'ko',ms=0.5)
        ax.plot(t,populations,c='tab:blue',lw=0.8)
        ax.set_xlim(14,19)
        ax.set_ylim(20.5,26.0)
        ax.set_xlabel('Measured time (s)')
        ax.set_ylabel('Superflid velocity (cm/s)')
        ax.plot([15.19,15.36],[21.4,21.4],'|-r',mew=2)
        # plt.annotate('', xy=(14.5,21.4), xytext=(15,21.4), arrowprops=dict(arrowstyle='<->'))
        ax.annotate('Threshold',(15.45,22.7))
        ax.annotate('$t_d$',(15.2,21),color='r')
        
        polish(fig, 1, name='images//datapoints', extension='.png', grid=True,width_to_height = 0.5,tight_layout=True)        

        
    name_lifetimes = fr'D:\OneDrive_copy\OneDrive - Univerzita Karlova\DATA\2023_4_Helmholtz_resonators_recal_2\Helmholtz_buffer_lifetimes\resonator_500B_{temp}.npy'

     
    d_lifetimes = np.load(name_lifetimes,allow_pickle=True).item()
    # lifetimes_up = d_lifetimes['series lifetimes up (s)']
    lifetimes_down = d_lifetimes['series lifetimes down (s)']
    mean_lifetimes = d_lifetimes['mean lifetimes down (s)']

    first =True
    for j,series in enumerate(lifetimes_down):

        if mean_lifetimes[j] > 0.3 and first:
            ax_dist.hist(series,bins=30,histtype='step',label=f'T = {temp}K')
            break
    
    
    # d = np.load(files[indices[2*i+1]],allow_pickle=True).item()
    # t = d['time (s)']
    # vs = d['Velocity (m/s)']    
    # P = d['Pressure (Pa/m)']

    # print('P2',P)
    # ax[1].plot(t,vs,'ko-',ms=0.5,lw=0.1)
    # fig.supxlabel('Time (s)')
    # fig.supylabel('Superfluid velocity')
    # fig.suptitle(f'Temperature = {temp}K')
ax_dist.semilogy()
ax_dist.legend()
ax_dist.set_xlim(0,2)
ax_dist.set_xlabel('Lifetime of dissipative state $t_{d}$ (s)')
ax_dist.set_ylabel('Number of occurances')
polish(fig_dist, 1, name='images//histogram', extension='.png', grid=True,width_to_height = 0.8,tight_layout=False)        


#%%
plt.close('all')
fig,ax = plt.subplots(1,2,dpi = 400)
vs0s = [27.476,25.745,20.625]



# Make a horizontal slider to control the frequency.

use_latex()
# plt.rcdefaults()
nu = -1.33


fig.tight_layout()
# 
colors = ['blue','orange','green']
styles = ['.','x','+']
for i,temp in enumerate(temps[1:]):
    ax[0].axvline(vs0s[i],c=f'tab:{colors[i]}',ls='--')
for i,temp in enumerate(temps[1:]):
    def tau(vs):
        return np.abs(vs)**nu
    filename = fr'D:\OneDrive_copy\OneDrive - Univerzita Karlova\DATA\2023_4_Helmholtz_resonators_recal_2\Helmholtz_buffer_lifetimes\resonator_500B_{temp}.npy'
    d = np.load(filename,allow_pickle = True).item()
    vs_down = d['V down (m/s)']*100 #cm/s
    tau_d   = d['mean lifetimes down (s)']
    
    indices = vs_down > 0
    vs_down = vs_down[indices]
    tau_d   = tau_d[indices]
    
    # par,sig = curve_fit(tau,vs_down,np.log(tau_d),p0 = [0])
    t = np.linspace(-30,0,500)
    # print(par)plt.
    
    # ax[0].semilogy(t,np.exp(tau(t,*par)),'k--')    
    ax[0].semilogy(vs_down,tau_d,      styles[i],ms=3,mew=0.5,label=f'{temp}K')
    ax[1].loglog(vs0s[i]-vs_down,tau_d,styles[i],ms=3,mew=0.5,label=f'{temp}K')

    

    
ax[1].loglog(-t,tau(t)*0.2,'k--')#,label = '$(v_0-v_s)^{-\\eta}$')
ax[1].set_xlim(1e-3,45)
ax[1].annotate('$(v_c-v_s)^{-\\eta}$',xy=(0.05,7),xytext=(0.0012,0.12),arrowprops=dict(arrowstyle='->'))
ax[1].legend(loc='lower left',frameon=True)
ax[1].set_xlabel('$v_c - v_s$ (cm/s)')

ax[0].set_ylabel('Average state lifetime $\\tau_d$ (s)')
ax[0].set_xlabel('Superfluid velocity $v_s$ (cm/s)')
ax[0].set_ylim(1e-4,50)
ax[0].legend(loc='lower right')

polish(fig, 1, name='images//mean_lifetimes', extension='.png', grid=True,width_to_height = 1.2,tight_layout=True)        

#%%
plt.rcdefaults()
plt.close('all')
from matplotlib.widgets import Slider
for i,temp in enumerate(temps[1:]):
    def tau(vs):
        return np.abs(vs)**nu
    filename = fr'D:\OneDrive_copy\OneDrive - Univerzita Karlova\DATA\2023_4_Helmholtz_resonators_recal_2\Helmholtz_buffer_lifetimes\resonator_500B_{temp}.npy'
    d = np.load(filename,allow_pickle = True).item()
    vs_down = d['V down (m/s)']*100 #cm/s
    tau_d   = d['mean lifetimes down (s)']
    
    fig_slider,ax_slider = plt.subplots()
    fig_slider.subplots_adjust(bottom=0.25)
    axvs0 = fig_slider.add_axes([0.25, 0.1, 0.65, 0.03])
    vs_slider = Slider(
        ax=axvs0,
        label='Velocity',
        valmin=20,
        valmax=30,
        valinit=20,
    )

    def update(val):
        ax_slider.clear()
        ax_slider.loglog(val-vs_down,tau_d,'.')
        ax_slider.loglog(-t,tau(t),'k--')
        ax_slider.set_title(val)
        fig_slider.canvas.draw_idle()
        
    vs_slider.on_changed(update)

    while not plt.waitforbuttonpress():
        pass
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 21 20:36:56 2023

@author: Marek
"""

import numpy as np
import matplotlib.pyplot as plt
import os.path as path
from glob import glob
import scipy.signal as sig
from scipy.optimize import curve_fit

plt.close('all')


def function(t,a,b,c,d):
    return d*t**3 + b*t**2 + c*t + d

switching_sweeps_names = {'1.35':'T1350_3','1.45':'T1450_3','1.65':'T1650_3','1.85':'T1850_2'}
switching_ampsweeps_names = {'1.35':'T1350_2','1.45':'T1450_2','1.65':'T1650_2','1.85':'T1850_1'}

date = '*'
letter = 'B'
distance = 500
temp = '1.65'

plot_hist = False

folder = switching_sweeps_names[temp]
folder_ampsweeps = switching_ampsweeps_names[temp]



basedir = r'D:\Github\Masters-thesis-analysis\filip_calibration\Helmholtz_res_drive_dep\{}\buffer'.format(folder)

ampsweeps_base = r'D:\OneDrive\OneDrive - Univerzita Karlova\DATA\2023-03-Helmholtz_resonators\Helmholtz_res_drive_dep\{}\ampsweeps_fast'.format(folder_ampsweeps)

resdir = path.join(basedir,r'resonator_{}{}'.format(distance,letter))
files = glob(path.join(resdir, fr'*{date}*.npy'))

ampsweeps = path.join(ampsweeps_base,r'resonator_{}{}_updowndown'.format(distance,letter))
files_ampsweeps = glob(path.join(ampsweeps, fr'*{date}*.npy'))

# files = files[:-53]
files.sort(key=lambda f: np.load(f, allow_pickle=True).item()['App_gen (V)'])



fig_ampsweeps, ax_ampsweeps = plt.subplots()


cmap = plt.get_cmap('viridis')
cmap_hist = plt.get_cmap('Reds')

#%% plot ampsweeps
for file in files_ampsweeps[:]:#files[:20] and files[39:]:
    # print(file)
    d = np.load(file, allow_pickle=True).item()
    try:
        if not isinstance(d['x (A)'][0],list):    
            f = np.array(d['freq (Hz)'])
            x = np.array(d['x (A)'])
            y = np.array(d['y (A)'])
            A = np.array(d['App (V)'])
            A+=5.0
            A /= 10
            base = d['App_base (V)']
            N = len(x)
            leg1 = np.arange(0, int(N/3))
            leg2 = np.arange(int(N/3), int(2*N/3))
            leg3 = np.arange(int(2*N/3), N)
            ax_ampsweeps.plot(A[leg1], np.sqrt(x**2 + y**2)[leg1], '-o',ms=0.5, lw=0.5, color='tab:blue',alpha=1)
            ax_ampsweeps.plot(A[leg2], np.sqrt(x**2 + y**2)[leg2], '-o',ms=0.5, lw=0.5, color='tab:orange',alpha=1)
            ax_ampsweeps.plot(A[leg3], np.sqrt(x**2 + y**2)[leg3], '-o',ms=0.5, lw=0.5, color='tab:green',alpha=1)
    except:
        pass


#%% analyse switching
amplitudes = []
lifetimes_up=[]
lifetimes_down=[]

higher_1 = []
lower_1 = []

amp_last = np.load(files[0], allow_pickle=True).item()['App_gen (V)']
for file in files: #and files[39:]:
    # print(file)
    d = np.load(file, allow_pickle=True).item()
    if not isinstance(d['x (A)'][0],list):
        
        f = np.array(d['res freq (Hz)'])
        x = d['x (A)']
        y = d['y (A)']
        t = d['time']
        A = d['App_gen (V)']
        tc = d['time constant (s)']
        sample_rate = d['sample rate (Hz)']
        print(f'Date: {file[-27:-8]} \t Drive amplitude: {A:.4f}Vrms')
    
        r = np.abs(x+1j*y)
        
        qty = x
        N = len(qty)
        qty0 = 0.5*(np.max(qty)+np.min(qty))
        population = np.where(qty>qty0,np.ones_like(qty),np.zeros_like(qty))
        
        higher = []
        lower = []
        last = -1
        for pop in population[:]:
            if pop==1:
                if last == 1:
                    higher[-1]+=1
                else:
                    higher.append(1)
            else:
                if last == 0:
                    lower[-1]+=1
                else:
                    lower.append(1)
            last = pop
            

                

        
        if A==amp_last:
            higher_1 = higher
            lower_1 = lower

        else:

            amplitudes.append(A)
            higher_2 = np.concatenate([np.array(higher),np.array(higher_1)])
            lower_2 = np.concatenate([np.array(lower),np.array(lower_1)])
            
            lifetime_up = np.mean(higher_2)/sample_rate
            lifetime_down = np.mean(lower_2)/sample_rate
            lifetimes_up.append(lifetime_up)
            lifetimes_down.append(lifetime_down)
            
            
            if plot_hist:
                fig_hist, ax_hist = plt.subplots(2,1)
                
                ax_hist[0].hist(higher_2,bins=60, color='tab:blue')
                ax_hist[1].hist(lower_2,bins=60, color='tab:orange')
    
                ax_hist[0].set_title('Up')
                ax_hist[1].set_title('Down')
                ax_hist[1].set_xlabel('Lifetime (s)')
    
                fig_hist.suptitle(f'HISTOGRAM\n{letter}{distance}|folder: {folder}|Amp: {A:.3f}Vrms|Relative drive: {A/3.5:.3f}')
                fig_hist.tight_layout()
        
        amp_last = A
        
        counts,edges = np.histogram(r,2)
        # amplitudes_plot = np.full_like(r,A)
        # for k,count in enumerate(counts):
            # ax_ampsweeps.plot(amplitudes_plot/base,r, '-o',ms=0.5, color='k')
            # ax_ampsweeps.plot(A/3.5,(edges[k]+edges[k+1])/2, '-o',ms=2, color='k',alpha=count/np.sum(counts))
            # print(f'Len up: {len(higher_2)}')
            # print(f'Len down: {len(lower_2)}')
        ax_ampsweeps.plot(np.full_like(r,A/3.5),r, '-o',ms=1, color='k')
        
        # fig, ax = plt.subplots(3, 1, sharex= True)
        # ax[0].plot(t, x, '-o',ms=0.5, lw=0.5, color='tab:blue')
        # ax[2].plot(t, (np.max(r)-np.min(r))*population + np.min(r), 'o',ms=0.5, lw=0.5, color='tab:orange')
        # ax[2].axhline(qty0, color='tab:red')
        # ax[1].plot(t,  y, '-o',ms=0.5, lw=0.5, color='tab:orange')
        # ax[2].plot(t, r,'-o',ms=0.5, lw=0.5, color='tab:green')
        # fig.suptitle(f'{letter}{distance}|folder: {folder}|type: buffer|SIGNAL|Amp: {A:.3f}Vrms|tc: {tc}s')
        # ax[2].set_xlabel("time (s)")
        # ax[1].set_ylabel("response amplitude (a.u.)")
        # ax[1].set_title('Y component')
        # ax[0].set_title('X component')
        # ax[2].set_title('R component')
        # fig.tight_layout()
        
# ax_ampsweeps.plot(A,np.mean())
    



#%% plot lifetimes 



fig_lifetime, ax_lifetime = plt.subplots()
fig_lifetime_log, ax_lifetime_log = plt.subplots()
crit_drive = 0#amplitudes[np.argmax(lifetimes_down)]/base

ax_lifetime_log.loglog(np.array(amplitudes)/base-crit_drive,lifetimes_up, 'o-', lw=1, color='tab:red',label='up')
ax_lifetime_log.loglog(np.array(amplitudes)/base-crit_drive,lifetimes_down, 'o-', lw=1, color='k',label='down')

ax_lifetime.semilogy(np.array(amplitudes)/base-crit_drive,lifetimes_up, 'o-', lw=1, color='tab:red',label='up')
ax_lifetime.semilogy(np.array(amplitudes)/base-crit_drive,lifetimes_down, 'o-', lw=1, color='k',label='down')

ax_lifetime.set_xlabel("Drive (a.u.)")
ax_lifetime.set_ylabel("Mean state lifetime (s)")
# ax_lifetime.set_title(f'MEAN LIFETIME\n{letter}{distance}|folder: {folder}|type: buffer|tc: {tc}s')
# ax_lifetime.legend()
fig_lifetime.tight_layout()

ax_lifetime_log.set_xlabel("Drive (a.u.)")
ax_lifetime_log.set_ylabel("Mean state lifetime (s)")
# ax_lifetime_log.set_title(f'MEAN LIFETIME\n{letter}{distance}|folder: {folder}|type: buffer|tc: {tc}s')
# ax_lifetime_log.legend()
fig_lifetime_log.tight_layout()












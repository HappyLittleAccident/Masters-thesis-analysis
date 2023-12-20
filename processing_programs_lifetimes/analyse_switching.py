# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 17:03:05 2023

@author: Marek
"""

import numpy as np
import matplotlib.pyplot as plt
import os.path as path
from glob import glob
from scipy import signal
from scipy.optimize import curve_fit
from scipy import fft
from scipy import signal
from scipy.interpolate import Rbf

old_165_B=True

plot = False
plot_dataset = False

def gaussian(x,x0,w,A):
    return A*np.exp(-0.5*((x-x0)/w)**2)/(w*np.sqrt(2*np.pi))

def two_gaussian(x,x0,w,A,x01,w1,A1):
    return gaussian(x,x0,w,A) + gaussian(x,x01,w1,A1)


plt.close('all')

switching_sweeps_names = {'1.35':'T1350_3','1.45':'T1450_3','1.65':'T1650_3','1.85':'T1850_2'}
switching_ampsweeps_names = {'1.35':'T1350_2','1.45':'T1450_2','1.65':'T1650_2','1.85':'T1850'}

temps = ['1.35','1.45','1.65','1.85']
letters = ['A','B','C','D']
if plot:
    fig_hist,ax_hist = plt.subplots()

if plot_dataset:
    fig,ax = plt.subplots()
    fig_corr, ax_corr = plt.subplots()  
cmap = plt.get_cmap('viridis')
cutoffs = [200]
for letter in letters:
    print(f'Resonator: {letter}')
    if letter == 'B':
        for temp in temps[2:3]:
            fig_lifetime, ax_lifetime = plt.subplots()
            # fig_lifetime_log,ax_lifetime_log = plt.subplots()
            for cutoff in cutoffs:

                print(f'Temperature: {temp}')
                date = '*'
                distance = 500
                   
                folder = switching_sweeps_names[temp]
                folder_ampsweeps = switching_ampsweeps_names[temp]
                
                basedir = r'..\2023-03-Helmholtz_resonators\Helmholtz_res_drive_dep\{}\buffer'.format(folder)
                
                ampsweeps_base = r'..\2023-03-Helmholtz_resonators\Helmholtz_res_drive_dep\{}\ampsweeps_fast'.format(folder_ampsweeps)
                
                resdir = path.join(basedir,r'resonator_{}{}'.format(distance,letter))
                files = glob(path.join(resdir, r'*.npy'))
                
                ampsweeps = path.join(ampsweeps_base,r'resonator_{}{}_updowndown'.format(distance,letter))
                files_ampsweeps = glob(path.join(ampsweeps, r'*.npy'))
                
                # files = files[:-53]
                files.sort(key=lambda f: np.load(f, allow_pickle=True).item()['App_gen (V)'])
                
                
                folder = switching_sweeps_names[temp]
                basedir = r'..\2023-03-Helmholtz_resonators\Helmholtz_res_drive_dep\{}\buffer'.format(folder)
                
                params = np.load(f'lifetimes_params\params_{letter}_T{temp}_finished.npy', allow_pickle=True).item()
                y1_list = params['y1']
                y2_list = params['y2']
                dy1_list = params['dy1']
                dy2_list = params['dy2']
                
                
                if temp=='1.65' and letter == 'B' and old_165_B:
                    params = np.load(f'params_{letter}_T{temp}_manual.npy', allow_pickle=True).item()
                    y1_list = params['y1']
                    y2_list = params['y2']
                    dy1_list = params['dy']
        
        
                
                #%% analyse switching
                amplitudes = []
                lifetimes_up=[]
                lifetimes_down=[]
                
                amplitudes_f = []
                lifetimes_up_f=[]
                lifetimes_down_f=[]
                
                higher_stack = []
                lower_stack = []
                
                higher_stack_f = []
                lower_stack_f = []
                mean_rs = []
                y1s = []
                y2s =[]
                
                amp_last = -1
                amp_last_f = -1
                y2_last = 0
                y1_last = 0
                r_last = 0
                for i,file in enumerate(files): #and files[39:]:
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
                        # print(f'Date: {file[-27:-8]} \t Drive amplitude: {A:.4f}Vrms')
                    
                        r = np.abs(x+1j*y)
                        
                        n=50
                        b = [1.0 / n] * n
                        a = 1
                        qty_f = np.copy(r)
                        # rbf = Rbf(t, r, function = 'thin_plate', smooth = 10)
                        # qty_f = rbf(t)
                        qty_f[n:] = signal.lfilter(b,a,r)[n:]
                        w=1
                        # qty_f = np.convolve(r,np.ones(w),mode='same')[w:-w]/w
                        # qty = signal.savgol_filter(r,200,6)
                        sos = signal.butter(8, 2*cutoff/sample_rate, output='sos')
        
                        qty = signal.sosfiltfilt(sos, r)
                        plot_fft = False
        
                        # qty = r
                        if 0:
                            hist,edges = np.histogram(qty,80)
                            
                            qty_values=0.5*(edges[1:]+edges[:-1])
                            if plot:
                                ax_hist.clear()
                                ax_hist.plot(qty_values,hist,'.-k',ms=1)
                            if plot_dataset:
                                ax.clear()
                                ax.plot(t,r,'.-k',ms=2,lw=0.2)
                                ax.plot(t,qty_f,'-r',lw=0.5)
                                ax.axhline(np.mean(r),c='r',ls='-')
                                ax.axhline(np.mean(r)+np.std(r),c='r',ls='--')
                                ax.axhline(np.mean(r)-np.std(r),c='r',ls='--')
                                # try:
                                #     ax.axhline(y1_list[i],c='r',ls='-')
                                #     ax.axhline(y1_list[i]+3*dy1_list[i],c='r',ls='--')
                                #     ax.axhline(y1_list[i]-3*dy1_list[i],c='r',ls='--')
                                #     ax.axhline(y2_list[i]+dy2_list[i],c='b',ls='--')
                                #     ax.axhline(y2_list[i]-dy2_list[i],c='b',ls='--')
                                #     ax.axhline(y2_list[i],c='b',ls='-')
                                # except:
                                    # pass
                                fig.canvas.draw()
                            
                            if plot_fft:
                                ()
                                transformed = fft.fft(x + 1j*y)
                                freq = fft.fftfreq(len(qty),1/sample_rate)
                                
                                
                                for fq1,fq2 in [(25,np.inf)]:
                                    index = np.logical_and(np.abs(freq)>fq1,np.abs(freq)<fq2)
                                    transformed[index]=np.zeros_like(transformed)[index]
                                inverse = np.abs(fft.ifft(transformed))
                                qty=inverse    
                                if plot_dataset:
                                    ax_corr.clear
                                    ax_corr.plot(freq[1:],np.abs(transformed[1:]),'.-k',ms=1)
                                    
        
                                    ax.plot(t,inverse,'.-',c='tab:orange',ms=1)
                                    fig.canvas.draw()
                                    fig_corr.canvas.draw()
                                    
                            """
                            Fit gaussian on data
                            """
                            try:
                                if plot:
                                    p0s = [y1_list[i],(qty_values[0]-qty_values[-1])/3,np.max(hist)]
                                    par,sig = curve_fit(gaussian,qty_values,hist,p0=p0s)
                                
                                    chi = np.sum((gaussian(qty_values,*par)/(hist+1) -1)**2)/len(hist)
                                    skewness= np.mean((qty-np.mean(qty))**3)
                                    if np.abs(skewness)>1e-31:
                                        p0s = [y1_list[i],(qty_values[0]-qty_values[-1])/3,np.max(hist)]
                                        p0s1 = [y2_list[i],(qty_values[0]-qty_values[-1])/3,np.max(hist)]
                                        par1,sig1 = curve_fit(two_gaussian,qty_values,hist,p0=[*p0s,*p0s1])
                                        
                                        ax_hist.plot(qty_values,two_gaussian(qty_values,*par1),'b-')
                                        ax_hist.plot(qty_values,gaussian(qty_values,*par1[:3]),'b--')
                                        ax_hist.plot(qty_values,gaussian(qty_values,*par1[3:]),'b--')
            
                                    
                                    ax_hist.plot(qty_values,gaussian(qty_values,*par),'r-')
                                    ax_hist.set_title(f'chi: {np.sum(chi) :.2e}'+
                                                      # f'\n mean fit - hist: {par[0] - np.mean(qty):.2e}'+
                                                      f'\n skew fit - hist {skewness:.2e}')
                                    fig_hist.canvas.draw()
                                
        
                            except Exception as e:
                                print(e)
                                
                            if plot or plot_dataset:
                                while not plt.waitforbuttonpress():
                                    pass
                        """
                        Binarisation
                        """
                        
                        N = len(qty)
                        qty0 = 0.5*(np.max(qty)+np.min(qty)) 
                        
                        if old_165_B and letter == 'B' and temp == '1.65':
                            valid_data = np.logical_or((dy1_list[i]>=qty-y1_list[i]) == (qty-y1_list[i] >= -dy1_list[i]),
                                                        (dy1_list[i]>=qty-y2_list[i]) == (qty-y2_list[i] >= -dy1_list[i]) )
                            qty = qty[valid_data]
                            population = np.where(qty>np.mean(np.array([y1_list[i],y2_list[i]])),np.ones_like(qty),np.zeros_like(qty))
                        
                            valid_data = np.logical_or((dy1_list[i]>=qty_f-y1_list[i]) == (qty_f-y1_list[i] >= -dy1_list[i]),
                                                        (dy1_list[i]>=qty_f-y2_list[i]) == (qty_f-y2_list[i] >= -dy1_list[i]) )
                            qty_f = qty_f[valid_data]
                            population_f = np.where(qty_f>np.mean(np.array([y1_list[i],y2_list[i]])),np.ones_like(qty_f),np.zeros_like(qty_f))
                        
                            
                        else:
                            try:
                                
                                
                                # mean = np.mean(qty)
                                # stdev = np.std(qty)
                                # if qty[0]>0.5*(np.max(qty) + np.min(qty)):
                                #     population = [True]
                                # else:
                                #     population = [False]
                                
                                # for j in range(len(qty)-4):
                                #     # plt.waitforbuttonpress()
                                #     first = True
        
                                    
                                #     n=2
                                #     if population[j]:
                                #         if j<n+1:
                                #             mean_state = np.mean(qty[:j+1][population])
                                #         else:
                                #             mean_state = np.mean(qty[-n+j+1:j+1][population[:-n-1:-1]])
                                        
                                #         if qty[j+4]-mean_state<-3*stdev:
                                #             population.append(False)
                                #         else:
                                #             population.append(True)
                                #     else:
                                #         if j<n+1:
                                #             mean_state = np.mean(qty[:j+1][not population])
                                #         else:
                                #             mean_state = np.mean(qty[-n+j+1:j+1][not population[:-n-1:-1]])
                                #         if qty[j+4]-mean_state>3*stdev:
                                #             population.append(True)
                                #         else:
                                #             population.append(False)
                                       
                                population = np.where(qty>=y1_list[i],np.ones_like(qty),np.zeros_like(qty))
                                population_f = np.where(qty_f>=y1_list[i],np.ones_like(qty_f),np.zeros_like(qty_f))
                            except Exception as e:
                                print(e)
                                population = np.where(qty>qty0,np.ones_like(qty),np.zeros_like(qty))
                                population_f = np.where(qty_f>qty0,np.ones_like(qty),np.zeros_like(qty))
        
                        """
                        End of binarisation
                        """
                        
                        higher = []
                        lower = []
                        last = -1
                        for pop in population:
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
                        
                        if len(higher)==0:
                            higher = [0]
                        if len(lower) ==0:
                            lower = [0]
                
                        if A==amp_last:
                            higher_stack = np.concatenate([np.array(higher),np.array(higher_stack)])
                            lower_stack = np.concatenate([np.array(lower),np.array(lower_stack)])
                
                        if (A != amp_last or file == files[-1]) and i >0:
                            amplitudes.append(amp_last/7)
                            y1s.append(y1_last)
                            y2s.append(y2_last)
                            mean_rs.append(r_last)
                            
                            lifetime_up = np.mean(higher_stack)/sample_rate
                            lifetime_down = np.mean(lower_stack)/sample_rate
                            lifetimes_up.append(lifetime_up)
                            lifetimes_down.append(lifetime_down)
                            
                            higher_stack = higher
                            lower_stack = lower
                            
                        
                        higher_f = []
                        lower_f = []
                        last_f = -1
                        for pop_f in population_f[:]:
                            if pop_f:
                                if last_f == 1:
                                    higher_f[-1]+=1
                                else:
                                    higher_f.append(1)
                            else:
                                if last_f == 0:
                                    lower_f[-1]+=1
                                else:
                                    lower_f.append(1)
                            last_f = pop_f
                        if len(higher_f)==0:
                            higher_f = [0]
                        if len(lower_f) ==0:
                            lower_f = [0]                    
                    
                        if A==amp_last_f:
                            higher_stack_f = np.concatenate([np.array(higher_f),np.array(higher_stack_f)])
                            lower_stack_f = np.concatenate([np.array(lower_f),np.array(lower_stack_f)])
                
                        if (A != amp_last_f or file == files[-1]) and i >0:
                            amplitudes_f.append(amp_last_f/3.5)
                            y1s.append(y1_last)
                            y2s.append(y2_last)
                            mean_rs.append(r_last)
                            
                            lifetime_up_f = np.mean(higher_stack_f)/sample_rate
                            lifetime_down_f = np.mean(lower_stack_f)/sample_rate
                            lifetimes_up_f.append(lifetime_up_f)
                            lifetimes_down_f.append(lifetime_down_f)
                            
                            higher_stack_f = higher_f
                            lower_stack_f = lower_f
                        
                        amp_last_f = A    
                        amp_last = A
                        y1_last = y1_list[i]
                        y2_last = y2_list[i]
                        r_last=np.mean(r)
        
                        # auto_up = sig.correlate(population,population)/np.std(population)**2
                        # x_ax = sig.correlation_lags(len(population),len(population))
                        
                        
                        # ax_corr.plot(x_ax,auto_up,c = cmap(amp_last/3.5))
                cmap_down = plt.get_cmap('Reds')
                cmap_up = plt.get_cmap('inferno')                          
        
                crit_drive = 0
                crit_r_up = 0#1.642e-8
                crit_r_down = 0#1.183e-8
                ax_lifetime.semilogy(np.array(amplitudes),lifetimes_down, '.-', lw=1, color='k',label='down')
                ax_lifetime.semilogy(np.array(amplitudes),lifetimes_up, '.-', lw=1, color='b',label='up')
          
                # ax_lifetime.semilogy(np.array(amplitudes),lifetimes_down_f, 'o-', lw=1, color='tab:orange',label='down')
                # ax_lifetime.semilogy(np.array(amplitudes),lifetimes_up_f, 'o-', lw=1, color='b',label='up')
                
        
                # ax_lifetime_log.loglog(np.abs(np.array(amplitudes)-crit_r_down),lifetimes_down, 'o-', lw=1, color=cmap((cutoff-cutoffs[0])/(cutoffs[-1]-cutoffs[0])),label='down')
                # ax_lifetime_log.loglog(np.abs(np.array(amplitudes)-crit_r_up),lifetimes_up, 'o-', lw=1, color=cmap((cutoff-cutoffs[0])/(cutoffs[-1]-cutoffs[0])),label='up')
                ax_lifetime.set_xlabel("Drive (a.u.)")
                ax_lifetime.set_ylabel("Mean state lifetime (s)")
                ax_lifetime.set_title(f'MEAN LIFETIME\n{letter}{distance}|folder: {folder}|type: buffer|tc: {tc}s')
                # ax_lifetime.legend()
                # fig_lifetime.tight_layout()
                # while not plt.waitforbuttonpress():
                #     pass
        
        

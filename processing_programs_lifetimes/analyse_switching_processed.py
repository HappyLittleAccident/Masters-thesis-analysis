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
from scipy import interpolate
import os

old_165_B=False
save = False
plot = False
plot_dataset = False
apply_notch = False

def gaussian(x,x0,w,A):
    return A*np.exp(-0.5*((x-x0)/w)**2)/(w*np.sqrt(2*np.pi))

def two_gaussian(x,x0,w,A,x01,w1,A1):
    return gaussian(x,x0,w,A) + gaussian(x,x01,w1,A1)

filter_order = 4
plt.close('all')

switching_sweeps_names = {'1.35':'T1350','1.45':'T1450','1.65':'T1650','1.85':'T1850'}
switching_ampsweeps_names = {'1.35':'T1350','1.45':'T1450','1.65':'T1650','1.85':'T1850'}

temps = ['1.35','1.45','1.65','1.85']
letters = ['A','B','C','D']
if plot:
    fig_hist,ax_hist = plt.subplots()


if plot_dataset:
    fig,ax = plt.subplots()
    fig_corr, ax_corr = plt.subplots()  
cmap = plt.get_cmap('viridis')
cutoffs = np.linspace(25,200,2)
for letter in letters:
    print(f'Resonator: {letter}')
    if letter == 'D':
        for temp in temps[0:]:
            
            fig_ampsweeps,ax_ampsweeps = plt.subplots(2,1,sharex=True)
            fig_lifetime, ax_lifetime = plt.subplots()
            # fig_lifetime_log,ax_lifetime_log = plt.subplots()
        
            parameters_to_save = {}
            print(f'Temperature: {temp}')
            date = '*'
            distance = 500
               
            folder = switching_sweeps_names[temp]
            folder_ampsweeps = switching_ampsweeps_names[temp]
            
            basedir = r'2023_4_Helmholtz_resonators_recal\Helmholtz_res_drive_dep\{}\buffer'.format(folder)
            
            ampsweeps_base = r'2023_4_Helmholtz_resonators_recal\Helmholtz_res_drive_dep\{}\ampsweeps_fast'.format(folder_ampsweeps)
            
            resdir = path.join(basedir,r'resonator_{}{}'.format(distance,letter))
            files = glob(path.join(resdir, r'*.npy'))
            
            ampsweeps = path.join(ampsweeps_base,r'resonator_{}{}_updowndown'.format(distance,letter))
            files_ampsweeps = glob(path.join(ampsweeps, r'*.npy'))
            
            # files = files[:-53]
            files.sort(key=lambda f: np.load(f, allow_pickle=True).item()['App_gen (V)'])

            
            params = np.load(f'lifetimes_params_processed\params_{letter}_T{temp}_3.npy', allow_pickle=True).item()
            y1_list = params['y1']
            y2_list = params['y2']
            dy1_list = params['dy1']
            dy2_list = params['dy2']
            p_points = params['p points']
            v_points = params['v points']
            
            # def interpolated(x):
            #     x_points = p_points
            #     y_points = v_points
            
            #     tck = interpolate.splrep(x_points, y_points,k=1)
            #     return interpolate.splev(x, tck)
            
            def interpolated(x):
                x_points = p_points
                y_points = v_points
            
                value = np.interp(x,x_points,y_points)
                return value
            
            
            if temp=='1.65' and letter == 'B' and old_165_B:
                params = np.load(f'params_{letter}_T{temp}_manual.npy', allow_pickle=True).item()
                y1_list = params['y1']
                y2_list = params['y2']
                dy1_list = params['dy']
    
    
            
            #%% analyse switching
            amplitudes = []
            lifetimes_up=[]
            lifetimes_down=[]
            
            series_up = []
            series_down = []
            
            amplitudes_f = []
            lifetimes_up_f=[]
            lifetimes_down_f=[]
            
            higher_stack = []
            lower_stack = []
            
            higher_stack_f = []
            lower_stack_f = []
            mean_rs = []
            mean_rs_down = []
            mean_rs_up = []
            
            amp_last = -1
            P_last = -1
            pressures = []
            y2_last = 0
            y1_last = 0
            r_mean_last = 0
            r_down_last=0
            r_up_last = 0
            for i,file in enumerate(files): #and files[39:]:
                # print(file)
                d = np.load(file, allow_pickle=True).item()
            
                r = np.array(d['Velocity (m/s)'])
                P = np.array(d['Pressure (Pa/m)'])
                t = d['time (s)']
                A = np.array(d['App_gen (V)'])
                sample_rate = d['sample rate (Hz)']
                # print(f'Date: {file[-27:-8]} \t Drive amplitude: {A:.4f}Vrms')
                
                
                # if apply_notch:
                #     transformed = np.abs(fft.fft(r))
                #     freq = fft.fftfreq(len(r),1/sample_rate)
                    
                #     indices,_s = signal.find_peaks(transformed,prominence=0.2)
                #     peaks = transformed[indices]
    
                #     freq_select = freq[indices]
                #     dummy_peaks = np.copy(peaks)
                #     r_b = r
                #     for j in range(len(peaks)-1):
                #         if np.abs(freq_select[j])>40 and np.abs(freq_select[j]-freq_select[i+1]) > 8 or j==len(peaks)-2:
                            
                #             f_out= freq_select[np.argmax(dummy_peaks[:i+1])]
                    
                #             dummy_peaks[:j+1]*=0
                            
                #             b_n,a_n = signal.iirnotch(np.abs(f_out),40,sample_rate)
                            
                #             r_b = np.copy(signal.filtfilt(b_n, a_n,r_b))
                #     qty = r_b
                # else:
                qty=r
                # if 0:
                #     hist,edges = np.histogram(qty,80)
                    
                #     qty_values=0.5*(edges[1:]+edges[:-1])
                #     if plot:
                #         ax_hist.clear()
                #         ax_hist.plot(qty_values,hist,'.-k',ms=1)
                #     if plot_dataset:
                #         ax.clear()
                #         ax.plot(t,r,'.-k',ms=2,lw=0.2)
                #         ax.axhline(np.mean(r),c='r',ls='-')
                #         ax.axhline(np.mean(r)+np.std(r),c='r',ls='--')
                #         ax.axhline(np.mean(r)-np.std(r),c='r',ls='--')
                #         # try:
                #         #     ax.axhline(y1_list[i],c='r',ls='-')
                #         #     ax.axhline(y1_list[i]+3*dy1_list[i],c='r',ls='--')
                #         #     ax.axhline(y1_list[i]-3*dy1_list[i],c='r',ls='--')
                #         #     ax.axhline(y2_list[i]+dy2_list[i],c='b',ls='--')
                #         #     ax.axhline(y2_list[i]-dy2_list[i],c='b',ls='--')
                #         #     ax.axhline(y2_list[i],c='b',ls='-')
                #         # except:
                #             # pass
                #         fig.canvas.draw()
                    
                            
                #     """
                #     Fit gaussian on data
                #     """
                #     try:
                #         if plot:
                #             p0s = [y1_list[i],(qty_values[0]-qty_values[-1])/3,np.max(hist)]
                #             par,sig = curve_fit(gaussian,qty_values,hist,p0=p0s)
                        
                #             chi = np.sum((gaussian(qty_values,*par)/(hist+1) -1)**2)/len(hist)
                #             skewness= np.mean((qty-np.mean(qty))**3)
                #             if np.abs(skewness)>1e-31:
                #                 p0s = [y1_list[i],(qty_values[0]-qty_values[-1])/3,np.max(hist)]
                #                 p0s1 = [y2_list[i],(qty_values[0]-qty_values[-1])/3,np.max(hist)]
                #                 par1,sig1 = curve_fit(two_gaussian,qty_values,hist,p0=[*p0s,*p0s1])
                                
                #                 ax_hist.plot(qty_values,two_gaussian(qty_values,*par1),'b-')
                #                 ax_hist.plot(qty_values,gaussian(qty_values,*par1[:3]),'b--')
                #                 ax_hist.plot(qty_values,gaussian(qty_values,*par1[3:]),'b--')
    
                            
                #             ax_hist.plot(qty_values,gaussian(qty_values,*par),'r-')
                #             ax_hist.set_title(f'chi: {np.sum(chi) :.2e}'+
                #                               # f'\n mean fit - hist: {par[0] - np.mean(qty):.2e}'+
                #                               f'\n skew fit - hist {skewness:.2e}')
                #             fig_hist.canvas.draw()
                        

                #     except Exception as e:
                #         print(e)
                        
                #     if plot or plot_dataset:
                #         while not plt.waitforbuttonpress():
                #             pass
                """
                Binarisation
                """
                
                N = len(qty)
                qty0 = 0.5*(np.max(qty)+np.min(qty)) 
                
                if old_165_B and letter == 'B' and temp == '1.65':
                    valid_data = np.logical_or((dy1_list[i]>=qty-y1_list[i]) == (qty-y1_list[i] >= -dy1_list[i]),
                                                (dy1_list[i]>=qty-y2_list[i]) == (qty-y2_list[i] >= -dy1_list[i]) )
                    print(valid_data)
                    qty = qty[valid_data]
                    population = np.where(qty>np.mean(np.array([y1_list[i],y2_list[i]])),np.ones_like(qty),np.zeros_like(qty))
            
                    
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
                               
                        population = np.where(qty>=interpolated(P),np.ones_like(qty),np.zeros_like(qty))
                    except Exception as e:
                        print(e)
                        population = np.where(qty>qty0,np.ones_like(qty),np.zeros_like(qty))

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
                    pressures.append(P_last)
                    amplitudes.append(amp_last)

                    mean_rs.append(r_mean_last)
                    mean_rs_up.append(r_up_last)
                    mean_rs_down.append(r_down_last)
                    
                    lifetime_up = np.mean(higher_stack)/sample_rate
                    lifetime_down = np.mean(lower_stack)/sample_rate
                    lifetimes_up.append(lifetime_up)
                    lifetimes_down.append(lifetime_down)
                    
                    series_up.append(np.array(higher_stack)/sample_rate)
                    series_down.append(np.array(lower_stack)/sample_rate)
                    
                    higher_stack = higher
                    lower_stack = lower
                    
                P_last = P
                amp_last = A
                
                # y1_last = y1_list[i]
                # y2_last = y2_list[i]
                r_mean_last = np.mean(r)
                r_down_last=np.mean(r[population==0])
                r_up_last = np.mean(r[population==1])

                # auto_up = sig.correlate(population,population)/np.std(population)**2
                # x_ax = sig.correlation_lags(len(population),len(population))
                
                
                # ax_corr.plot(x_ax,auto_up,c = cmap(amp_last/3.5))
            cmap_down = plt.get_cmap('Reds')
            cmap_up = plt.get_cmap('inferno')                          
                
            
            
            crit_drive = 0
            crit_r_up = 0#1.642e-8
            crit_r_down = 0#1.183e-8
            ax_lifetime.semilogy(mean_rs_down,lifetimes_down, 'r.-', lw=1,label='down')
            ax_lifetime.semilogy(mean_rs_up,lifetimes_up, 'k.-', lw=1,label='up')
      
            # ax_lifetime.semilogy(np.array(amplitudes),lifetimes_down_f, '.-', lw=1, color='tab:orange',label='down')
            # ax_lifetime.semilogy(np.array(amplitudes),lifetimes_up_f, '.-', lw=1, color='b',label='up')
            
    
            # ax_lifetime_log.loglog(np.abs(np.array(amplitudes)-crit_r_down),lifetimes_down, 'o-', lw=1, color=cmap((cutoff-cutoffs[0])/(cutoffs[-1]-cutoffs[0])),label='down')
            # ax_lifetime_log.loglog(np.abs(np.array(amplitudes)-crit_r_up),lifetimes_up, 'o-', lw=1, color=cmap((cutoff-cutoffs[0])/(cutoffs[-1]-cutoffs[0])),label='up')
            ax_lifetime.set_xlabel("Drive (a.u.)")
            ax_lifetime.set_ylabel("Mean state lifetime (s)")
            ax_lifetime.set_title(f'MEAN LIFETIME\n{letter}{distance}|folder: {folder}|type: buffer')
            # ax_lifetime.legend()
            # fig_lifetime.tight_layout()
            # while not plt.waitforbuttonpress():
            #     pass
            d = np.load(files_ampsweeps[0], allow_pickle=True).item()
            
                    
            upP = d['upPressure (Pa/m)']
            downP = d['downPressure (Pa/m)']
            jumpP = d['jumpPressure (Pa/m)']
            upV=d['upVelocity (m/s)']
            downV=d['downVelocity (m/s)']
            jumpV=d[ 'jumpVelocity (m/s)']
            
            ax_ampsweeps[0].plot(upP, upV, '-o',ms=0.5, lw=0.5, color='tab:blue',alpha=1)
            ax_ampsweeps[0].plot(downP,downV, '-o',ms=0.5, lw=0.5, color='tab:orange',alpha=1)
            ax_ampsweeps[0].plot(jumpP, jumpV, '-o',ms=0.5, lw=0.5, color='tab:green',alpha=1)
            ax_ampsweeps[1].semilogy(np.array(pressures),lifetimes_down, 'r.-', lw=1,label='down')
            
            
            if save:
                if apply_notch:
                    notch = '_notch_filter'
                else:
                    notch = ''
                parameters_to_save={'series lifetimes up (s)':series_up,'series lifetimes down (s)':series_down,
                                    'mean lifetimes up (s)':np.array(lifetimes_up),
                                    'mean lifetimes down (s)':np.array(lifetimes_down),
                                    'P (Pa/m)':np.array(pressures),'sample rate (Hz)':sample_rate,
                                    'App_gen (V)':np.array(amplitudes),'V mean (m/s)':np.array(mean_rs),
                                    'V up (m/s)':np.array(mean_rs_up),'V down (m/s)':np.array(mean_rs_down)}
                path_save = r'..\lifetimes'
                os.makedirs(path_save,exist_ok=True)
                path_save = os.path.join(path_save,f'resonator_500{letter}_{temp}{notch}_3.npy')
                np.save(path_save, parameters_to_save)

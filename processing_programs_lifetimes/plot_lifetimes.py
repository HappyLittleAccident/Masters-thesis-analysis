# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 19:47:08 2023

@author: Marek
"""

import numpy as np
import matplotlib.pyplot as plt
from os import path
from glob import glob

plt.close('all')

plot_histogram= False

switching_ampsweeps_names = {'1.35':'T1350','1.45':'T1450','1.65':'T1650','1.85':'T1850'}

temps = ['1.35','1.45','1.65','1.85']
letters = ['A','B','C','D']

fig_hist,ax_hist = plt.subplots()
fig_cdf,ax_cdf = plt.subplots()

for letter in letters:
    print(f'Resonator: {letter}')
    if letter == 'B':
        for temp in temps[1:4]:
            fig,ax = plt.subplots()
            fig_lifetime, ax_lifetime = plt.subplots()
            for spline in ['piecewise']:
                # fig_ampsweeps,ax_ampsweeps = plt.subplots(2,1,sharex=True)
                
    
                # fig_lifetime_log,ax_lifetime_log = plt.subplots()
            
                parameters_to_save = {}
                print(f'Temperature: {temp}')
                date = '*'
                distance = 500
                   
                folder_ampsweeps = switching_ampsweeps_names[temp]
                
                ampsweeps_base = r'..\2023_4_Helmholtz_resonators_recal\Helmholtz_res_drive_dep\{}\ampsweeps_fast'.format(folder_ampsweeps)
                
                
                lifetimes_parameters = fr'D:\OneDrive_copy\OneDrive - Univerzita Karlova\DATA\lifetimes\{spline}_spline'
                files = glob(path.join(lifetimes_parameters, r'*.npy'))
                
                probabilities_parameters = r'D:\OneDrive_copy\OneDrive - Univerzita Karlova\DATA\lifetimes'
                files_probabilities = glob(path.join(probabilities_parameters,r'*.npy'))
                
                ampsweeps = path.join(ampsweeps_base,r'resonator_{}{}_updowndown'.format(distance,letter))
                files_ampsweeps = glob(path.join(ampsweeps, r'*.npy'))
                
                
        
                files_t = [file for file in files if temp in file]
                files_ampsweeps_t = [file for file in files_probabilities if switching_ampsweeps_names[temp] in file]
                for i,file in enumerate(files_t): #and files[39:]:
                    # print(file)
                    d = np.load(file, allow_pickle=True).item()
                    V_down = d['V down (m/s)']
                    V_up = np.array(d['V up (m/s)'])
                    V_mean = d['V mean (m/s)']
                    P = np.array(d['P (Pa/m)'])
                    A = np.array(d['App_gen (V)'])
                    lifetimes_up = d['mean lifetimes up (s)']
                    lifetimes_down = d['mean lifetimes down (s)']
                    sample_rate = d['sample rate (Hz)']
                    populations_down = d['series lifetimes down (s)']
                    populations_up = d['series lifetimes up (s)']
                    # print(f'Date: {file[-27:-8]} \t Drive amplitude: {A:.4f}Vrms')
                    
                    if not int(file.split('_')[-1][-5]) == 2 or 'notch' in file:
                        continue
                    
                    try:
                         d = np.load(files_ampsweeps_t[i],allow_pickle=True).item()
                         prob_up = d['Prob up']
                         prob_down = d['Prob down']
                         P_amp = d['P (Pa/m)']
                         V_amp_down = d['V down (m/s)']
                         V_amp_up = d['V up (m/s)']
                         ax.semilogy(V_amp_down,prob_down,'r.-',label='down ampsweep')
                         ax.semilogy(V_amp_up,prob_up,'g.-',label='up ampsweep')
                         ax.set_xlabel('Pressure (Pa/m)')
                         ax.set_ylabel('Prob of finding system in state')
                         ax.set_title(f'Temperature: {temp}K')
                         
                    except Exception as e:
                        print(e)
                    cmap_down = plt.get_cmap('Reds')
                    cmap_up = plt.get_cmap('inferno')                          
                        
                    
                    
                    crit_drive = 0
                    crit_r_up = 0#1.642e-8
                    crit_r_down = 0#1.183e-8
                    
                    show_notch = True
                    if spline == 'lin':
                        c='r'
                    elif spline == 'cubic':
                        c='b'
                    else:
                        c='k'
                    lab = ''
                    if 'notch' in file:
                        lab = ' notch'
                        c='k'
                        show_notch = False
                    
                    lab += '_'+spline+'_'+file.split('_')[-1][-5]

                    
                    T_mean_up_array = []
                    T_ups = []
                    T_downs = []
                    if show_notch:
                        # ax_lifetime.semilogy(P,lifetimes_down, '.-'+c, lw=1,label='down' + lab)
                        ax_lifetime.semilogy(P,lifetimes_up, '.-', lw=1,label='up')
                        First_nonzero = False
                        for i,T_mean_down in enumerate(lifetimes_down):
                            T_up= np.sum(populations_up[i])
                            T_down=np.sum(populations_down[i])
                            T_mean_up = T_mean_down * T_up/T_down 
                            T_mean_up_array.append(T_mean_up)
                            T_ups.append(T_up)
                            T_downs.append(T_down) 
                            if not First_nonzero and T_mean_up == 0:
                                First_nonzero = True
                                P0 = P[i]
                                V0 = V_down[i]
                            
                        ax_lifetime.semilogy(P,T_mean_up_array,'.-m')
                    
                        T_ups = np.array(T_ups)
                        T_downs = np.array(T_downs)
                        
                         
                        ax.semilogy(V_down,T_downs/(T_ups + T_downs),'.-k',label='down')
                        
                        ax.semilogy(V_up,T_ups/(T_ups + T_downs),'.-b',label='up')
                        ax.legend()
                    
            
                    # ax_lifetime_log.loglog(np.abs(np.array(amplitudes)-crit_r_down),lifetimes_down, 'o-', lw=1, color=cmap((cutoff-cutoffs[0])/(cutoffs[-1]-cutoffs[0])),label='down')
                    # ax_lifetime_log.loglog(np.abs(np.array(amplitudes)-crit_r_up),lifetimes_up, 'o-', lw=1, color=cmap((cutoff-cutoffs[0])/(cutoffs[-1]-cutoffs[0])),label='up')
                    ax_lifetime.set_xlabel("State velocity (m/s)")
                    ax_lifetime.set_ylabel("Mean state lifetime (s)")
                    ax_lifetime.set_title(f'MEAN LIFETIME\n{letter}{distance}|type: buffer|T: {temp}')
                    ax_lifetime.legend()
                    
                    fig_lifetime.savefig(f'{temp}.png')
                    
                    if plot_histogram: #and lab == ' notch_1':
                        
                        
                        Ts_hists = []
                        max_Ts= 2#np.max(np.concatenate(populations_down))
                        min_Ts = np.min(np.concatenate(populations_down))
                        for times in populations_down:
                            
                            hist,edges = np.histogram(times,bins=np.linspace(min_Ts,max_Ts,200))
                            Ts_values=0.5*(edges[1:]+edges[:-1])
                            Ts_hists.append(hist)
                        
                    
                    
                        Ps,Ts = np.meshgrid(P,Ts_values)
                        # ax_hist.pcolormesh(Ps,Ts,np.array(Ts_hists).T,cmap = 'viridis',shading = 'nearest',norm='log')
                        
                        count_index = 0
                        for pop,pop1 in zip(populations_down,populations_up):
                            ax_hist.clear()
                            ax_hist.hist(pop,20,label='down',weights=np.ones_like(pop)/len(pop))
                            ax_hist.hist(pop1,20,alpha=0.5,label='up',weights=np.ones_like(pop1)/len(pop1))
                            ax_hist.set_yscale('log')
                            ax_hist.set_xlabel('State lifetime (s)')
                            ax_hist.set_ylabel('Lifetime probability')
                            ax_hist.set_title(f'Temperature {temp}K \nPressure: {1e-3*P[count_index]:.2f}kPa/m')
                            ax_hist.legend()
                            fig_hist.canvas.draw()
                            
                            ax_cdf.clear()
                            cdf_up = (np.arange(len(pop1))+1)/len(pop1)
                            cdf_down = (np.arange(len(pop))+1)/len(pop)
                            # ax_cdf.hist(pop,20,cumulative=True)
                            # ax_cdf.hist(pop1,20,alpha=0.5,cumulative=True)
                            ax_cdf.plot(np.sort(pop),1-cdf_down,'-',label='down')
                            ax_cdf.plot(np.sort(pop1),1-cdf_up,'-',label='up')
                            ax_cdf.set_yscale('log')
                            ax_cdf.set_xlabel('State lifetime (s)')
                            ax_cdf.set_ylabel('1 - CDF')
                            ax_cdf.set_title(f'Temperature {temp}K\nPressure: {1e-3*P[count_index]:.2f}kPa/m')
                            ax_cdf.legend()
                            fig_cdf.canvas.draw()
                            count_index +=1
                            
                            while not plt.waitforbuttonpress():
                                pass
                    
                    
                    # fig_lifetime.tight_layout()
                    # while not plt.waitforbuttonpress():
                    #     pass
                    # d = np.load(files_ampsweeps[0], allow_pickle=True).item()
                    
                            
                    # upP = d['upPressure (Pa/m)']
                    # downP = d['downPressure (Pa/m)']
                    # jumpP = d['jumpPressure (Pa/m)']
                    # upV=d['upVelocity (m/s)']
                    # downV=d['downVelocity (m/s)']
                    # jumpV=d[ 'jumpVelocity (m/s)']
                    
                    # ax_ampsweeps[0].plot(upP, upV, '-o',ms=0.5, lw=0.5, color='tab:blue',alpha=1)
                    # ax_ampsweeps[0].plot(downP,downV, '-o',ms=0.5, lw=0.5, color='tab:orange',alpha=1)
                    # ax_ampsweeps[0].plot(jumpP, jumpV, '-o',ms=0.5, lw=0.5, color='tab:green',alpha=1)
                    # ax_ampsweeps[1].semilogy(P,lifetimes_down, 'r.-', lw=1,label='down')
    
    

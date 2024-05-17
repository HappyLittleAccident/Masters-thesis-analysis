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

save = True
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
    if letter == 'B':
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
            
            """
            ########################################################################################################
            Given temp and letter, get:
                filename of file with selected points
                directory of buffer data
                directory of correspoding amspweep data
            
            Load filenames and sort them               
            """
            params_name = f'lifetimes_energylevels\params_{letter}_T{temp}.npy'
            basedir = r'..\..\2023_4_Helmholtz_resonators_recal_2\Helmholtz_buffer_measurements\{}\buffer'.format(folder)
            ampsweeps_base = r'..\..\2023_4_Helmholtz_resonators_recal_2\Helmholtz_res_drive_dep\{}\ampsweeps_fast'.format(folder_ampsweeps)            

            resdir = path.join(basedir,r'resonator_{}{}'.format(distance,letter))
            files = glob(path.join(resdir, r'*.npy'))
            files.sort(key=lambda f: np.load(f, allow_pickle=True).item()['App_gen (V)'])

            ampsweeps = path.join(ampsweeps_base,r'resonator_{}{}_updowndown'.format(distance,letter))
            files_ampsweeps = glob(path.join(ampsweeps, r'*.npy'))
            """
            ########################################################################################################
            Load selected points and make an interpolated curve
            """
            params = np.load(params_name, allow_pickle=True).item()
            p_points = params['p points']
            v_points = params['v points']
            
            def interpolated(x):
                x_points = p_points
                y_points = v_points
            
                value = np.interp(x,x_points,y_points)
                return value            


            
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

            r_mean_last = 0
            r_down_last = 0
            r_up_last   = 0
            for i,file in enumerate(files): #and files[39:]:
                # print(file)
                d = np.load(file, allow_pickle=True).item()
            
                r = np.array(d['Velocity (m/s)'])
                P = np.array(d['Pressure (Pa/m)'])
                t = d['time (s)']
                A = np.array(d['App_gen (V)'])
                sample_rate = d['sample rate (Hz)']
                qty=r

                """
                Binarisation
                """
                
                N = len(qty)
                qty0 = 0.5*(np.max(qty)+np.min(qty)) 
                    
                try:                           
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
                

                r_mean_last = np.mean(r)
                r_down_last=np.mean(r[population==0])
                r_up_last = np.mean(r[population==1])

            cmap_down = plt.get_cmap('Reds')
            cmap_up = plt.get_cmap('inferno')                          
                
            ax_lifetime.semilogy(mean_rs_down,lifetimes_down, 'r.-', lw=1,label='down')
            ax_lifetime.semilogy(mean_rs_up,lifetimes_up, 'k.-', lw=1,label='up')
      
            # ax_lifetime.semilogy(np.array(amplitudes),lifetimes_down_f, '.-', lw=1, color='tab:orange',label='down')
            # ax_lifetime.semilogy(np.array(amplitudes),lifetimes_up_f, '.-', lw=1, color='b',label='up')
            
            ax_lifetime.set_xlabel("Drive (a.u.)")
            ax_lifetime.set_ylabel("Mean state lifetime (s)")
            ax_lifetime.set_title(f'MEAN LIFETIME\n{letter}{distance}|folder: {folder}|type: buffer')
            ax_lifetime.legend()
            fig_lifetime.tight_layout()

            
            d_amp = np.load(files_ampsweeps[0], allow_pickle=True).item()
                    
            upP = d_amp['upPressure (Pa/m)']
            downP = d_amp['downPressure (Pa/m)']
            jumpP = d_amp['jumpPressure (Pa/m)']
            upV=d_amp['upVelocity (m/s)']
            downV=d_amp['downVelocity (m/s)']
            jumpV=d_amp[ 'jumpVelocity (m/s)']
            
            ax_ampsweeps[0].plot(upP, upV, '-o',ms=0.5, lw=0.5, color='tab:blue',alpha=1)
            ax_ampsweeps[0].plot(downP,downV, '-o',ms=0.5, lw=0.5, color='tab:orange',alpha=1)
            ax_ampsweeps[0].plot(jumpP, jumpV, '-o',ms=0.5, lw=0.5, color='tab:green',alpha=1)
            ax_ampsweeps[1].semilogy(np.array(pressures),lifetimes_down, 'r.-', lw=1,label='down')
            

            if save:
                parameters_to_save={'series lifetimes up (s)':series_up,'series lifetimes down (s)':series_down,
                                    'mean lifetimes up (s)':np.array(lifetimes_up),
                                    'mean lifetimes down (s)':np.array(lifetimes_down),
                                    'P (Pa/m)':np.array(pressures),'sample rate (Hz)':sample_rate,
                                    'App_gen (V)':np.array(amplitudes),'V mean (m/s)':np.array(mean_rs),
                                    'V up (m/s)':np.array(mean_rs_up),'V down (m/s)':np.array(mean_rs_down)}
                path_save = r'lifetimes_processed'
                os.makedirs(path_save,exist_ok=True)
                path_save = os.path.join(path_save,f'resonator_500{letter}_{temp}.npy')
                np.save(path_save, parameters_to_save)

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

show_hist = True
save  = False
plt.close('all')

switching_sweeps_names = {'1.35':'T1350_3','1.45':'T1450_3','1.65':'T1650_3','1.85':'T1850_2'}
switching_ampsweeps_names = {'1.35':'T1350','1.45':'T1450','1.65':'T1650','1.85':'T1850'}

temps = ['1.35','1.45','1.65','1.85']
letters = ['A','B','C','D']


cmap = plt.get_cmap('inferno')
cmap_up = plt.get_cmap('viridis')
cmap_down = plt.get_cmap('viridis')
cutoffs = np.linspace(25,200,2)
for letter in letters:
    print(f'Resonator: {letter}')
    if letter == 'B':
        ampsweeps_folders=glob(r'D:\OneDrive_copy\OneDrive - Univerzita Karlova\DATA\2023_4_Helmholtz_resonators_recal\Helmholtz_res_drive_dep\*')
        ampsweeps_temps = [x.split('\\')[-1][1:] for x in ampsweeps_folders if not '_' in x.split('\\')[-1]]
        fig_lifetime_up, ax_lifetime_up = plt.subplots()
        fig_lifetime_down, ax_lifetime_down = plt.subplots()
        v1s=[]
        v2s=[]
        P1s=[]
        P2s=[]
        stdevs_integrated= []
        for temp in ampsweeps_temps[0:]:

            
            # fig_lifetime,ax_lifetime = plt.subplots()
        
            parameters_to_save = {}
            print(f'Temperature: {temp}')
            distance = 500
            
            
            folder_ampsweeps = fr'D:\OneDrive_copy\OneDrive - Univerzita Karlova\DATA\2023_4_Helmholtz_resonators_recal\Helmholtz_res_drive_dep\T{temp}'
            files_ampsweeps = glob(folder_ampsweeps + fr'\ampsweeps_fast\resonator_{distance}{letter}_updowndown\*.npy')
            
            params = np.load(f'lifetimes_params_ampsweeps_processed\params_{letter}_TT{temp}_2.npy', allow_pickle=True).item()
            y1_list = params['y1']
            y2_list = params['y2']
            dy1_list = params['dy1']
            dy2_list = params['dy2']
            p_points = params['p points']
            v_points = params['v points']
            
            def interpolated(x):
                x_points = p_points
                y_points = v_points
            
                value = np.interp(x,x_points,y_points)
                return value
    
    
            
            #%% analyse switching
            vs = []
            ps = []
            mean_vs_down= []
            mean_vs_up=[]
            mean_vs=[]
            stdevs=[]
            
            # plt.figure()
            for file in files_ampsweeps: #and files[39:]:
                data = np.load(file, allow_pickle=True).item()
                
                upP = data['upPressure (Pa/m)']
                downP = data['downPressure (Pa/m)']
                jumpP = data['jumpPressure (Pa/m)']
                upV=data['upVelocity (m/s)']
                downV=data['downVelocity (m/s)']
                jumpV=data[ 'jumpVelocity (m/s)']
                
                # plt.cla()
                # plt.plot(downP,interpolated(downP),'r')
                # plt.plot(downP,downV,'g.')
                # plt.draw()
                # while not plt.waitforbuttonpress():
                #     pass
                
                # self.vs.append(upV)
                vs.append(downV)
                vs.append(jumpV)
                
                
                # self.ps.append(upP)
                ps.append(downP)
                ps.append(jumpP)
            
            """
            Put data in histogram
            """
            vs = np.concatenate(vs)
            ps = np.concatenate(ps)
            
            hist,p_bins,v_bins = np.histogram2d(ps,vs,bins = [500,500])    
            p_values = 0.5*(p_bins[1:]+p_bins[:-1])
            v_values = 0.5*(v_bins[1:]+v_bins[:-1])

            """
            calculate probabilities
            """
            up_count = np.zeros_like(p_values)
            down_count = np.zeros_like(p_values)
            

            for i,p in enumerate(p_values):
                up_state = v_values > interpolated(p)
                down_state = np.invert(up_state)
                
                
                up = np.sum(hist[i,up_state])
                down = np.sum(hist[i,down_state])
                both = up+down
                
                up_count[i]= up
                down_count[i]= down
                
                mean_vs_up.append(np.sum(hist[i,up_state]*v_values[up_state])/up)
                mean_vs_down.append(np.sum(hist[i,down_state]*v_values[down_state])/down)
                mean_vs.append(np.sum(hist[i]*v_values)/both)
                
                stdevs.append(np.sqrt(np.sum((hist[i,:]*v_values**2)/np.sum(hist[i,:])) - (np.sum(hist[i]*v_values)/np.sum(hist[i,:]))**2))
            
            
            both_count = up_count + down_count
            prob_up = up_count/both_count
            prob_down = down_count/both_count
            
            
            
            if temp == '1950':
                index_v2 = -1
                index_v1 = -1
            else:
                index_v2 = np.nonzero(prob_up==0)[0][0]
                index_v1 = np.nonzero(prob_down)[0][0]
                
                v1s.append(mean_vs[index_v1])
                v2s.append(mean_vs[index_v2])
                P1s.append(p_values[index_v1])
                P2s.append(p_values[index_v2])
                stdevs_integrated.append(np.sum(stdevs))
            
            cmap_down = plt.get_cmap('Reds')
            cmap_up = plt.get_cmap('inferno')   
                       
            """
            Show relevant stuff on 3D hist
            """
            if show_hist:
                mean_vs= np.array(mean_vs)
                mean_vs_up= np.array(mean_vs_up)
                mean_vs_down= np.array(mean_vs_down)
                fig_hist3D,ax_hist3D = plt.subplots()
                P,V = np.meshgrid(p_bins,v_bins)
                ax_hist3D.pcolormesh(P,V,hist.T,norm='log',cmap=cmap,alpha=0.8)
                # ax_hist3D.plot(p_values,interpolated(p_values),'b--',label='cut line')
                ax_hist3D.plot(p_values,mean_vs_up-mean_vs_down,'k-',label='V up state')
                # ax_hist3D.plot(p_values,mean_vs_down,'r-',label='V mean')
                # ax_hist3D.plot(p_values,mean_vs,'g-',label='V down state')
                ax_hist3D.legend()
                
                
            """
            Plot probabilities
            """
            crit_drive = 0
            crit_r_up = 0#1.642e-8
            crit_r_down = 0#1.183e-8
            ax_lifetime_down.semilogy(p_values,prob_down, '-',c=cmap_up((float(temp)-1325)/575), lw=1,label='down')
            ax_lifetime_up.semilogy(p_values,prob_up,'-',c=cmap_up((float(temp)-1325)/575), lw=1,label='up')

            # ax_lifetime.semilogy(p_values,prob_down, '-',c=cmap_up((float(temp)-1325)/575), lw=1,label='down')
            # ax_lifetime.semilogy(p_values,prob_up,'-',c=cmap_up((float(temp)-1325)/575), lw=1,label='up')
            # ax_lifetime.semilogy(p_values,stdevs)
            
            # plt.figure()
            # plt.plot(p_values,mean_vs_up,'k')
            # plt.plot(p_values,mean_vs_down,'r')
            # plt.plot(p_values,mean_vs,'b')
            # ax_lifetime.semilogy(np.array(amplitudes),lifetimes_down_f, '.-', lw=1, color='tab:orange',label='down')
            # ax_lifetime.semilogy(np.array(amplitudes),lifetimes_up_f, '.-', lw=1, color='b',label='up')
            
    
            # ax_lifetime_log.loglog(np.abs(np.array(amplitudes)-crit_r_down),lifetimes_down, 'o-', lw=1, color=cmap((cutoff-cutoffs[0])/(cutoffs[-1]-cutoffs[0])),label='down')
            # ax_lifetime_log.loglog(np.abs(np.array(amplitudes)-crit_r_up),lifetimes_up, 'o-', lw=1, color=cmap((cutoff-cutoffs[0])/(cutoffs[-1]-cutoffs[0])),label='up')
            ax_lifetime_up.set_xlabel("Drive (a.u.)")
            ax_lifetime_up.set_ylabel("State probability (s)")
            ax_lifetime_up.set_title(f'{letter}{distance}\n temp: {temp}, type: ampsweeps')
            
            ax_lifetime_down.set_xlabel("Drive (a.u.)")
            ax_lifetime_down.set_ylabel("State probability (s)")
            ax_lifetime_down.set_title(f'{letter}{distance}\n temp: {temp}, type: ampsweeps')
            # ax_lifetime.legend()
            # fig_lifetime.tight_layout()
            # while not plt.waitforbuttonpress():
            #     pass

            
            if save:
                parameters_to_save={'P (Pa/m)':np.array(p_values),'V mean (m/s)':np.array(mean_vs),'Prob up':prob_up,'Prob down':prob_down,
                                    'V up (m/s)':np.array(mean_vs_up),'V down (m/s)':np.array(mean_vs_down)}
                path_save = r'..\lifetimes'
                os.makedirs(path_save,exist_ok=True)
                path_save = os.path.join(path_save,f'resonator_500{letter}_{temp}_ampsweeps_2.npy')
                np.save(path_save, parameters_to_save)
        
        # fig_crit,ax_crit = plt.subplots()
        # ampsweeps_temps = [float(x)*1e-3 for x in ampsweeps_temps]
        # ax_crit.plot(ampsweeps_temps[:-1],stdevs_integrated,'r.-',label = 'Std integrated')
        # # ax_crit.plot(ampsweeps_temps[:-1],v2s,'k.-',label = 'Switching end')
        # ax_crit.set_xlabel('T (K)')
        # ax_crit.set_ylabel('V (m/s)')
        # ax_crit.legend()        
        
        # fig_critP,ax_critP = plt.subplots()
        # ax_critP.plot(ampsweeps_temps[:-1],-np.array(P1s)+np.array(P2s),'r.-',label = 'Switching duration')
        # # ax_critP.plot(ampsweeps_temps[:-1],P2s,'k.-',label = 'Switching end')
        # ax_critP.set_xlabel('T (K)')
        # ax_critP.set_ylabel('P (Pa/m)')
        # ax_critP.legend()
        
        
        

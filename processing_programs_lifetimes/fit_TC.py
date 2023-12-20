# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 12:08:21 2023

@author: Marek
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from glob import glob
from os import path


def up(x,T,A,x0,A0):
    return A*(1-np.exp(-(x-x0)/T)) +A0

def down(x,T,A,x0,A0):
    return A*np.exp(-(x-x0)/T) +A0


plt.close('all')

switching_sweeps_names = {'1.35':'T1350_3','1.45':'T1450_3','1.65':'T1650_3','1.85':'T1850_2'}
switching_ampsweeps_names = {'1.35':'T1350_2','1.45':'T1450_2','1.65':'T1650_2','1.85':'T1850'}

temps = ['1.35','1.45','1.65','1.85']
letters = ['A','B','C','D']


fig,ax = plt.subplots()
fig_fit,ax_fit = plt.subplots()

cmap = plt.get_cmap('viridis')
for letter in letters:
    print(f'Resonator: {letter}')
    if letter == 'D':
        for temp in temps[1:2]:
            print(f'Temperature: {temp}')
            date = '*'
            distance = 500
               
            folder = switching_sweeps_names[temp]
            folder_ampsweeps = switching_ampsweeps_names[temp]
            
            basedir = r'..\2023-03-Helmholtz_resonators\Helmholtz_res_drive_dep\{}\buffer'.format(folder)
            
            ampsweeps_base = r'D:..\2023-03-Helmholtz_resonators\Helmholtz_res_drive_dep\{}\ampsweeps_fast'.format(folder_ampsweeps)
            
            resdir = path.join(basedir,r'resonator_{}{}'.format(distance,letter))
            files = glob(path.join(resdir, r'*.npy'))
            
            ampsweeps = path.join(ampsweeps_base,r'resonator_{}{}_updowndown'.format(distance,letter))
            files_ampsweeps = glob(path.join(ampsweeps, r'*.npy'))
            
            # files = files[:-53]
            files.sort(key=lambda f: np.load(f, allow_pickle=True).item()['App_gen (V)'])
            
            
            folder = switching_sweeps_names[temp]
            basedir = r'..\2023-03-Helmholtz_resonators\Helmholtz_res_drive_dep\{}\buffer'.format(folder)
            
            # params = np.load(f'lifetimes_params\params_{letter}_T{temp}.npy', allow_pickle=True).item()
            # y1_list = params['y1']
            # y2_list = params['y2']
            # dy_list = params['dy']
            
            #%% analyse switching
            amplitudes = []
            lifetimes_up=[]
            lifetimes_down=[]
            
            higher_stack = []
            lower_stack = []
            mean_rs = []
            y1s = []
            y2s =[]
            
            amp_last = -1
            y2_last = 0
            y1_last = 0
            r_last = 0
            k=0
            for i,file in enumerate(files[:]): #and files[39:]:
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
                    datapoints = ax.plot(t,r,'.-',ms=1)
                    fig.canvas.draw()
                    
                    if i==51:
                        t11=20.846
                        t12=21.06
                        
                        t21=21.123
                        t22=21.198
                        select1=np.logical_and(t>t11,t<t12)
                        select2=np.logical_and(t>t21,t<t22) 
                        
                            #T,A,x0,A0
                        p01 = [0.02,5e-10,t11,1.52e-8]
                        p02 = [0.03,5e-10,t21,1.52e-8]
                        # ax.plot(t[select1],up(t[select1],*p01),'r--')
    
                        par1,sig1 = curve_fit(up,t[select1],r[select1],p0 = p01)
                        par2,sig2 = curve_fit(down,t[select2],r[select2],p0 = p02)
                        
                        ax_fit.plot(t,r,'.-',ms=1)
                        ax_fit.plot(t[select1],up(t[select1],*par1),'r-')
                        ax_fit.plot(t[select2],down(t[select2],*par2),'r-')
                        ax_fit.set_xlim(np.min([t11,t22,t12,t21]),np.max([t11,t22,t12,t21]))
                        fig_fit.canvas.draw()

                        
                        
                        print(f'T_up: {par1[0]:.2e}s\nT_down: {par2[0]:.2e}s')
                        
                    # while not plt.waitforbuttonpress():
                    #     pass
                    # print(i)
                    datapoints[0].remove()
                    

                    
                
                    

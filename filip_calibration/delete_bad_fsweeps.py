# -*- coding: utf-8 -*-
"""
Created on Wed May 17 18:41:25 2023

@author: Admin
"""
    
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
from glob import glob
import scipy as sc
import sys 
from scipy.optimize import curve_fit
import helium_old as he

from fitting_class import MicrowavePeak
from fitting_class import microwave_peak

from Bridge_losses_class import bridgelosses
from Circuit_losses_class import circuitlosses

temps = ['T1325', 'T1350_4', 'T1375', 'T1400', 'T1425', 'T1450_2','T1475', 'T1500', 'T1550', 'T1650_2', 'T1700', 'T1750', 'T1850', 'T1950']


plt.rcdefaults()
resonators_geometry = {
    'B':{        
        'w':  500e-6,
        'l':  750e-6,
        'C0': 335.6978986554392e-12,
        'kp': 1.4546169638969617e7
        
        },
    'D':{
        'w': 500e-6, 
        'l': 950e-6,
        'C0':466.4158697787019e-12,
        'kp':1.6174502699593267e7
        }}

plt.close('all')
save =True
for temp in temps[10:]:
    
    T = float(temp[1:].split('_')[0])*1e-3
    print(f'\ntemp: {T}\n')   

    
    letter = 'B'
    height = 500
    Temp = temp
    
    Nbot = 2
    Ntop = 9
    

    fitx = True
    fity = True
    
    folder = temp
    
    resonator = f'{height:}{letter:}'
    
    
    fsweeps =       rf'D:\Github\Masters-thesis-analysis\filip_calibration\Helmholtz_res_drive_dep\{temp}\fsweeps\resonator_{height:}{letter:}\*.npy'
    freqsweep_all = rf'D:\Github\Masters-thesis-analysis\filip_calibration\Helmholtz_res_drive_dep\{temp}\fsweeps\resonator_{height:}{letter:}\*.npy'
                    
    
    Ffiles = glob(fsweeps)
    Freqfiles = glob(freqsweep_all)
    
    '''background estimation'''
    par = []
    a1s = []
    b1s = []
    a2s = []
    b2s = []
    drs = []
    correct_drives = []
    correct_ys = []
    
    fig,ax = plt.subplots(2,1)
    
    
    # Ffiels = Ffiles.sort(key=lambda x : np.load(x,allow_pickle=True).item()['App (V)'])
    Freqfiles = np.sort(Freqfiles)
    
    for num, file in enumerate(Ffiles[102:105]):
        #print(num)
        d = np.load(file, allow_pickle=True).item()
        
        f = np.array(d['freq (Hz)'])
        x = np.array(d['x (A)'])
        y = np.array(d['y (A)'])
        r = np.sqrt(x**2 + y**2)
        amp = d['App (V)']
        
        if num == 0:
            A= np.max(np.sqrt(x**2 + y**2))
            ix1 = np.argmax(x)
            ix2 = np.argmin(x)
            gamma_est = abs(f[ix1] - f[ix2])/2
            ixf = np.argmax(np.sqrt(x**2 + y**2))
            f0_est = f[ixf]
            A= np.max(np.sqrt(x**2 + y**2))*gamma_est
            init_p = [f0_est, A, gamma_est,0,np.mean(y),0,np.mean(x),0]
        else:
            init_p = [par['f0'],par['A'],par['gamma'],par['a1'],par['b1'],par['a2'],par['b2'],par['phi']]
        
        """
        fit frequency sweep with linear background
        """      
        # ax[0].plot(f, np.real(microwave_peak(f,*init_p)), 'r--')
        # ax[1].plot(f, np.imag(microwave_peak(f,*init_p)), 'r--')
        peak = MicrowavePeak(f, x, y, init_p)
        result = peak.fit(plot=True)
        par = result.best_values
        res_fit = result.best_fit
        ax[0].plot(f, x, 'o',ms = 3)
        ax[1].plot(f, y, 'o', ms = 3)
        index1 = f < (f[0] + 50)
        index2 = f>f[-1]-10

        array1 = y[index1]
        std1 = np.std(array1)
        
        array2 = y[index2]
        std2 = np.std(array2)
        print(std1,'\n',std2)
        if std1 > 1e-10 or std2 > 1e-10:
            plt.waitforbuttonpress()
            print('remove?')
            decision = input()
            
            if decision == 'y':
                os.remove(file)
                
        # ax[0].plot(f, np.real(res_fit), 'k--')
        # ax[1].plot(f, np.imag(res_fit), 'k--')
        """
        get linear background params
        a1 -> lin y
        b1 -> const y
        a2 -> lin x
        b2 -> const x
        """
        # a1s.append(par['a1'])
        # b1s.append(par['b1'])
        # a2s.append(par['a2'])
        # b2s.append(par['b2'])
        # drs.append(amp) #drives
        
            
    # fig2,ax2 = plt.subplots(2,2)
    # ax2[0,0].plot(drs,a1s,'o')
    # ax2[1,0].plot(drs,b1s,'o')
    # ax2[0,1].plot(drs,a2s,'o')
    # ax2[1,1].plot(drs,b2s,'o')
    while not plt.waitforbuttonpress():
        pass
    
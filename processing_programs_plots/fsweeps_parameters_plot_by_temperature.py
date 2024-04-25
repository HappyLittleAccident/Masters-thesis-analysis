# -*- coding: utf-8 -*-
"""
Created on Fri Jan 26 16:46:56 2024

@author: Marek
"""

import numpy as np
import matplotlib.pyplot as plt
import helium_old as he
from glob import glob
import os
from scipy.interpolate import CubicSpline
import sympy as sp

plt.close('all')

figres,axres = plt.subplots(1,1,sharex=True)    

figdamp,axdamp = plt.subplots(1,1,sharex=True)    

figgraph,axgraph = plt.subplots(2,1,sharex=True)

folder = r'D:\OneDrive_copy\OneDrive - Univerzita Karlova\DATA\2023_4_Helmholtz_resonators_recal_2\resonance_freq'

files = glob(os.path.join(folder,'*.npy'))



def compres (T):
    a,b,c,d = ( 3.65456274e-10, -1.19112443e-09,  1.72705353e-09, 1.12838569e-08)
    chi = a*T**3 + b*T**2 + c*T + d
    return chi*10


T_spline = [0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0]
alpha_spline = np.array([4.5, 3.3, 0.8, -3.2, -9.5, -18.0, -26.8, -38.0, -52.2, -71.3, -97.6, -129.6])*1e-4
expand = CubicSpline(T_spline,alpha_spline)


resonators_geometry = {
    'A':{        
        'w':  1000e-6,
        'l':  1.601*1050e-6,
        'C0': 338.63213961685864e-12,
        'kp': 3*1.2400632746323976e7
        },

    
    'B':{        
        'w':  1000e-6,
        'l':  1.55*750e-6,
        'C0': 335.6978986554392e-12,
        'kp': 2*1.4546169638969617e7
        },
    
    'C':{        
        'w':  950e-6,
        'l':  1.03*1000e-6,
        'C0': 312.8732100355237e-12,
        'kp': 2*1.1909096177342793e7
        },
    'D':{
        'w': 1000e-6, #orig= 500e-6
        'l': 1.405*950e-6,
        'C0':466.4158697787019e-12,
        'kp':2*1.6174502699593267e7
        }}



A = np.pi*(2.5e-3)**2
colors = ['r','b','tab:orange','k']
for i,file in enumerate(files[:]):
    letter = file.split('\\')[-1].split('.')[0][-1]

    d = np.load(file,allow_pickle=True).item()
    temps = np.array(d['temperature (K)'])
    amps = d['drive (Vrms)']
    f0s = d['frequency (Hz)']
    ws = d['width (Hz)']


    resonator = resonators_geometry[letter]
    k = resonator['kp']
    l = resonator['l']
    for j,T in enumerate(temps[:]):
        axres.clear()
        axres.plot(amps[j],ws[j],'o',c=colors[i],label='Measured')
        figres.canvas.draw()
        while not plt.waitforbuttonpress():
            pass
    

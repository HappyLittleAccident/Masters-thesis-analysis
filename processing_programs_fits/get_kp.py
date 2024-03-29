# -*- coding: utf-8 -*-
"""
Created on Fri Feb 23 13:02:15 2024

@author: Marek
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from glob import glob


folder  = r'D:\OneDrive_copy\OneDrive - Univerzita Karlova\DATA\2023_4_Helmholtz_resonators_recal_2\chip_calibration\resonator_500A'

def lin(t,a,b):
    return a*t + b

resonators_geometry = {
    'A':{        
        'w':  1000e-6,
        'l':  1050e-6,
        'C0': 338.63213961685864e-12,
        'kp': 1.2254e+07 #1.2400632746323976e7
        },

    
    'B':{        
        'w':  1000e-6,
        'l':  55750e-6,
        'C0': 335.6978986554392e-12,
        'kp': 1.4546169638969617e7
        },
    
    'C':{        
        'w':  950e-6,
        'l':  1000e-6,
        'C0': 312.8732100355237e-12,
        'kp': 1.1909096177342793e7
        },
    'D':{
        'w': 1000e-6, #orig= 500e-6
        'l': 950e-6,
        'C0':466.4158697787019e-12,
        'kp':1.6174502699593267e7
        }}
D = 500e-9
epsilon = 8.854e-12

plt.close('all')
files = glob(folder+'\\*.npy')

biases = np.zeros(len(files[6:]))
C_corrs = np.zeros(len(files[6:]))
fig,ax=plt.subplots()
for i,file in enumerate(files[6:]):
    d = np.load(file,allow_pickle=True).item()
    bias = d['bias (V)']
    f = np.array(d['freq (Hz)'][1:])
    C = d['capacitance (F)'][1:]
    C_s = d['C serial (F)']
    I = np.abs(np.array(d['x (A)'][1:])+1j*np.array(d['y (A)'][1:]))
    U = d['amplitude (Vrms)']
    C_ch = resonators_geometry['A']['w']*resonators_geometry['A']['l']*epsilon/D
    print('file no.',i)

    par_C,sig_C = curve_fit(lin,2*np.pi*f,I)
    C_fit = par_C[0]/U
    C_corr = 1/(1/C_fit - 1/C_s) - 2*C_ch
    
    biases[i]=bias
    C_corrs[i]=C_corr
    
    ax.plot(bias**2,C_fit,'bo')#,label = f)    
    ax.plot(bias**2,C_corr,'ro')
    
par_k,sig_k = curve_fit(lin,biases**2,C_corrs)

kp=par_k[1]**2/par_k[0]/D**2
print(f'{kp:.4e}')
# ax.legend()
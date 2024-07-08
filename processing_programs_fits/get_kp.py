# -*- coding: utf-8 -*-
"""
Created on Fri Feb 23 13:02:15 2024

@author: Marek
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from glob import glob
from format_figure import use_latex,polish

plt.close('all')
use_latex()

def current(omega,CU,UdivR,fi):
    curr = (1j*omega*CU + UdivR)*np.exp(1j*fi)
    return np.concatenate([
                           np.real(curr),
                           np.imag(curr)
                           ])

def lin(Ubsq,sigma,C0):
    return Ubsq*sigma + C0

resonators_geometry = {
    'A':{        
        'w':  1000e-6,
        'l':  1050e-6,
        'C0': 338.63213961685864e-12,
        'kp': 1.2254e+07 #1.2400632746323976e7
        },

    
    'B':{        
        'w':  1480e-6,
        'l':  650e-6,
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
fig,ax=plt.subplots(dpi=200)
fig_c,ax_c = plt.subplots()
resonators = ['A','B','C']
geometry = ['C','J','R']
colors = ['tab:blue','tab:orange','tab:green']
for j,resonator in enumerate(resonators[:]):
    
    folder  = rf'D:\OneDrive_copy\OneDrive - Univerzita Karlova\DATA\2023_4_Helmholtz_resonators_recal_2\chip_calibration\resonator_500{resonator}'
    
    files = glob(folder+'\\*.npy')
    
    biases = []
    C_corrs = []

    if resonator == 'C':
        start = 104
        stop = 200
    elif resonator == 'B':
        start = 4
        stop = -1
    else:
        start = 6
        stop = -1
    for i,file in enumerate(files[start:stop]):
        d = np.load(file,allow_pickle=True).item()
        bias = d['bias (V)']
        f = np.array(d['freq (Hz)'][10:])
        C = d['capacitance (F)'][10:]
        C_s = d['C serial (F)']
        I = np.array(d['x (A)'][10:])+ 1j*np.array(d['y (A)'][10:])
        U = d['amplitude (Vrms)']
        C_ch = resonators_geometry[resonator]['w']*resonators_geometry[resonator]['l']*epsilon/D

    
        par_C,sig_C = curve_fit(current,2*np.pi*f,np.concatenate([I.real,I.imag]))
        C_fit = par_C[0]/U
        C_corr = 1/(1/C_fit - 1/C_s) - 2*C_ch


        ax_c.clear()
        ax_c.plot(f,I.real,'ko')
        ax_c.plot(f,I.imag,'bo')
        ax_c.plot(f,current(2*np.pi*f,*par_C)[:len(f)])
        ax_c.plot(f,current(2*np.pi*f,*par_C)[len(f):])
        fig_c.canvas.draw()
        errs = np.abs(np.diag(sig_C))
        # print(f'sig fit: {errs}')

        if (errs[0]>3e-30 or errs[0]<1e-31) and resonator == 'C':
                
            # while (not plt.waitforbuttonpress()):
                pass
        else:
        
            biases.append(bias)
            C_corrs.append(C_corr)
        
            # ax.plot(bias**2,C_fit,'bo')#,label = f)    
            ax.plot(bias**2,C_corr*1e12,'o',c=colors[j])
            
            fig.canvas.draw()
    biases = np.array(biases)
    C_corrs = np.array(C_corrs)
    par_k,sig_k = curve_fit(lin,biases**2,C_corrs)
    
    U_bs = np.linspace(0,100,30)
    ax.plot(U_bs,lin(U_bs,*par_k)*1e12,'k--')
    ax.plot([],[],'o',c=colors[j],label = 'Channel geometry '+geometry[j])

    
    kp=3*par_k[1]**2/par_k[0]/D**2
    print(f'{kp:.4e} filip: {3*resonators_geometry[resonator]["kp"]:.4e}')
    #%%
ax.plot([],[],'k--',label=r'fit to $C_M = C_0 + \sigma U_B^2$')
ax.legend()
ax.set_xlabel(r'$U_B^2 \mathrm{(V^2)}$')
ax.set_ylabel(r'$C_M$ (pF)')
#%%
polish(fig, 1, name=r'D:\Github\Masters-thesis-analysis\processing_programs_plots\images\calibrate_kp', extension='.png', grid=True,width_to_height = 0.7,tight_layout=True)        



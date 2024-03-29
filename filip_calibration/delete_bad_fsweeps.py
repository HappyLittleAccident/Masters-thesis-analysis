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
import scipy as sp

from fitting_class import MicrowavePeak
from fitting_class import microwave_peak

from Bridge_losses_class import bridgelosses
from Circuit_losses_class import circuitlosses

temps = ['T1325', 'T1350_4', 'T1375', 'T1400', 'T1425', 'T1450_2','T1475',
         'T1500', 'T1550', 'T1650_2', 'T1700', 'T1750', 'T1850', 'T1950']
resonators = ['D']

"""
last temp finished 1.3750K
"""

plt.rcdefaults()
# resonators_geometry = {
#     'B':{        
#         'w':  500e-6,
#         'l':  750e-6,
#         'C0': 335.6978986554392e-12,
#         'kp': 1.4546169638969617e7
        
#         },
#     'D':{
#         'w': 500e-6, 
#         'l': 950e-6,
#         'C0':466.4158697787019e-12,
#         'kp':1.6174502699593267e7
#         }}

plt.close('all')
save =True
fig_remove,ax_remove = plt.subplots(2,1)
fig,ax = plt.subplots(2,1)
cmap = plt.get_cmap('viridis')

rightclick = False

def onclick(event):
    global rightclick
    if event.button == 3:
        rightclick = True
    else:
        rightclick = False

for resonator in resonators[:]:
    for temp in temps[:]:
        folder = temp
        T = float(temp[1:].split('_')[0])*1e-3
        print(f'\ntemp: {T}\n')   
    
        if T <0:
            continue
        
        letter = resonator
        height = 500
        
        
        
        fsweeps = rf'D:\Github\Masters-thesis-analysis\filip_calibration\Helmholtz_res_drive_dep\{temp}\fsweeps_ampshift\resonator_{height:}{letter:}\*.npy'
        """          D:\Github\Masters-thesis-analysis\filip_calibration\Helmholtz_res_drive_dep\T1350_4\fsweeps\resonator_500A
        """
        Ffiles = glob(fsweeps)

        if Ffiles == []:
            continue

                    
        batch_start = 0
        amp_last = 0
        for num, file in enumerate(Ffiles[:]):
            
            d = np.load(file, allow_pickle=True).item()
            
            f = np.array(d['freq (Hz)'])
            x = np.array(d['x (A)'])
            y = np.array(d['y (A)'])
            r = np.sqrt(x**2 + y**2)
            amp = d['App (V)']

            ax[0].plot(f, x, 'o-',alpha=0.5,ms = 1,c=cmap(amp/7))
            ax[1].plot(f, y, 'o-',alpha=0.5, ms = 1,c=cmap(amp/7))
            
            ax[0].plot(f[0],x[0],'ro')
            ax[0].plot(f[-1],x[-1],'go')
            ax[1].plot(f[0],y[0],'ro')
            ax[1].plot(f[-1],y[-1],'go')
            
        print('all files')
        fig.canvas.draw()
        fig.canvas.flush_events()
        
        while not plt.waitforbuttonpress():
            pass
        
        ax[0].clear()
        ax[1].clear()
        
        for num, file in enumerate(Ffiles[:]):
            
            d = np.load(file, allow_pickle=True).item()
            
            f = np.array(d['freq (Hz)'])
            x = np.array(d['x (A)'])
            y = np.array(d['y (A)'])
            r = np.sqrt(x**2 + y**2)
            amp = d['App (V)']
            
            fin = num == len(Ffiles)-1
            if amp_last > amp or fin:
                if fin:
                    ax[0].plot(f, x, 'o-',alpha=0.5,ms = 1,c=cmap(amp/7))
                    ax[1].plot(f, y, 'o-',alpha=0.5, ms = 1,c=cmap(amp/7))
                    ax[0].plot(f[0],x[0],'ro')
                    ax[0].plot(f[-1],x[-1],'go')
                    ax[1].plot(f[0],y[0],'ro')
                    ax[1].plot(f[-1],y[-1],'go')
                print('one batch')
                fig.canvas.draw()
                fig.canvas.flush_events()
                
                while not plt.waitforbuttonpress():
                    pass
                
                print('rightlick for manual selection')

                plt.waitforbuttonpress()
                fig.canvas.mpl_connect('button_press_event', onclick)
                plt.waitforbuttonpress()
                
                if rightclick:
                    print('press to continue')
                    plt.waitforbuttonpress()
                    print('pressed')
                    end = num
                    if fin:
                        end = num+1
                    for file_remove in Ffiles[batch_start:end]:
                        ax_remove[0].clear()
                        ax_remove[1].clear()
                        
                                    
                        dr = np.load(file_remove, allow_pickle=True).item()
                        
                        fr = np.array(dr['freq (Hz)'])
                        xr = np.array(dr['x (A)'])
                        yr = np.array(dr['y (A)'])
                        ax_remove[0].plot(fr, xr, 'o-',ms = 1)
                        ax_remove[1].plot(fr, yr, 'o-', ms = 1)
                        

                        
                        fig_remove.canvas.draw()
                        fig_remove.canvas.flush_events()
                        
                        button = plt.waitforbuttonpress()
                        if button:
                            os.remove(file_remove)
                            print('removed')
                    ax[0].clear()
                    ax[1].clear()
                    print('removing finished')
                else:
                    ax[0].clear()
                    ax[1].clear()
                batch_start = num
            amp_last = amp
            ax[0].plot(f, x, 'o-',alpha=0.5,ms = 1,c=cmap(amp/7))
            ax[1].plot(f, y, 'o-',alpha=0.5, ms = 1,c=cmap(amp/7))
            ax[0].plot(f[0],x[0],'ro')
            ax[0].plot(f[-1],x[-1],'go')
            ax[1].plot(f[0],y[0],'ro')
            ax[1].plot(f[-1],y[-1],'go')
        ax[0].clear()
        ax[1].clear()
        fig.canvas.draw()

        
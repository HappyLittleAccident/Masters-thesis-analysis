# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 16:55:10 2023

@author: Marek
"""

import numpy as np
import matplotlib.pyplot as plt
import os
from glob import glob
from format_figure import polish, use_latex
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter
import matplotlib as mpl

plt.close('all')
folder = r'..\filip_calibration\2023_4_Helmholtz_resonators_recal_2\Helmholtz_fsweeps'

temps = ['T1325', 'T1350', 'T1375', 'T1400', 'T1425', 'T1450', 'T1475',
         'T1500', 'T1550', 'T1650', 'T1700', 'T1750', 'T1850', 'T1950']
resonators = ['500A', '500B', '500C']
names = ['C', 'J', 'R']

cmap = plt.get_cmap('viridis')


use_latex()
# plt.rcdefaults()
def lin(x,a,b):
    return a*x + b

colors = ['tab:blue', 'k', 'tab:orange']
for i, resonator in enumerate(resonators[2:3]):
    figamp, axamp = plt.subplots()
    figwidth, axwidth = plt.subplots(2,1,dpi=250)
    cmaps = ['winter', 'winter', 'winter']
    inferno = plt.get_cmap(cmaps[i])
    if resonator == '500C':
        mult = 1
    else:
        mult = 1
    log = 0
    for log in [0,1]:
        for temp in temps[0:-1:1+11*log]:
            print(temp)
            filenames = glob(os.path.join(folder, temp, f'{resonator}', '*.npy'))
            A_max = np.max([np.load(file_temp, allow_pickle=True).item()[
                           'App_gen (V)'] for file_temp in filenames])
            widths = np.zeros(len(filenames))
            amps = np.zeros(len(filenames))
            drives = np.zeros(len(filenames))
            index = False
            for j, file in enumerate(filenames[0::]):
    
                data = np.load(file, allow_pickle=True).item()
                f = data['Freq (Hz)']
                r = data['Velocity (m/s)']*100
                A = data['App_gen (V)']*mult
                
                if f[r>np.max(r)/2][-1] - f[r>np.max(r)/2][0] > 90:
                    continue
                
                if j > 5:
                    if A < drives[j-1]:
                        break
    
                drives[j] = A
                amps[j] = np.max(r)
    
                widths[j] = f[r > amps[j]/2][-1]-f[r > amps[j]/2][0]
    
    
    
                if (not index) and widths[j] > 1.08*np.mean(widths[0:j]) and j > 5:
                    index = j
                
    
            widths = widths[widths != 0]
            amps = amps[amps != 0]
            drives = drives[drives != 0]
    
            
            amps = amps/12.041535654092046 if resonator == '500C' else amps
    
            
            def poly(x, b, d, x0):
                return ((-b/2/x0)*x**3 + b*x**2 + (-x0*b/2)*x)*np.heaviside(x-x0, 1) + d
            if resonator == '500C':
                try:
                    par, sig = curve_fit(lin, amps[amps<0.45][3:], widths[amps < 0.45][3:], p0=[
                                          0, np.mean(widths[amps<0.4])])
                except:
                    print(np.shape(widths), np.shape(amps))
    
    
            # widths_hat = savgol_filter(widths, len(widths)//3, 3, mode='mirror')
            # axwidth.plot(amps,poly(amps,*par),'k--')
            # axwidth.axhline(np.mean(widths[:index]))
    
            scale = (float(temp[1:])*1e-3 - 1.325)/(1.95 - 1.325)
    
            if log:
                axwidth[1].semilogy(amps,
                    widths-(
                        np.mean(widths[:index]) if resonator != '500C' else np.mean(widths[amps<0.4])#lin(amps,*par)
                            )+1, '.', lw=0.5, c=inferno(scale))
                # axwidth.plot(amps,np.log(widths_hat-np.mean(widths[:index])),c=inferno(scale))
    
            else:
                axwidth[0].plot(amps, widths, '.', c=inferno(scale))
                # axwidth.plot(amps,lin(amps,*par),'r-')
                # axwidth.plot(amps,widths_hat,c=inferno(scale))
            # axwidth.plot(amps[index],widths[index],'ko')
            # axwidth.plot(par[2],par[1],'ko')
    
            axamp.plot(drives, amps, '-', c=inferno(scale), label=temp)
            # axwidth.text(amps[-1]+0.01,widths[-1]+2,temp[1:2]+'.'+temp[2:] + ' K',rotation = 45)
            #
    
            # axamp.set_title(resonator + 'amp')
            # axwidth.set_title(resonator + 'width')
        if resonator == '500C':
            axwidth[1].set_xlabel('Superfluid velocity (Arb.)')
        else:
            axwidth[1].set_xlabel('Superfluid velocity (cm/s)')
        if log:
            axwidth[1].set_ylabel('FWHM (Hz)')
            axwidth[1].set_ylim(0.4,None)
        else:
            axwidth[0].set_ylabel('FWHM (Hz)')

        axamp.legend()
    figwidth.tight_layout()
    figwidth.colorbar(plt.cm.ScalarMappable(norm=mpl.colors.Normalize(
        vmin=1.325, vmax=1.85, clip=False), cmap=inferno), ax=axwidth, label='Temperature (K)')

    polish(figwidth, 1, name=f'images//widths{names[i]}'+(
        'log' if False else ''), extension='.png', grid=True)

#%%
for i, resonator in enumerate(resonators[0:3]):
    cmaps = ['winter', 'winter', 'winter']
    inferno = plt.get_cmap(cmaps[i])
    fig, ax = plt.subplots(1, 1, sharex=True,dpi=250)
    if resonator == '500C':
        mult = 1
    else:
        mult = 1
    for temp in temps[0:1]:
        print(temp)
        filenames = glob(os.path.join(folder, temp, f'{resonator}', '*.npy'))
        A_max = np.max([np.load(file_temp, allow_pickle=True).item()[
                       'App_gen (V)'] for file_temp in filenames])
        widths = np.zeros(len(filenames))
        amps = np.zeros(len(filenames))
        drives = np.zeros(len(filenames))
        index = False
        for j, file in enumerate(filenames[0::2]):

            data = np.load(file, allow_pickle=True).item()
            f = data['Freq (Hz)']
            r = data['Velocity (m/s)']*100
            A = data['App_gen (V)']*mult
            ax.plotw(f,r,'.',c=cmap(A/A_max))
    
    ax.set_xlabel('Driving frequency (Hz)')
    ax.set_ylabel('Superfluid velocity (cm/s)')
    
    fig.tight_layout()
    fig.colorbar(plt.cm.ScalarMappable(norm=mpl.colors.Normalize(
        vmin=1.325, vmax=1.85, clip=False), cmap=cmap), ax=ax, label='Pressure gradient (Pa/m)')
    
    polish(fig, 1, name=f'images//fsweeps{names[i]}', extension='.png', grid=True,width_to_height = 1.5)

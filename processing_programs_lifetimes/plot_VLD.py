# -*- coding: utf-8 -*-
"""
Created on Thu Jul 27 11:22:38 2023

@author: Marek
"""
import numpy as np
import matplotlib.pyplot as plt
import os
from glob import glob
from scipy.optimize import curve_fit

plt.close('all')

folder = r'D:\OneDrive\OneDrive - Univerzita Karlova\DATA\2023_4_Helmholtz_resonators_recal\Helmholtz_res_VLD'

v=2
def function(x,a,b,c,d,e):
    # return a*x**2 + b*x + c
    # return a+ d*x - b*np.heaviside(x-c,0)*(x-e)
    # return np.sqrt(a + d*(x + (b*x/(c**2 + (x)**2))**2))
    return np.sqrt(d*x + (b*x/(c**2 + (x/e)**2))**2)+a


temperature_folders = glob(folder+'\\*')
temps_strings = [name[-4:] for name in temperature_folders]

for temp in temps_strings[:]:
    filenames = glob(folder + f'\\T{temp}\\500D\\*.npy')
    
    # fig,ax = plt.subplots()
    fig_hist,ax_hist = plt.subplots()
    
    upVs   = []
    downVs = []
    jumpVs = []
    
    upLs   = []
    downLs = []
    jumpLs = []
    first = True
    for filename in filenames[:]:
        
        data = np.load(filename,allow_pickle= True).item()

        upV = data['upV (m/s)']
        downV = data['downV (m/s)']
        jumpV = data['jumpV (m/s)']
        
        upVs.append(upV)
        downVs.append(downV)
        jumpVs.append(jumpV)
        
        upL = data['upL (m^-2)']
        downL = data['downL (m^-2)']
        jumpL = data['jumpL (m^-2)']
        
        upLs.append(upL)        
        downLs.append(downL)
        jumpLs.append(jumpL)
        
        # if first:        
        #     ax.plot(upV,upL,'.',c='tab:blue',ms=1)
        #     # ax.plot(downV,downL,c='tab:orange')
        #     # ax.plot(jumpV,jumpL,c='tab:green')
        #     ax.set_xlabel('V ($m s^{-1}$)')
        #     ax.set_ylabel('VLD ($m^{-2}$)')
        #     ax.set_title(f'Temp: {temp}K')
        #     fig.savefig(fr'D:\_personal\_lab\talk_3_8_2023\ampsweeps_VLD_{temp}K.png')
        #     first = False
    
    upVs = np.concatenate(upVs)
    downVs = np.concatenate(downVs)
    jumpVs = np.concatenate(jumpVs)
    
    upLs = np.concatenate(upLs)
    downLs = np.concatenate(downLs)
    jumpLs = np.concatenate(jumpLs)
    
    Vs = np.concatenate([upVs,downVs,jumpVs])
    Ls = np.concatenate([upLs,downLs,jumpLs])
    
    # Ls/=np.max(Ls)
    # Vs/=np.max(Vs)
    
    hist = ax_hist.hist2d(Vs,Ls,bins = 600,cmap='inferno',norm='log',alpha=1)
    fig_hist.colorbar(hist[3], ax=ax_hist)
    
    ax_hist.set_xlabel('V ($m s^{-1}$)')
    ax_hist.set_ylabel('VLD ($m^{-2}$)')
    ax_hist.set_title(f'Temp: {temp}K')
    # ax_hist.plot(Vs,Ls)
    fig_hist.tight_layout()
    fig_hist.savefig(fr'D:\_personal\_lab\talk_3_8_2023\ampsweeps_VLD_hist{temp}K.png')
    p0s = [0.38,  1.1 , 1 ,0.38,5 ]
    t = np.linspace(0,np.max(Ls))
    # plt.figure()
    # plt.plot(function(t,*p0s),t,'r--')
    # ax_hist.plot(function(t,*p0s),t,'b--')
    # try:
    #     par,sig = curve_fit(function,Ls,Vs,p0 = p0s)
    #     ax_hist.plot(function(t,*par),t,'b-')
    #     print(par)
    # except:
    #     print('fit fail')
    

        
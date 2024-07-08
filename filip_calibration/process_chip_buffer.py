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

from Bridge_losses_class import bridgelosses
from Circuit_losses_class import circuitlosses

temps = ['1.35','1.45','1.65','1.85']
switching_sweeps_names = {'1.35':'T1350_3','1.45':'T1450_3','1.65':'T1650_3','1.85':'T1850_2'}
Temps = {'1.35':'T1350','1.45':'T1450','1.65':'T1650_2','1.85':'T1850'}

plt.rcdefaults()
resonators_geometry = {
    'B':{        
        'w':  500e-6,
        'l':  750e-6,
        'C0': 335.6978986554392e-12,
        'kp': 3*1.4546169638969617e7
        
        },
    'D':{
        'w': 500e-6, 
        'l': 950e-6,
        'C0':466.4158697787019e-12,
        'kp':3*1.6174502699593267e7
        }}

plt.close('all')
save =True
for temp in temps[:]:
    
    T = float(temp)
    print(T)
    def linfit(x,a,b):
        y = a*x + b
        return y
    
    def cubicfit(x,a,b,c,d):
        y = a*x**3 + b*x**2 + c*x + d
        return y
    
    def compres (T):
        a,b,c,d = ( 3.65456274e-10, -1.19112443e-09,  1.72705353e-09, 1.12838569e-08)
        chi = a*T**3 + b*T**2 + c*T + d
        return chi*10
      
    

    
    letter = 'B'
    height = 500
    Temp = Temps[temp]
    
    Nbot = None
    Ntop = 9
    

    fitx = True
    fity = True
    
    folder = switching_sweeps_names[temp]
    
    resonator = f'{height:}{letter:}'
    
    
    fsweeps = rf'..\..\2023-03-Helmholtz_resonators\Helmholtz_res_drive_dep\{Temp}\fsweeps\resonator_{height:}{letter:}\*.npy'
    buffer = rf'..\..\2023-03-Helmholtz_resonators\Helmholtz_res_drive_dep\{folder}\buffer\resonator_{height:}{letter:}\*.npy'
    
    
    Ffiles = glob(fsweeps)
    BufferFiles = glob(buffer)
    
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
    
    
    Ffiles = np.sort(Ffiles)
    BufferFiles = np.sort(BufferFiles)
    for num, file in enumerate(Ffiles[10:20]):
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
            init_p = [f0_est, A, gamma_est,0,0,0,0,0]
        else:
            init_p = [par['f0'],par['A'],par['gamma'],par['a1'],par['b1'],par['a2'],par['b2'],par['phi']]
        """
        fit frequency sweep with linear background
        """           
        peak = MicrowavePeak(f, x, y, init_p)
        result = peak.fit(plot=True)
        par = result.best_values
        res_fit = result.best_fit
        ax[0].plot(f, x, 'o',ms = 3)
        ax[1].plot(f, y, 'o', ms = 3)
        ax[0].plot(f, np.real(res_fit), 'k--')
        ax[1].plot(f, np.imag(res_fit), 'k--')
        """
        get linear background params
        a1 -> lin y
        b1 -> const y
        a2 -> lin x
        b2 -> const x
        """
        a1s.append(par['a1'])
        b1s.append(par['b1'])
        a2s.append(par['a2'])
        b2s.append(par['b2'])
        drs.append(amp) #drives
            
    fig2,ax2 = plt.subplots(2,2)
    ax2[0,0].plot(drs,a1s,'o')
    ax2[1,0].plot(drs,b1s,'o')
    ax2[0,1].plot(drs,a2s,'o')
    ax2[1,1].plot(drs,b2s,'o')
    
    def background (f,a,b):
        y = a*f + b
        return y
    
    """
    Perform linear fit of background parameters as function of drive amplitude,
    or chose background params as mean
    """
    if fitx:
       drs = np.array(drs)
       a2s = np.array(a2s)
       b2s = np.array(b2s)
       para,resa = curve_fit(linfit, drs, a2s) 
       parb,resb = curve_fit(linfit, drs, b2s) 
       ax2[0,1].plot(drs, linfit(drs,para[0],para[1]),'k--')
       ax2[1,1].plot(drs, linfit(drs,parb[0],parb[1]),'k--')  
    else:
        a2mean = np.mean(a2s[:])
        b2mean = np.mean(b2s[:])
        
    if fity:
        drs = np.array(drs)
        a1s = np.array(a1s)
        b1s = np.array(b1s)
        para2,resa2 = curve_fit(linfit, drs, a1s) 
        parb2,resb2 = curve_fit(linfit, drs, b1s) 
        ax2[0,0].plot(drs, linfit(drs,para2[0],para2[1]),'k--')
        ax2[1,0].plot(drs, linfit(drs,parb2[0],parb2[1]),'k--')  
    else:
        a1mean = np.mean(a1s[:3])
        b1mean = np.mean(b1s[:3])    
        print(a1mean,b1mean)
        ax[1].plot(f,background(f,a1mean,b1mean),'r--')    
    
    
    '''background + losses subtraction + velocity and force calculation'''
    #sys.path.append(r'D:\python\bridge_losses')
    
    def Pgrad (bias,drive,T):
        R = 2.5e-3
        h = 500e-9
        A = np.pi*R**2
        VB = A*h
        l = resonators_geometry[letter]['l']
        kp = resonators_geometry[letter]['kp']
        C0 = resonators_geometry[letter]['C0']
        chi = compres(T)
        Pgrad = (2*A*C0*bias)/(chi*VB*kp + 2*A**2)*drive/(h*l)
        return Pgrad
    
    def Velocity (bias, I,T):
        R = 2.5e-3
        h = 500e-9
        A = np.pi*R**2
        VB = A*h
        rho = 145.1281
        er = 1.0565
        e0 = 8.854e-12
        w = resonators_geometry[letter]['w']
        a = h*w
        chi = compres(T)
        rhos = he.superfluid_density(T)
        #print(rhos)
        #print(chi)
        kp = resonators_geometry[letter]['kp']
        C0 = resonators_geometry[letter]['C0']
        
        SIG = chi*h*kp/(2*A)
        #print(SIG)
        g = (2*a*rhos)/(VB*rho)*(1 - 2*(er - 1)/er*SIG)/(1 + 2*SIG)
        v = I/(g*bias*C0)
        return v
        
    # fig3,ax3 = plt.subplots(1,1)
    # fig4,ax4 = plt.subplots(1,1)
    for file in BufferFiles: 
        d = np.load(file,allow_pickle=True).item()
        amplitude = np.array(d['App_gen (V)'])
        
        x = np.array(d['x (A)'])
        y = np.array(d['y (A)'])
        f0 = d['res freq (Hz)']
        sample_rate = d['sample rate (Hz)']
        t = d['time']
        bias = d['bias (V)']
        
        
        Amp = amplitude
        
        if fitx:
            bgx = background(f0, linfit(Amp, para[0], para[1]),linfit(Amp, parb[0], parb[1]))
        else:
            bgx = background(f0,a2mean,b2mean)
            
        if fity:
            bgy = background(f0, linfit(Amp, para2[0], para2[1]),linfit(Amp, parb2[0], parb2[1]))
        else:
            bgy = background(f0,a1mean,b1mean)
            
        #print(bgx,bgy)
        #bgx, bgy = (0,0)
        R = np.sqrt((x - bgx)**2 + (y - bgy)**2) 

        amploss = bridgelosses(letter, f0, Amp).get_loss()        
        cirloss = circuitlosses(f0).get_loss()
        Amp = Amp*amploss*np.sqrt(2)
        
        R = R/cirloss*np.sqrt(2)
       
        P = Pgrad(bias, Amp,T)*1e-3
        
        # print(P)
        
        V = Velocity(bias,R,T)
        
        # ax4.plot(t,V,'-o',ms = 0.2,lw = 0.2, color = 'tab:blue')
     
        
        data_to_save = {'Pressure (Pa/m)': P*1e3,'Velocity (m/s)': V,'Temperature (K)': Temp,
                        'App_gen (V)':amplitude,'time (s)':t,'sample rate (Hz)':sample_rate}
        
        if save:
            path_save = fr'D:\OneDrive_copy\OneDrive - Univerzita Karlova\DATA\2023_4_Helmholtz_resonators_recal\Helmholtz_res_drive_dep\{folder}\buffer\resonator_{height:}{letter:}'
            os.makedirs(path_save,exist_ok=True)
            path_save = os.path.join(path_save,file.split('\\')[-1])
            np.save(path_save, data_to_save)
        
        
    # ax4.set_xlabel('pressure gradient (kPa/m)')
    # ax4.set_ylabel('velocity (m/s)')
    # ax4.set_title(f'T = {T:} K')
        

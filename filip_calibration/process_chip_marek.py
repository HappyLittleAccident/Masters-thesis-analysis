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
    'A':{        
        'w':  1000e-6,
        'l':  1050e-6,
        'C0': 338.63213961685864e-12,
        'kp': 1.2400632746323976e7
        },

    
    'B':{        
        'w':  500e-6,
        'l':  750e-6,
        'C0': 335.6978986554392e-12,
        'kp': 3*1.4546169638969617e7
        },
    
    'C':{        
        'w':  950e-6,
        'l':  1000e-6,
        'C0': 312.8732100355237e-12,
        'kp': 3*1.1909096177342793e7
        },
    'D':{
        'w': 500e-6, 
        'l': 950e-6,
        'C0':466.4158697787019e-12,
        'kp':3*1.6174502699593267e7
        }}

plt.close('all')
save = True

letters = ['A','B','C','D']
fig2,ax2 = plt.subplots(2,2,sharex=True)
fig,ax = plt.subplots(2,1,sharex=True)
 
fig3,ax3 = plt.subplots(1,1)
fig4,ax4 = plt.subplots(1,1)

for letter in letters[2:3]:
    print('letter: '+letter)
    for temp in temps[:]:
        ax[0].clear()
        ax[1].clear()
        
        T = float(temp[1:].split('_')[0])*1e-3
        print(f'\ntemp: {T}\n')

        def linfit(x,a,b):
            y = a*x + b
            return y
        
        def cubicfit(x,a,b,c,d):
            y = a*x**3 + b*x**2 + c*x + d
            return y
        
        def compres (T):
            a,b,c,d = ( 3.65456274e-10, -1.19112443e-09,  1.72705353e-09, 1.12838569e-08)
            chi = a*T**3 + b*T**2 + c*T + d
            return chi
          
        
        # plt.close('all')
        height = 500
        
        Nbot = 0
        Ntop = -2
        
        fitx = True
        fity = True
        
        resonator = f'{height:}{letter:}'
        
        ampsweeps = rf'D:\Github\Masters-thesis-analysis\filip_calibration\Helmholtz_res_drive_dep/{temp:}/ampsweeps_fast/resonator_{resonator}_updowndown/*.npy'        
        fsweeps   = rf'fsweeps_sorted_raw\{temp}\fsweeps\resonator_{height:}{letter:}\*.npy'
        
        Afiles = glob(ampsweeps)
        Ffiles = glob(fsweeps)
        
        '''background estimation'''
        par = []
        a1s = []
        b1s = []
        a2s = []
        b2s = []
        drs = []
        # fig,ax = plt.subplots(2,1)
        
        Afiles = np.sort(Afiles)
        Ffiles = np.sort(Ffiles)
        for num, file in enumerate(Ffiles[1:20]):
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
                       
            peak = MicrowavePeak(f, x, y, init_p)
            result = peak.fit(plot=True)
            par = result.best_values
            res_fit = result.best_fit
            

            ax[0].plot(f, x, 'o',ms = 3)
            ax[1].plot(f, y, 'o', ms = 3)
            ax[0].plot(f, np.real(res_fit), 'k--')
            ax[1].plot(f, np.imag(res_fit), 'k--')
            fig.canvas.draw()
            
            a1s.append(par['a1'])
            b1s.append(par['b1'])
            a2s.append(par['a2'])
            b2s.append(par['b2'])
            drs.append(amp)
                
        # fig2,ax2 = plt.subplots(2,2)
        for axis in np.ravel(ax2):
            axis.clear()
        
        ax2[0,0].plot(drs,a1s,'o')
        ax2[1,0].plot(drs,b1s,'o')
        ax2[0,1].plot(drs,a2s,'o')
        ax2[1,1].plot(drs,b2s,'o')

        
        def background (f,a,b):
            y = a*f + b
            return y
        
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
        
        fig2.canvas.draw()
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
            
        first = True
        for file in Afiles[Nbot:-2]: 
            d = np.load(file,allow_pickle=True).item()
            Amp_base = d['App_base (V)']
            Mod = np.array(d['App (V)'])
            x = np.array(d['x (A)'])
            y = np.array(d['y (A)'])
            depth = d['depth']
            f0 = d['freq (Hz)']
            bias = d['bias (V)']
            
            Amp = 0.5*Amp_base*(1 + depth*Mod/500)
            
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
            
            N = len(Amp)
            Ndiv = int(N/3)
            upAmp = Amp[:Ndiv]
            downAmp = Amp[Ndiv:2*Ndiv]
            jumpAmp = Amp[2*Ndiv:]
            upR = R[:Ndiv]
            downR = R[Ndiv:2*Ndiv]
            jumpR = R[2*Ndiv:]
            
            uploss   = bridgelosses(letter, f0, upAmp  ).get_loss()
            downloss = bridgelosses(letter, f0, downAmp).get_loss()
            jumploss = bridgelosses(letter, f0, jumpAmp).get_loss()
            
            cirloss = circuitlosses(f0).get_loss()
            
            upR = upR/cirloss*np.sqrt(2)
            downR = downR/cirloss*np.sqrt(2)
            jumpR = jumpR/cirloss*np.sqrt(2)
            

            
            upAmp = upAmp*uploss*np.sqrt(2)
            downAmp = downAmp*downloss*np.sqrt(2)
            jumpAmp = jumpAmp*jumploss*np.sqrt(2)
        
            
            upP = Pgrad(bias, upAmp,T)*1e-3
            downP = Pgrad(bias, downAmp,T)*1e-3
            jumpP = Pgrad(bias, jumpAmp,T)*1e-3
            
            
            
            upV = Velocity(bias,upR,T)
            downV = Velocity(bias,downR,T)
            jumpV = Velocity(bias,jumpR,T)
            
            
    
            if first:
                ax3.plot(upAmp,upR,'-o',ms = 0.2,lw = 0.2, color = 'tab:blue')
                ax3.plot(downAmp[:-5],downR[:-5],'-o',ms = 0.2,lw = 0.2, color = 'tab:orange')
                ax3.plot(jumpAmp[5:],jumpR[5:],'-o',ms = 0.2,lw = 0.2, color = 'tab:green')
                point3=ax3.plot(upAmp[-1],upR[-1],'ro')
                
                ax4.plot(upP,upV,'-o',ms = 0.2,lw = 0.2, color = 'tab:blue')
                ax4.plot(downP[:-5],downV[:-5],'-o',ms = 0.2,lw = 0.2, color = 'tab:orange')
                ax4.plot(jumpP[5:],jumpV[5:],'-o',ms = 0.2,lw = 0.2, color = 'tab:green')
                
                point = ax4.plot(upP[-1],upV[-1],'ro')
                first = False
            data_to_save = {'upPressure (Pa/m)': upP*1e3, 'downPressure (Pa/m)': downP*1e3, 'jumpPressure (Pa/m)': jumpP*1e3,
                            'upVelocity (m/s)': upV, 'downVelocity (m/s)': downV, 'jumpVelocity (m/s)': jumpV, 'Temperature (K)': T}
            
            fig3.canvas.draw()
            fig4.canvas.draw()
            if save:
                folder = temp
                name = file.split('\\')[-1]
                print(name)
                path_save = rf'D:\Github\Masters-thesis-analysis\processing_programs_plots\Helmholtz_res_drive_dep\{folder.split("_")[0]}\ampsweeps_fast\resonator_500{letter}_updowndown'
                os.makedirs(path_save,exist_ok=True)
                np.save(path_save + '\\' +name,data_to_save,allow_pickle=True)

        plt.figure()
        plt.plot(upAmp,bgx[:Ndiv])
        
        plt.plot(upAmp,bgy[:Ndiv])
        
        ax4.set_xlabel('pressure gradient (kPa/m)')
        ax4.set_ylabel('velocity (m/s)')
        ax4.set_title(f'T = {T:} K')
        plt.waitforbuttonpress()
        point[0].remove()
        point3[0].remove()
    ax3.clear()
    ax4.clear()

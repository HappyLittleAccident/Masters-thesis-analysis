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
        'kp': 1.4546169638969617e7
        },
    
    'C':{        
        'w':  950e-6,
        'l':  1000e-6,
        'C0': 312.8732100355237e-12,
        'kp': 1.1909096177342793e7
        },
    'D':{
        'w': 500e-6, 
        'l': 950e-6,
        'C0':466.4158697787019e-12,
        'kp':1.6174502699593267e7
        }}

plt.close('all')
save = False

letters = ['A','B','C','D']
fig2,ax2 = plt.subplots(2,2,sharex=True)
fig,ax = plt.subplots(2,1,sharex=True)
for letter in letters[2:3]:
    print('letter: '+letter)
    for temp in temps[:]:
        
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
            return chi*10
          
        
        height = 500
        Temp = temp
        
    
        fitx = True
        fity = True
        
        folder = temp
        
        resonator = f'{height:}{letter:}'
        
        
        fsweeps =       rf'Helmholtz_res_drive_dep\{temp}\fsweeps\resonator_{height:}{letter:}\*.npy'
        freqsweep_all = rf'Helmholtz_res_drive_dep\{temp}\fsweeps\resonator_{height:}{letter:}\*.npy'
                        
        
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
        
    
        
        
        # Ffiels = Ffiles.sort(key=lambda x : np.load(x,allow_pickle=True).item()['App (V)'])
        Freqfiles = np.sort(Freqfiles)
        
        for num, file in enumerate(Ffiles[1:10]):
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
    
    
            

        ax2[0,0].plot(drs,a1s,'o')
        ax2[1,0].plot(drs,b1s,'o')
        ax2[0,1].plot(drs,a2s,'o')
        ax2[1,1].plot(drs,b2s,'o')
        
        def background (f,a,b):
            y = a*(f-2000) + b
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
            print(g)
            v = I/(g*bias*C0)
            return v
            
        # fig3,ax3 = plt.subplots(1,1)
        # fig4,ax4 = plt.subplots(1,1)
        for file in Freqfiles: 
            d = np.load(file,allow_pickle=True).item()
            amplitude = np.array(d['App (V)'])
            
            f = np.array(d['freq (Hz)'])
            x = np.array(d['x (A)'])
            y = np.array(d['y (A)'])
            bias = d['bias (V)']
            
            
            Amp = amplitude
            
            if fitx:
                bgx = background(f, linfit(Amp, para[0], para[1]),linfit(Amp, parb[0], parb[1]))
            else:
                bgx = background(f,a2mean,b2mean)
                
            if fity:
                bgy = background(f, linfit(Amp, para2[0], para2[1]),linfit(Amp, parb2[0], parb2[1]))
            else:
                bgy = background(f,a1mean,b1mean)
                
            #print(bgx,bgy)
            #bgx, bgy = (0,0)
            
            R = np.sqrt((x - bgx)**2 + (y - bgy)**2) 
            
            V = np.zeros_like(x)
            for i,f0 in enumerate(f):
                        
                cirloss = circuitlosses(f0).get_loss()
                
                R[i] = R[i]/cirloss*np.sqrt(2) 
                v = Velocity(bias,R[i],T)[0]
    
                V[i] = v
                # ax4.plot(t,V,'-o',ms = 0.2,lw = 0.2, color = 'tab:blue')
            amploss = bridgelosses(letter, f[np.argmax(V)], Amp).get_loss()
            print(amploss)
            Amp = Amp*amploss*np.sqrt(2)
            P = Pgrad(bias, Amp,T)*1e-3    
            # print(f'pressure={P:.2f}')
            data_to_save = {'Pressure (Pa/m)': P*1e3,'Velocity (m/s)': V,'Temperature (K)': Temp,'Freq (Hz)':f,
                            'App_gen (V)':amplitude}
            
            if save:
                path_save = fr'2023_4_Helmholtz_resonators_recal_2\Helmholtz_fsweeps\{folder.split("_")[0]}\{height:}{letter:}'
                os.makedirs(path_save,exist_ok=True)
                path_save = os.path.join(path_save,file.split('\\')[-1])
                np.save(path_save, data_to_save)
        
        
        ax2[0,0].set_title('lin y')
        ax2[1,0].set_title('con y')
        ax2[0,1].set_title('lin x')
        ax2[1,1].set_title('con x')
        
        fig2.canvas.draw()
        fig.canvas.draw()
        
        if not save:
            while not plt.waitforbuttonpress():
                pass
        ax[0].clear()
        ax[1].clear()
        for axis in np.ravel(ax2):
            axis.clear()
        
    
            
        # ax4.set_xlabel('pressure gradient (kPa/m)')
        # ax4.set_ylabel('velocity (m/s)')
        # ax4.set_title(f'T = {T:} K')
            

# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 23:32:51 2024

@author: Marek
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 19:40:09 2024

@author: Marek
"""

import numpy as np
import matplotlib.pyplot as plt
import helium_old as he
from glob import glob
import os
from scipy.interpolate import CubicSpline
import sympy
from sympy import I, init_printing

plt.close('all')

#i love this line of code
# k,A,rho,A,D,chi,alpha,a,rhos,rhon,T0,cp,R,s,l,B,kappa,L,eta,omega = sympy.symbols('k A rho A D chi alpha a rhos rhon T0 cp R s l B kappa L eta omega',positive = True)
omega = sympy.symbols('omega')
F,x,vs,vn,dT,dP = sympy.symbols('F x vs vn dT dP',complex = True)

def compres (T):
    a,b,c,d = ( 3.65456274e-10, -1.19112443e-09,  1.72705353e-09, 1.12838569e-08)
    chi = a*T**3 + b*T**2 + c*T + d
    return chi*10


def res_freq_simple(k,rhos,l,a,rho,A,D,chi,T0,c,s):
    return np.sqrt(2*rhos*a/(l*rho**2) *k/(2*A**2 + k*D*A*chi))/2/np.pi 
                   

T_spline = [0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0]
alpha_spline = np.array([4.5, 3.3, 0.8, -3.2, -9.5, -18.0, -26.8, -38.0, -52.2, -71.3, -97.6, -129.6])*1e-4
alpha_function = CubicSpline(T_spline,alpha_spline)

const_a = 1.6
const_b = 1.12
const_c = 1.49
const_d = 1.79
ka = 3
kb = 3
kd = 3
const_k = 1
cut = 0


resonators_geometry = {
    'A':{        
        'w':  1000e-6,
        'l':  const_a*(1050e-6),
        'C0': 338.63213961685864e-12,
        'kp': ka*1.2400632746323976e7,
        'D':  500e-9
        },

    
    'B':{        
        'w':  500e-6,
        'l':  const_b*(750e-6-cut),
        'C0': 335.6978986554392e-12,
        'kp': kb*1.4546169638969617e7,
        'D':  500e-9
        },
    
    'C':{        
        'w':  950e-6,
        'l':  const_c*(1000e-6-cut),
        'C0': 312.8732100355237e-12,
        'kp': 3*1.1909096177342793e7,
        'D':  500e-9
        },
    'D':{
        'w': 1000e-6,#500e-6,
        'l': const_d*(950e-6-cut),
        'C0':466.4158697787019e-12,
        'kp':kd*1.6174502699593267e7,
        'D':  500e-9
        }}


A = np.pi*(2.5e-3)**2
colors = ['k','b','r','tab:orange','tab:green','tab:pink','m','tab:grey']

folder = r'D:\OneDrive_copy\OneDrive - Univerzita Karlova\DATA\2023_4_Helmholtz_resonators_recal_2\resonance_freq'
files = glob(os.path.join(folder,'*.npy'))


figres,axres = plt.subplots(1,1,sharex=True)    
figrest,((axvn,axx),(axT,axP)) = plt.subplots(2,2,sharex=True)

figdamp,axdamp = plt.subplots(1,1,sharex=True) 
figgraph,axgraph = plt.subplots(2,1,sharex=True)
cmap = plt.get_cmap('viridis')

options = [1,1,1,1,1,1,1e-25,1,1,1,1,1,1,1,1]
labels = ['no approx','alpha=0','vn T grad zero','vn P grad zero','dot vn zero','entropy dP zero','mass dT zero','full']

for i in range(2):

    f0_theory = []
    w_theory=[]
    f0_simple = []

    # fig, ax = plt.subplots()
    letter = files[0].split('\\')[-1].split('.')[0][-1]
    
    tempstop = -1
    d = np.load(files[0],allow_pickle=True).item()
    temperatures = d['temperature (K)'][0:tempstop]
    f0s = d['frequency (Hz)'][0:tempstop]
    f0 = [np.mean(temp) for temp in f0s]
    ws = d['width (Hz)'][0:tempstop]
    w = [np.mean(temp) for temp in ws]

    # ax.plot(temperatures,f0,'o',c=colors[i],label='Measured')
    

    print(int(i==0))
    for T in temperatures:
        resonator = resonators_geometry[letter]
        C0 = resonator['C0']
        k = resonator['kp']
        l = resonator['l']
        D = resonator['D']
        a = D*resonator['w']
        rhos = he.superfluid_density(T)[0]
        rho = he.density(T)[0]
        rhon = rho-rhos
        chi = compres(T)
        T0 = T
        cp = he.heat_capacity(T)[0]/0.0040026
        s = he.entropy(T)[0]
        R = 550*T**(-3.6) #60e-4*T**(-3) #
        B = he.mutual_friction_B(T)[0]
        kappa = he.quantized_circulation
        L = 0e9
        alpha = alpha_function(T).item()
        eta = he.viscosity(T)[0]
        
        c1 = he.first_sound_velocity(T)[0]
        c2 = he.second_sound_velocity(T)[0]
        c4 = np.sqrt(c1**2*rhos/rho + c2**2 * (1-rhos)/rho)
    
        eqnF = sympy.Eq(x*k + dP*A,F)
        
        eqnM = sympy.Eq(I*omega*rho*A*D*(chi*dP - alpha*dT*int(i==0))-2*I*omega*x*A*rho,
                        -2*a*(rhos*vs + 2*rhon*vn/3))
        eqnS = sympy.Eq(I*omega*(dT-T0*alpha*dP/(cp*rho)*int(i==0)),
                        (-dT +2*s*a*T0*R*rhos*(vs-2*vn/3*int(i==0)))/(R*rho*D*A*cp))
        eqnVs = sympy.Eq(I*omega*rhos*vs,
                         rhos*dP/(rho*l) - rhos*s*dT/l - B*kappa*rhon*rhos*L*(vs-2*vn/3*int(i==0))/rho)
        eqnVn = sympy.Eq(I*omega*rhon*vn*int(i==0),
                         rhon*dP/(rho*l) + rhos*s*dT/l*int(i==0) + B*kappa*rhon*rhos*L*(vs-2*vn/3*int(i==0))/rho - eta*8*vn/(D**2))

        solution = sympy.solve([eqnF,eqnM,eqnS,eqnVs,eqnVn],x,dP,dT,vs,vn,dict = True)

        
        f0_basic = res_freq_simple(k,rhos,l,a,rho,A,D,chi,T0,cp,s)
        freqs = 2*np.pi*np.linspace(f0_basic-400,f0_basic + 400,10000)
        vs_plot =sympy.lambdify(omega,sympy.simplify(solution[0][vs]/F))
        vn_plot =sympy.lambdify(omega,sympy.simplify(solution[0][vn]/F))
        T_plot =sympy.lambdify(omega,sympy.simplify(solution[0][dT]/F))
        P_plot =sympy.lambdify(omega,sympy.simplify(solution[0][dP]/F))
        x_plot =sympy.lambdify(omega,sympy.simplify(solution[0][x]/F))
        
        w_theory.append(np.abs(freqs[np.argmax(np.imag(vs_plot(freqs)))]-
                               freqs[np.argmin(np.imag(vs_plot(freqs)))])/2/np.pi)
            
        f0_theory.append(freqs[np.argmax(vs_plot(freqs))]/2/np.pi)
        f0_simple.append(res_freq_simple(k,rhos,l,a,rho,A,D,chi,T0,cp,s))

        
        Fes = (C0*10*8e-3)/D
        
        color_scale = (T-temperatures[0])/(temperatures[-1]-temperatures[0])     
        axgraph[0].plot(freqs/np.pi/2,np.real(vs_plot(freqs)*Fes),c=cmap(color_scale))
        axgraph[1].plot(freqs/np.pi/2,np.imag(vs_plot(freqs)*Fes),c=cmap(color_scale))
        
        axx.plot(freqs/np.pi/2,np.abs(x_plot(freqs)*Fes),c=cmap(color_scale))
        axvn.plot(freqs/np.pi/2,np.abs(vn_plot(freqs)*Fes),c=cmap(color_scale))
        axP.plot(freqs/np.pi/2,np.abs(P_plot(freqs)*Fes),c=cmap(color_scale))
        axP.axhline(Fes*2*A/(2*A**2 + k*D*A*chi))
        axT.plot(freqs/np.pi/2,np.abs(T_plot(freqs)*Fes),c=cmap(color_scale))
        # print('c4: ',c4)
    print(letter + ': ',f0[0],'\ntheory: ',f0_theory[0], '\ntheory simple: ',f0_simple[0])

    axx.set_title('x')
    axvn.set_title('vn')
    axT.set_title('dT')
    axP.set_title('dP')
    figrest.tight_layout()
    # ax.plot(temperatures,f0_theory,'k--')
    

    axres.plot(temperatures,f0_theory,'--',c=colors[i],label='Theory ' + labels[i])    
    
    axdamp.plot(temperatures,w_theory,'--',c=colors[i],label='Theory ' + labels[i])
    
# axres.plot(temperatures,f0,'o',c='k',label='Measured '+letter)
axres.set_xlabel('Temperature (K)')
axres.set_ylabel('Resonance frequency (Hz)')
axres.legend()  
    
# axdamp.plot(temperatures,w,'o',c='k',label='Measured '+letter)  
axdamp.set_xlabel('Temperature (K)')
axdamp.set_ylabel('Linewidth (Hz)')
axdamp.legend()


    
    










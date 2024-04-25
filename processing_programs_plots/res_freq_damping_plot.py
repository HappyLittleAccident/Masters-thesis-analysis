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


def res_freq(k,rhos,l,a,rho,A,D,chi,T0,c,s):
    return np.sqrt(2*rhos*a/(l*rho**2) *k/(2*A**2 + k*D*A*chi) 
                   +2*a*s**2 * T0*rhos/(l*c*rho*A*D))/2/np.pi
def width(a,s,T0,rhos,R,l,B,kappa,rho,L):
    return (2*a*s**2 * T0 * rhos * R/l + B*kappa*(rho-rhos)*L/rho)/np.pi/2 

def res_freq_full(omega0,tau,OMEGA):
    return np.sqrt((
                    omega0**2 *tau**2 + OMEGA*tau -1 +
                    np.sqrt(4*omega0**2 * tau**2 + (1-omega0**2*tau**2-OMEGA*tau)**2)
                    )/2
                   )/tau/2/np.pi

def width_full(omega0,tau,OMEGA,B,kappa,rho,rhos,L):
    return (2*OMEGA/(
                    omega0**2 *tau**2 + OMEGA*tau + 1 +
                    np.sqrt(4*omega0**2 * tau**2 + (1-omega0**2*tau**2-OMEGA*tau)**2)
                    )
        + B*kappa*(rho-rhos)*L/rho)/2/np.pi

def response(omega,omega0,tau,OMEGA):
    return 1j*omega/(-omega**2+omega0**2 + 1j*omega*OMEGA/(1j*tau * omega +1))

D = 500e-9
A = np.pi*(2.5e-3)**2
colors = ['r','b','tab:orange','k']
for i,file in enumerate(files[:]):
    letter = file.split('\\')[-1].split('.')[0][-1]

    d = np.load(file,allow_pickle=True).item()
    T = np.array(d['temperature (K)'])

    f0s = d['frequency (Hz)']
    f0 = [np.mean(x) for x in f0s]
    
    ws = d['width (Hz)']
    w = [np.mean(x) for x in ws]
    
    resonator = resonators_geometry[letter]
    k = resonator['kp']
    l = resonator['l']
    a = D*resonator['w']
    
    rhos = he.superfluid_density(T)
    rho = he.density(T)
    chi = compres(T)
    T0 = T
    c = he.heat_capacity(T)/0.0040026
    s = he.entropy(T)
    R = 300e-4*T**(-3)/A #89.1*T**(-3.6)
    B = he.mutual_friction_B(T)
    kappa = he.quantized_circulation
    L = 0e9
    c4 = he.fourth_sound_velocity(T)
    alpha = expand(T)
    
    tau = (R*c*D*A*rho)
    OMEGA = 2*s**2 * a * T0*R*rhos/l
    omega0 = np.sqrt(2*rhos*a/(l*rho**2) *k/(2*A**2 + k*D*A*chi))
    
    csi = 2e-3*1.5**3
    rhosi = 2634
    # tau = s/s*1e-4

    er = 1.0565
    sigma = chi*D*k/(2*A)
    g = 2*a*rhos/(D*A*rho*(1 + 2*sigma))
    print(sigma)
    # print((2*chi/(A*(1+sigma)))*(resonator['C0']**2)*30*30*1*f0*2*np.pi/D)
    print(g)
    print(30**2*resonator['C0']*chi/(D*A*(1 + 2*sigma)))
    print(4e-9/g/30/resonator['C0'])
    # print(np.array([csi*rhosi*1*9*7e-9*R,tau]).T)
    # print(np.array([1*tau/tau,tau*OMEGA,tau**2*omega0**2]).T)
    # print(OMEGA*tau)
    # print('Heat expansion \t entropy flow \t SiO2 heatflow')
    # print(np.array([omega0**2 * alpha*rho*l*A*D/(2*rhos*a*c),2*s*a*c4*rhos*R*l/c4/tau,l/c4/tau]).T)
    # print(np.array([2*a*rhos*l/D/rho/A,alpha*T0]).T)
    
    # print((T0*2*s*a*rhos*1e-1)/(rho*A*D*c*1e4))
    
    # print([chi,alpha*1e-4])

    
    # theory_res_freq = res_freq(k,rhos,l,a,rho,A,D,chi,T0,c,s)
    theory_res_freq = res_freq_full(omega0,tau,OMEGA)
    pure_theory= res_freq(k,rhos,l,a,rho,A,D,chi,T0,c,0)
    # print(theory_res_freq-omega0/2/np.pi)
    
    # theory_width = width(a,s,T0,rhos,R,l,B,kappa,rho,L)
    theory_width = width_full(omega0,tau,OMEGA,B,kappa,rho,rhos,L)
    

    
    axres.plot(T,f0,'o',c=colors[i],label='Measured')
    axres.plot(T,theory_res_freq,'--',c=colors[i],label='Theoretical prediction corrected')
    axres.plot(T,pure_theory,'-.',lw=0.5,c=colors[i],label='Theoretical prediction')
    
    axdamp.plot(T,w,'o',c=colors[i],label='Measured')
    axdamp.plot(T,theory_width,'--',c=colors[i],label='Theoretical prediction')
    
    omega0_graph = []
    damping_graph = []
    omega = 2*np.pi*np.linspace(1400,2400,5000)
    cmap = plt.get_cmap('viridis')
    for j in range(len(T)):
        res_curve = response(omega, omega0[j], tau[j], OMEGA[j])
        omega0_graph.append(omega[np.argmax(np.real(res_curve))]/2/np.pi)
        damping_graph.append(np.abs(omega[np.argmax(np.imag(res_curve))]-omega[np.argmin(np.imag(res_curve))])/2/np.pi)
        
        color_scale = (T[j]-T[0])/(T[-1]-T[0])        
        axgraph[0].plot(omega/np.pi/2,np.real(res_curve),c=cmap(color_scale))
        axgraph[1].plot(omega/np.pi/2,np.imag(res_curve),c=cmap(color_scale))
        figgraph.canvas.draw()
        
    # axres.plot(T,omega0_graph,'-',c=colors[i])
    # axdamp.plot(T,damping_graph,'-',c=colors[i])
    axres.set_xlabel('Temperature (K)')
    axres.set_ylabel('Resonance frequency (Hz)')
    axres.legend()
    
    
    axdamp.set_xlabel('Temperature (K)')
    axdamp.set_ylabel('Linewidth (Hz)')
    axdamp.legend()
    
    figres.canvas.draw()
    figdamp.canvas.draw()
    print(np.log(3e-9/1e-10)*kappa/2/np.pi/D)
 
    break
    while not plt.waitforbuttonpress():
        pass
    axgraph[0].clear()
    axgraph[1].clear()
    axdamp.clear()
    axres.clear()

    
    
    
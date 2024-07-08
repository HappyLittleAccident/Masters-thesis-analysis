# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 16:55:10 2023

@author: Marek
"""

import numpy as np
import matplotlib.pyplot as plt
import os
from glob import glob
from format_figure import use_latex,polish
import matplotlib as mpl
import helium_old as he
from scipy.interpolate import CubicSpline


plt.close('all')

temps = ['T1325', 'T1350', 'T1375', 'T1400', 'T1425', 'T1450', 'T1475',
         'T1500', 'T1550', 'T1650', 'T1700', 'T1750', 'T1850', 'T1950']
resonators = ['500A', '500B', '500C','500D']
names = ['C', 'J', 'R','G']

use_latex()
# plt.rcdefaults()
plt.ion()

cmap = plt.get_cmap('winter')
viridis = plt.get_cmap('viridis')

# plt.ioff()
# plt.ion()
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
        'kp': 3*1.4546169638969617e7
        },
    
    'C':{        
        'w':  950e-6,
        'l':  1.03*1000e-6,
        'C0': 312.8732100355237e-12,
        'kp': 3*1.1909096177342793e7
        },
    'D':{
        'w': 1000e-6, #orig= 500e-6
        'l': 1.405*950e-6,
        'C0':466.4158697787019e-12,
        'kp':3*1.6174502699593267e7
        }}

D = 500e-9
A = np.pi*(2.5e-3)**2
colors = ['r','b','tab:orange','k']
plt.close('all')
velocity_corrections= [1.7,1.3,7]
fig,ax = plt.subplots(1,1,sharex=False,dpi=150)
for i,resonator in enumerate(resonators[0:2]):

    fig_w, ax_w = plt.subplots()
    for j,temp in enumerate(temps):
        scale = (float(temp[1:])*1e-3 - 1.325)/(1.95 - 1.325)
        

        filenames = glob(
            fr'Helmholtz_res_drive_dep\{temp}\ampsweeps_fast\resonator_{resonator}_updowndown\*.npy'
            )
        filenames_vld = glob(
            rf'D:\OneDrive_copy\OneDrive - Univerzita Karlova\DATA\2023_4_Helmholtz_resonators_recal_2\Helmholtz_res_VLD\{temp}\{resonator}\*.npy'
            )

        try:
            data = np.load(filenames[10],allow_pickle=True).item()
            data_vld = np.load(filenames_vld[10],allow_pickle=True).item()
        except:
            pass
        if i == 2 and float(temp[1:]) >1800:
            continue
        
        vs_treshold = 0.15*(i==0) + 0.095*(i==1) + 0.02*(i==2)
        friction_threshold = 20 #150 - 135*(i==1) - (i==2)*75
        
        upP = (data['upPressure (Pa/m)']/1.5)[data_vld['upV (m/s)']>vs_treshold]
        downP = data['downPressure (Pa/m)']
        jumpP = data['jumpPressure (Pa/m)']
        up = data['upVelocity (m/s)']*100
        down =  data['downVelocity (m/s)']*100
        jump = data['jumpVelocity (m/s)']*100
        
        
        upV =(data_vld['upV (m/s)']*velocity_corrections[i])[data_vld['upV (m/s)']>vs_treshold]
        upL = (data_vld['upL (m^-2)'])[data_vld['upV (m/s)']>vs_treshold]
        T = data_vld['Temp (K)']
        width = data_vld['Width (Hz)']
        letter = resonator[-1]

        
        
        resonator_pars = resonators_geometry[letter]
        k = resonator_pars['kp']
        l = resonator_pars['l']
        a = D*resonator_pars['w']

        
        rhos = he.superfluid_density(T)
        rho = he.density(T)
        chi = compres(T)
        T0 = T
        c = he.heat_capacity(T)/0.0040026
        s = he.entropy(T)
        B = he.mutual_friction_B(T)
        alpha = he.mutual_friction_alpha(T)
        kappa = he.quantized_circulation
        L = 0e9
        c4 = he.fourth_sound_velocity(T)
        # alpha = expand(T)
        
        

        ax_w.plot(upP/rho/upV/2/np.pi)
        ax_w.axhline(width)
        
        friction_scale = 1e3*upV/np.sqrt(2)/(alpha*upL*kappa)
        friction_scale_alternative = 1e3*upP/(np.sqrt(2)*rho*(2*np.pi*width+alpha*upL*kappa)*alpha*upL*kappa)
        # ax.plot((upV*100)[np.logical_and(friction_scale<5,friction_scale>0)], friction_scale[np.logical_and(friction_scale<5,friction_scale>0)],c=cmap(scale))
        ax.plot(
            (upV*100)[np.logical_and(friction_scale<friction_threshold,friction_scale>0)],
            friction_scale_alternative[np.logical_and(friction_scale<friction_threshold,friction_scale>0)]
            ,c=cmap(scale),marker='.' if i==0 else '.',ls='',ms=3)
        
        
        # print(T, (upP/rho/upV/2/np.pi) - width)
        # ax.plot(upV*100,upL + indent,c=cmap(scale))
t = np.linspace(25.1,55,1000)
ax.set_xlim(5,60)
ax.plot(t,1/(t-25) +0.011*(t-55)**2 +1,ls=':',c='k')
ax.annotate('Resonator J',xy=(0.02,0.25),xycoords='axes fraction')
ax.annotate('Resonator C',xy=(0.7,0.5),xycoords='axes fraction')
ax.axhline(1,c='k',ls = '--',label='channel dimension') #1mm
ax.set_xlabel('Superfluid velocity (cm/s)')
ax.set_ylabel(r'Friction scale $\mathcal{L}_{\mu}$ (mm)')
ax.legend()
fig.colorbar(plt.cm.ScalarMappable(norm=mpl.colors.Normalize(
vmin=1.325, vmax=1.95, clip=False), cmap=cmap), ax=ax, label='Temperature (K)')
# fig.colorbar(plt.cm.ScalarMappable(norm=mpl.colors.Normalize(
# vmin=1.325, vmax=1.95, clip=False), cmap=viridis), ax=ax, label='Temperature (K)',location='left')


# ax.set_ylabel(r'VLD $(m^{-2})$')
polish(fig, 1, name=f'images//friction{names[i]}', extension='.png', grid=True,width_to_height = 1.5,tight_layout=True)        

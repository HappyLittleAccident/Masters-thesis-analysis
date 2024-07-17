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
from format_figure import use_latex,polish
import matplotlib as mpl

plt.close('all')
# plt.rcdefaults()
use_latex()

#i love this line of code
# k,A,rho,A,D,chi,alpha,a,rhos,rhon,T0,cp,R,s,l,B,kappa,L,eta,omega = sympy.symbols('k A rho A D chi alpha a rhos rhon T0 cp R s l B kappa L eta omega',positive = True)
omega = sympy.symbols('omega')
F,x,vs,vn,dT,dP,ell = sympy.symbols('F x vs vn dT dP l',complex = True)

def compres (T):
    a,b,c,d = ( 3.65456274e-10, -1.19112443e-09,  1.72705353e-09, 1.12838569e-08)
    chi = a*T**3 + b*T**2 + c*T + d
    return chi*10


def res_freq_simple(k,rhos,l,a,rho,A,D,chi,T0,c,s):
    return np.sqrt(2*rhos*a/(l*rho**2) *k/(2*A**2 + k*D*A*chi))/2/np.pi 

def res_freq_approximate(k,rhos,l,a,rho,A,D,chi,T0,c,s,R):
    return np.sqrt(2*rhos*a/(l*rho**2) *k/(2*A**2 + k*D*A*chi) +2*a*s**2*rhos*T0/(c*A*D*rho*l))/2/np.pi       

def res_freq_full(omega0,tau,OMEGA):
    return np.sqrt((
                    omega0**2 *tau**2 + OMEGA*tau -1 +
                    np.sqrt(4*omega0**2 * tau**2 + (1-omega0**2*tau**2-OMEGA*tau)**2)
                    )/2
                   )/tau/2/np.pi
def width(a,s,T0,rhos,R,l,B,kappa,rho,L,omega0,tau,eta):
    return (2*a*s**2 * T0 * rhos * R/(omega0**2 * tau**2 + 1)/l + (rho-rhos)**2*D**2*omega0**2/12/eta/rhos)/np.pi/2 



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
colors = ['tab:blue','tab:orange','tab:green']
markers = ['o','x','+']

folder = r'D:\OneDrive_copy\OneDrive - Univerzita Karlova\DATA\2023_4_Helmholtz_resonators_recal_2\resonance_freq'
files = glob(os.path.join(folder,'*.npy'))


fig,(axres,axdamp) = plt.subplots(2,1,dpi=400,sharex = True)    
#%%
figrest,((axvn,axx),(axT,axP),(dummy1,dummy2)) = plt.subplots(3,2,sharex=True,dpi=400,height_ratios = [1,1,0.1])

gs = axx.get_gridspec()
dummy1.remove()
dummy2.remove()
axcbar = figrest.add_subplot(gs[2,:])  


fig_t,(axres_t,axdamp_t) = plt.subplots(2,1,dpi=400,sharex = True)    


#%%
figgraph,axgraph = plt.subplots(2,2,sharex=True,dpi=400,width_ratios=[1,0.05])

gs = axgraph[0,0].get_gridspec()
for a in axgraph[0:,1]:
    a.remove()
axcbargraph = figgraph.add_subplot(gs[:, 1])  
#%%


cmap = plt.get_cmap('viridis')
inferno = plt.get_cmap('Reds')

# fig_meas,ax_meas = plt.subplots()
plot_freqs = []
plot_widths = []
for i,file in enumerate(files[:1]):

    # fig, ax = plt.subplots()
    letter = file.split('\\')[-1].split('.')[0][-1]
    
    tempstop = -1
    d = np.load(file,allow_pickle=True).item()
    temperatures = d['temperature (K)'][0:tempstop]
    f0s = d['frequency (Hz)'][0:tempstop]
    f0 = [np.mean(temp) for temp in f0s]
    ws = d['width (Hz)'][0:tempstop]
    w = [np.mean(temp) for temp in ws]

    # ax.plot(temperatures,f0,'o',c=colors[i],label='Measured')

    Rs = [89.1,550]
    linestyles = ['--','-']
    for R0,linestyle in zip(Rs,linestyles):
        w_theory = []
        w_simple = []
        f0_theory = []
        f0_simple = []
        for T in temperatures:

            T = T
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
            # R =  100*81.1*T**(-3.6) #60e-4*T**(-3) #
            R =  R0*T**(-3.6)
            B = he.mutual_friction_B(T)[0]
            kappa = he.quantized_circulation
            L = 0e9
            alpha = alpha_function(T).item()
            eta = he.viscosity(T)[0]
            
            
            c1 = he.first_sound_velocity(T)[0]
            c2 = he.second_sound_velocity(T)[0]
            c4 = np.sqrt(c1**2*rhos/rho + c2**2 * (1-rhos)/rho)
        
            eqnF = sympy.Eq(x*k + dP*A,F)
            
            eqnM = sympy.Eq(I*omega*rho*A*D*(chi*dP - alpha*dT)-2*I*omega*x*A*rho,
                            -2*a*(rhos*vs + 2*rhon*vn/3))
            eqnS = sympy.Eq(I*omega*(dT-T0*alpha*dP/(cp*rho)),
                            (-dT +2*s*a*T0*R*rhos*(vs-2*vn/3))/(R*rho*D*A*cp))
            eqnVs = sympy.Eq(I*omega*rhos*vs,
                             rhos*dP/(rho*l) - rhos*s*dT/l - B*kappa*rhon*rhos*L*(vs-2*vn/3)/rho)
            eqnVn = sympy.Eq(I*omega*rhon*vn,
                             rhon*dP/(rho*l) + rhos*s*dT/l + B*kappa*rhon*rhos*L*(vs-2*vn/3)/rho - eta*8*vn/(D**2))
            
            print('solving')
            solution = sympy.solve([eqnF,eqnM,eqnS,eqnVs,eqnVn],x,dP,dT,vs,vn,ell,dict = True)
            print('solved')
            
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
            tau = (R*cp*D*A*rho)
            OMEGA = 2*s**2 * a * T0*R*rhos/l
            omega0 = 2*np.pi*f0_basic
            f0_simple.append(res_freq_approximate(k,rhos,l,a,rho,A,D,chi,T0,cp,s,R))
            w_simple.append(width(a,s,T0,rhos,R,l,B,kappa,rho,L,omega0,tau,eta))


            # f0_simple.append(res_freq_simple(k,rhos,l,a,rho,A,D,chi,T0,cp,s))
            # geom_factor_simple.append(f0_basic/np.sqrt(rhos/(chi*rho**2)))
            
            Fes = (C0*30*8e-3)/D/4.647
            
            
            color_scale = (T-temperatures[0])/(temperatures[-1]-temperatures[0]) 
            # axgraph[0].clear()
            # axgraph[1].clear()
            if R0==550:
                axgraph[0,0].plot(freqs/np.pi/2,np.real(vs_plot(freqs)*Fes)*100,c=cmap(color_scale))
                axgraph[1,0].plot(freqs/np.pi/2,np.imag(vs_plot(freqs)*Fes)*100,c=cmap(color_scale))
                # axgraph[0].plot(freqs/np.pi/2,np.real(vs_plot(freqs)*Fes + vn_plot(freqs)*Fes*rhon/rhos),c=inferno(color_scale))
            # axgraph[1].plot(freqs/np.pi/2,np.imag(vs_plot(freqs)*Fes+ vn_plot(freqs)*Fes*rhon/rhos),c=inferno(color_scale))
            if R0==550:
                axx.plot(freqs/np.pi/2,np.abs(x_plot(freqs)*Fes*1e9),c=cmap(color_scale))
                axvn.plot(freqs/np.pi/2,np.abs(vn_plot(freqs)*Fes)*100,c=cmap(color_scale))
                axP.plot(freqs/np.pi/2,np.abs(P_plot(freqs)*Fes),c=cmap(color_scale))
                axT.plot(freqs/np.pi/2,np.abs(T_plot(freqs)*Fes*1e6),c=cmap(color_scale))
            # figgraph.canvas.draw()
            # print('c4: ',c4)
            # while not plt.waitforbuttonpress():
            #     pass

        # print(letter + ': ',f0[0],'\ntheory: ',f0_theory[0], '\ntheory simple: ',f0_simple[0])
        
        
        # axx.set_title('x')
        # axvn.set_title('vn')
        # axT.set_title('dT')
        # axP.set_title('dP')
        # figrest.tight_layout()
        # # ax.plot(temperatures,f0_theory,'k--')
    
        axres.plot(np.array(temperatures),f0_theory,linestyle,c=colors[i],label=fr'$R_0={R0:.1f}$'+'$\mathrm{J^{-1} s^{-1}}$')
        # 
        axdamp.plot(temperatures,w_theory,linestyle,c=colors[i],label=fr'$R_0={R0:.1f}$'+'$\mathrm{J^{-1} s^{-1}}$')
        if R0 == 550 and letter=='A':
            axres_t.plot(np.array(temperatures),f0_theory, '-',c='tab:orange',lw=0.8  ,label='Numeric solution')
            axres_t.plot(np.array(temperatures),f0_simple,ls=(0, (5, 10)),c='k',lw=0.8,label='Approximate prediction')
            
            axdamp_t.plot(np.array(temperatures),w_theory, '-',c='tab:orange',lw=0.8  ,label='Numeric solution')
            axdamp_t.plot(np.array(temperatures),w_simple,ls=(0, (5, 10)),c='k',lw=0.8,label='Approximate prediction')
            
    dotres,=axres.plot(np.array(temperatures),f0,markers[i],c=colors[i])
    plot_freqs.append(dotres)
    # axres.legend()
    
    dotwidth,=axdamp.plot(temperatures,w,markers[i],c=colors[i])
    plot_widths.append(dotwidth)
    
    

# axres.set_xlabel('Temperature (K)')
#%%
print('Grad P =',Fes*2/(2*A + k*D*chi)/l*1e-3)

axres.set_ylabel('Resonance frequency (Hz)')
axdamp.set_xlabel('Temperature (K)')
axdamp.set_ylabel('Linewidth (Hz)')
    
linedamp1, = axdamp.plot([],[],'--k')
linedamp2, = axdamp.plot([],[],'-k')

axdamp_legend_numerics = axdamp.legend([linedamp1,linedamp2],[r'$R_0=89.1$'+r' $\mathrm{K J^{-1} s^{-1}}$' ,r'$R_0=550$'+r' $\mathrm{K J^{-1} s^{-1}}$'],loc='upper left')
axdamp.legend(plot_freqs,['C','J','R'],loc='lower right')
axdamp.add_artist(axdamp_legend_numerics)

lineres1, = axres.plot([],[],'--k')
lineres2, = axres.plot([],[],'-k')
axres_legend_numerics = axres.legend([lineres1,lineres2],[r'$R_0=89.1$'+r' $\mathrm{K J^{-1} s^{-1}}$' ,r'$R_0=550$'+r' $\mathrm{K J^{-1} s^{-1}}$'],loc='lower left')
axres.legend(plot_freqs,['C','J','R'],loc='upper right')
axres.add_artist(axres_legend_numerics)
# fig.tight_layout()

polish(fig, 1, name='images//numerics', extension='.png', grid=True,width_to_height = 0.6,tight_layout=True)        

  
axres_t.set_ylabel('Resonance frequency (Hz)')
axdamp_t.set_xlabel('Temperature (K)')
axdamp_t.set_ylabel('Linewidth (Hz)') 
axdamp_t.legend()
axres_t.legend() 
#%%
polish(fig_t, 0.7, name='images//numerics_theory', extension='.png', grid=True,width_to_height = 0.6,tight_layout=True)        


axgraph[1,0].set_xlabel('Frequency (Hz)')
axgraph[0,0].set_ylabel('In-phase $v_s$ (cm/s)')
axgraph[1,0].set_ylabel('Out-of-phase $v_s$ (cm/s)')
figgraph.colorbar(plt.cm.ScalarMappable(norm=mpl.colors.Normalize(
vmin=1.325, vmax=1.95, clip=False), cmap=cmap), ax=None,cax=axcbargraph, label='Temperature (K)',location = 'right')
figrest.colorbar(plt.cm.ScalarMappable(norm=mpl.colors.Normalize(
vmin=1.325, vmax=1.95, clip=False), cmap=cmap), ax=None,cax=axcbar, label='Temperature (K)',location = 'bottom')


axP.set_xlabel('Frequency (Hz)')
axT.set_xlabel('Frequency (Hz)')
axvn.set_ylabel('$v_n$ (cm/s)')
axx.set_ylabel('$x$ (nm)')
axP.set_ylabel('$\\Delta P$ (Pa)')
axT.set_ylabel('$\\Delta T$ ($\\mathrm{\\mu K}$)')
polish(figrest, 1, name='images//other_variables', extension='.png', grid=True,width_to_height = 0.5,tight_layout=True)        
polish(figgraph, 0.9, name='images//vs_simulation', extension='.png', grid=True,width_to_height = 1,tight_layout=True)        





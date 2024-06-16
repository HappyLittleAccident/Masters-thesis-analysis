# -*- coding: utf-8 -*-
"""
Created on Thu May 30 12:15:52 2024

@author: Marek
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp,quad,simpson
from time import time


gamma0   = 2*np.pi*1e-3
omega0   = 2*np.pi
gamma_l0 = 2*np.pi*1e-3
tau      = 1e3
fi = np.pi/2
L1s = [1,0]
L2s = [0,1]
L3s = [0,0]
colors = ['b','r']
fis = [0,np.pi/4,np.pi/2,3*np.pi/4,np.pi]
          
          
plt.close('all')




save_dict = {}
i=0
for L1,L2,L3,color in zip(L1s,L2s,L3s,colors):
    calculated_widths = []
    print(L1,L2,L3)

    for fi in fis:
        start = time()
        fig,ax = plt.subplots()
        def gamma_l(t):
            return(
                      gamma_l0 
                    + L1*gamma_l0*(np.sin(  omega0*t+fi))
                    + L2*gamma_l0*(np.sin(2*omega0*t+fi))
                    # + L3*gamma_l0*(2*np.random.rand()-1)
                    # + L3*gamma_l0*np.random.normal(scale=0.5)
                  )
        
        def function(t,y):
            vs  = y[0]
            dvs = y[1]
            
            rhs = np.array([
                   dvs,
                  -gamma0*dvs - gamma_l(t)*dvs - vs*omega0**2 - omega0*np.sin(omega0*t)             
                  ])
            return rhs
        
        
        t_span = [0,5000]
        
        y0 = np.array([1/(gamma0+gamma_l(0)),0])
        t_eval = np.linspace(0,t_span[-1],100*1000)
        result = solve_ivp(function, t_span,y0,t_eval=t_eval,vectorized=True,rtol = 1e-12)
        
        
        
        y = result['y']
        t = result['t']
        
    
        threshold = t_span[-1] - 60
        t_calc = t[t>threshold]
        X = simpson(y=(y[0][t>threshold])*np.exp(1j*omega0*t_calc),x=t_calc)/(t[-1]-t_calc[0])
        end = time()
        amp_vrms = 2*np.abs(X)
        print(f'Calculation took {end-start:.1f}s')
        
        start=time()
        amp = np.max(y[0][t>threshold])
        calculated_width =  1/amp_vrms - gamma0
        ax.plot(t,y[0],lw=0.5)
        # ax.axvline(t[t>threshold][0],c='k')
        # # plt.figure()
        
        # for t0 in t[:1000]:
        #     plt.plot(t0,gamma_l(t0),'o')
        

        
        print(f'fi= {fi}')
        print('\nvs from vrms:     ',amp_vrms,'\nvs from amplitude:',amp)
        print('\ngamma_L0 calculated:',1/amp_vrms - gamma0,'\ngamma_l0 set:\t\t',gamma_l0)
        print(f'Figure took {time()-start:.1f}s')
        print('-------------------------')
        calculated_widths.append(calculated_width)
        



    save_dict[f'L{i}'] = np.array([fis,calculated_widths])
    i+=1

save_dict['gamma_l'] = gamma_l0
np.save('parameters.npy',allow_pickle=True,arr=save_dict)

#%%
plt.close('all')
calculated_data = np.load('parameters_first_second.npy',allow_pickle=True).item()
L1 = calculated_data['L0']
L2 = calculated_data['L1']

fig_gamma,ax_gamma = plt.subplots()
ax_gamma.axhspan(gamma_l0*0.5,gamma_l0*1.5,color='silver',alpha = 0.5)

ax_gamma.plot(L1[0,:],L1[1,:],'o',c='b',label="First harmonic")
ax_gamma.plot(L2[0,:],L2[1,:],'s',c='r',label='Second harmonic')


ax_gamma.axhline(gamma_l0,c='k',label=r'Set damping $\alpha\gamma L_0$')
ax_gamma.axhline(gamma_l0*1.5,c='k',ls='--')
ax_gamma.axhline(gamma_l0*0.5,c='k',ls='--')
ax_gamma.text(1,gamma0*0.52,r'50% error')
ax_gamma.set_ylim(0,0.012)

ax_gamma.legend(loc = 'upper right')
ax_gamma.set_xlabel('Phase $\\varphi$ (rad)')
ax_gamma.set_ylabel('Calculated damping $\\alpha\kappa L_0$ (arb.)')




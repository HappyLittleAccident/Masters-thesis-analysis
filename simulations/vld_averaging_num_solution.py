# -*- coding: utf-8 -*-
"""
Created on Thu May 30 12:15:52 2024

@author: Marek
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp,quad,simpson


gamma0   = 2*np.pi*1e-3
omega0   = 2*np.pi
gamma_l0 = 2*np.pi*1e-3
tau      = 1e3
fi = np.pi/2

plt.close('all')



for fi in [0]:
    fig,ax = plt.subplots()
    def gamma_l(t):
        return(
                  gamma_l0 
                # + damping*gamma_l0*(np.sin(  omega0*t))
                + 0*gamma_l0*(np.cos(2*omega0*t+fi))
                # + damping*gamma_l0*np.random.rand()
              )
    
    def function(t,y):
        vs  = y[0]
        dvs = y[1]
        
        rhs = np.array([
               dvs,
              -gamma0*dvs - gamma_l0*dvs - vs*omega0**2 - omega0*np.sin(omega0*t+fi)             
              ])
        return rhs
    
    
    t_span = [0,1000]
    
    y0 = np.array([75.37,0])
    t_eval = np.linspace(0,t_span[-1],1*1000*1000)
    result = solve_ivp(function, t_span,y0,t_eval=t_eval,method='RK23',vectorized=True,rtol = 1e-10)
    
    
    
    y = result['y']
    t = result['t']
    

    threshold = 900
    t_calc = t[t>threshold]
    X = simpson(y=(y[0][t>threshold])*np.exp(1j*omega0*t_calc),x=t_calc)/(t[-1]-t_calc[0])
    
    amp_vrms = 2*np.abs(X)
    
    amp = np.max(y[0][t>threshold])
    ax.plot(t,y[0])
    ax.axvline(t[t>threshold][0])
    # plt.figure()
    
    # for t0 in t[:1000]:
    #     plt.plot(t0,gamma_l(t0),'o')
    
    print(f'fi= {fi}')
    print('\nvs from vrms:     ',amp_vrms,'\nvs from amplitude:',amp)
    print('\ngamma_L0 calculated:',1/amp_vrms - gamma0,'\ngamma_l0 set:\t\t',gamma_l0)
    print('-------------------------')






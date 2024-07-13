# -*- coding: utf-8 -*-
"""
Created on Wed Jul 10 11:17:29 2024

@author: Marek
"""

import numpy as np
import matplotlib.pyplot as plt
import helium_old as he
from scipy.interpolate import CubicSpline
from format_figure import use_latex,polish

plt.rcdefaults()
use_latex()
T_lambda = 2.1720
def lambda_line(T):
    x = (T-T_lambda)
    return (
            0.42800749 -95.0719*x -86.417*x**2 -103.341*x**3 -77.52175*x**4 
            -0.37827065*np.exp(42.2507*x)
            ) #atmospheres
def melting_curve(T):
    P=[]
    for t in T:
        if t<1.437:
            A= -24.8646 
            B= -3.12957 
            C=  2.19310 
        elif t>=1.437 and t<1.76:
            A= -25.1356
            B= -3.44479 
            C=  2.32876 
        else:
            A= 8.0016
            B= 0.98476 
            C= 0.33553
        P.append(10**(B+C*t) - A)
    return  P


plt.close('all')
T_crit = 5.1
T = np.linspace(0.001,2.5,500)
T_svp = np.linspace(0.001, T_crit,500)
T_lambda_span = np.linspace(1.763+0.001,T_lambda,500) 
svp = he.SVP(T_svp)/101325

fig,ax = plt.subplots(dpi=300)
ax.plot(T_svp,svp,'k')
ax.plot(T_crit,2.09,'ko',ms=5,fillstyle = 'full',markerfacecolor='white')
ax.plot(T_lambda_span,lambda_line(T_lambda_span),'k-')
ax.plot(T,melting_curve(T),'k-')
ax.set_xlim(0,5.5)
ax.set_ylim(0,40)
ax.annotate('Liquid He II',xy=(1,10))
ax.annotate('Liquid He I', xy=(3.5,20))
ax.annotate('Solid He',xy=(0.7,33))
ax.annotate('He gas',xy=(4.95,0.51))
ax.annotate('Critical point',xy=(T_crit-0.01,2.2),xytext=(4,5),arrowprops=dict(arrowstyle='->'))
ax.annotate(r'$\mathrm{\lambda}$-line',xy=(1.9,20),rotation=100+180)
ax.set_xlabel('Temperature (K)')
ax.set_ylabel('Pressure (atm.)')


polish(fig, 1, name='images//phase_diagram', extension='.png',width_to_height=0.85, grid=False,tight_layout=True)        


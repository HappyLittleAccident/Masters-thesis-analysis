# -*- coding: utf-8 -*-
"""
Created on Wed Jul 10 11:17:42 2024

@author: Marek
"""

import numpy as np
import matplotlib.pyplot as plt
import helium_old as he
from scipy.interpolate import CubicSpline
from format_figure import use_latex,polish


plt.close('all')
use_latex()
T = np.linspace(0.001,he.T_lambda-0.0001,500)


rho=he.density(T)
rho_s=he.superfluid_density(T)
rho_n=he.normal_fluid_density(T)


fig,ax = plt.subplots(dpi=300)
ax.plot(T,rho_s/rho,'b-')
ax.plot(T,rho_n/rho,'r-')
ax.hlines(1,he.T_lambda-0.0001,he.T_lambda+0.5,color='r')
ax.hlines(0,he.T_lambda-0.0001,he.T_lambda+0.5,color='b')
ax.axvline(he.T_lambda,ls='--',c='k')
ax.annotate(r'$\mathrm{\lambda}$-transition',xy=(he.T_lambda+0.05,0.4),rotation=90)
ax.annotate(r'$\rho_s / \rho$',xy=(1,0.9),c='b')
ax.annotate(r'$\rho_n / \rho$',xy=(1,0.08),c='r')
ax.set_xlabel('Temperature (K)')
ax.set_ylabel('Fraction')
polish(fig, 1, name='images//rhos_rho', extension='.png', grid=True,tight_layout=True)        

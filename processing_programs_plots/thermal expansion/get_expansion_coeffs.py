# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 15:56:28 2024

@author: Marek
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.interpolate import CubicSpline

plt.close('all')
d = np.loadtxt('thermal_expansion_svp.txt')
d = d[d[:,0].argsort()]
T = d[:,0]
alpha = d[:,1]*1e-4
plt.figure()
plt.plot(T,alpha,'o')

# def expansion(t,a,b,c):
#     return a*t**2 + b*t + c

# par,err = curve_fit(expansion,T,alpha)

expansion = CubicSpline(T,alpha)

print([temp for temp in d[:,0]])
print([temp for temp in d[:,1]])

x = np.linspace(T[0],T[-1],500)
plt.plot(x,expansion(x))
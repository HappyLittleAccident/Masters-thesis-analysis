# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 14:24:54 2024

@author: Marek
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import constants as con

h = con.h
kb = con.Boltzmann
T = 1.5
c2 = 6000

c4 = 145
D = 500e-9

nu = np.linspace(0.1,1e12,1000)



l = 50*1e-9
z1 = 15
z2 = 17
z3 = 0.034

alpha = 4*z1*z3/( (z1+z3)**2 * np.cos(2*np.pi*nu*l/c2)**2 + (z2 + z3*z1/z2)**2 * np.sin(2*np.pi*nu*l/c2)**2)
alpha_no_al = 4*z1*z3/( (z1+z3)**2 )


U = np.heaviside(nu-c4/D,0)*np.pi*h*nu**3 /(np.exp(h*nu/(kb*T))-1)/c2**2
plt.close("all")
plt.figure()
plt.plot(nu,U*alpha,label = 'Al')
plt.plot(nu,alpha_no_al*U,label = 'No Al')
plt.legend()
print(np.sum(alpha * U))
print(np.sum(alpha_no_al * U))
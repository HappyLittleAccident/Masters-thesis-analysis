# -*- coding: utf-8 -*-
"""
Created on Wed May 17 18:41:25 2023

@author: Admin
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
from glob import glob
import scipy as sc
import sys 
from scipy.optimize import curve_fit
import helium_old as he

from fitting_class import MicrowavePeak

"""
Cycle trough all temperatures
"""
temperatures = ['T1325',
 'T1350',
 'T1350_2',
 'T1350_4',
 'T1375',
 'T1400',
 'T1425',
 'T1450',
 'T1450_2',
 'T1475',
 'T1500',
 'T1550',
 'T1650',
 'T1700',
 'T1750',
 'T1850',
 'T1900',
 'T1950']

"""
For each letter, cycle trough all temperatures
"""
letters = ['A','B','C','D']

resonator_A = {'kp':1.2400632746323976e7,'C0':338.63213961685864e-12,'l':1050e-6,'w':1000e-6}
resonator_B = {'kp':1.4546169638969617e7,'C0':335.6978986554392e-12,'l':750e-6,'w':500e-6}
resonator_C = {'kp':1.1909096177342793e7,'C0':312.8732100355237e-12,'l':1000e-6,'w':950e-6}
resonator_D = {'kp':1.6174502699593267e7,'C0':466.4158697787019e-12,'l':950e-6,'w':500e-6}

resonators = {'A':resonator_A,'B':resonator_B,'C':resonator_C,'D':resonator_D}

def lin(x,a,b):
    y = a*x + b
    return y

def cub(x,a,b,c,d):
    y = a*x**3 + b*x**2 + c*x + d
    return y

def compres(T):
    """
    He II compresibility
    """
    a,b,c,d = (3.65456274e-10, -1.19112443e-09,  1.72705353e-09, 1.12838569e-08)
    chi = a*T**3 + b*T**2 + c*T + d
    return chi
  

def Pgrad (bias,drive,T,R,h,l,kp,C0):
    """
    Return pressure gradient from voltage on chip and bias
    """
    A = np.pi*R**2
    VB = A*h
    chi = compres(T)
    Pgrad = (2*A*C0*bias)/(chi*VB*kp + 2*A**2)*drive/(h*l)
    return Pgrad

def Velocity (bias,I,T,R,h,rho,w,kp,C0):
    """
    Return He II velocity from measured current
    """
    A = np.pi*R**2
    VB = A*h
    a = h*w
    chi = compres(T)
    rhos = he.superfluid_density(T)
    er = 1.0565
    SIG = chi*h*kp/(2*A)
    #print(SIG)
    g = (2*a*rhos)/(VB*rho)*(1 - 2*(er - 1)/er*SIG)/(1 + 2*SIG)
    v = I/(g*bias*C0)
    return v

plt.close('all')

for letter in letters:
    """
    Save params for each resonator
    """
    
    for temperature in temperatures:
        height = 500 #nm
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        

# -*- coding: utf-8 -*-
"""
Created on Tue May 16 19:06:47 2023

@author: Admin
"""

import numpy as np
import matplotlib.pyplot as plt
from glob import glob
from scipy.optimize import curve_fit

class circuitlosses:
    
    def __init__ (self, f0):
        self.f0 = f0
        self.par = [ 1.64991875e-14,  1.11321488e-09, -1.08292915e-06,  6.27560887e-01]
        
    def cubic_func (self,f,par):
        a,b,c,d = par
        #print(a,b,c,d)
        y = a*f**3 + b*f**2 + c*f + d
        return y
    
    def get_loss(self):
        self.loss = self.cubic_func(self.f0, self.par)
        return self.loss
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 23 14:08:24 2022

@author: ev
"""

import numpy as np
import lmfit as lm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from scipy.signal import savgol_filter

def microwave_peak(f, f0, A,gamma, a1, b1, a2, b2,phi):
    z = f0*A*gamma/(f0**2 - f**2 + 1j*f*gamma)*np.exp(1j*phi) + ((a1*f + b1)*1j + (a2*f + b2))
    return z


def microwave_peak_nophi(f, f0, A,gamma, a1, b1, a2, b2):
    z = f*A/(f0**2 - f**2 + 1j*f*gamma) + (a1*(f-2000) + b1)*1j + (a2*(f-2000) + b2)
    return z


class MicrowavePeak:
    def __init__(self, f, x, y, init_p):
        self.f = f#/1e9
        self.x = x
        self.y = y
        self.init_p = init_p
        
        self.z = self.x + 1j*self.y
        self.model = lm.Model(microwave_peak, independent_vars='f')

        self.ix1 = np.argmax(self.x)
        self.ix2 = np.argmin(self.x)
        self.gamma_est = abs(self.f[self.ix1] - self.f[self.ix2])/2
        #self.p0 = self.model.make_params(A=np.max(np.sqrt(self.x**2 + self.y**2)), f0=np.mean(self.f), gamma=self.gamma_est, a1=1, b1=1, a2=1, b2=1)     
        #self.p0 = self.model.make_params(A=self.init_p[1], f0=self.init_p[0], gamma=self.init_p[2], a1=self.init_p[3], b1=self.init_p[4], a2=self.init_p[5], b2=self.init_p[6], phi = self.init_p[7])
        self.p0 = self.model.make_params(A=self.init_p[1], f0=self.init_p[0], gamma=self.init_p[2], a1=self.init_p[3], b1=self.init_p[4], a2=self.init_p[5], b2=self.init_p[6],phi = self.init_p[7])
    def fit(self, plot=False):
        self.p0['a1'].set(vary=False)
        self.p0['a2'].set(vary=True)
        self.out= self.model.fit(self.z, f=self.f, params=self.p0,weights = 1/np.sqrt(self.x**2 + self.y**2))
            
        return self.out

# =============================================================================
# if __name__ == '__main__':
#     import os.path as path
#     from glob import glob
#     
#     letter = 'A'
#     height = 500
#     
#     files = glob(rf'D:\2023_4_Helmholtz_resonators\Helmholtz_res_drive_dep\T1950\fsweeps\resonator_{height:}{letter:}'+ '\\' + '*.npy')
#     par = []
#     a1s = []
#     b1s = []
#     a2s = []
#     b2s = []
#     drs = []
#     fig,ax = plt.subplots(2,1)
#     for num, file in enumerate(files[10:20]):
#         #print(num)
#         d = np.load(file, allow_pickle=True).item()
#         
#         f = np.array(d['freq (Hz)'])
#         x = np.array(d['x (A)'])
#         y = np.array(d['y (A)'])
#         r = np.sqrt(x**2 + y**2)
#         amp = d['App (V)']
#         
#         if num == 0:
#             A= np.max(np.sqrt(x**2 + y**2))
#             ix1 = np.argmax(x)
#             ix2 = np.argmin(x)
#             gamma_est = abs(f[ix1] - f[ix2])/2
#             ixf = np.argmax(np.sqrt(x**2 + y**2))
#             f0_est = f[ixf]
#             A= np.max(np.sqrt(x**2 + y**2))*gamma_est
#             #b1est = np.mean(y[:10] + y[-10:])
#             #b2est = np.mean(x[:10] + x[-10:])
#             #print(f0_est,A,gamma_est)
#             init_p = [f0_est, A, gamma_est,0,0,0,0]
#         else:
#             #print(par)
#             init_p = [par['f0'],par['A'],par['gamma'],par['a1'],par['b1'],par['a2'],par['b2']]
#          
#         #A= np.max(np.sqrt(x**2 + y**2))             
#         peak = MicrowavePeak(f, x, y, init_p)
#         result = peak.fit(plot=True)
#         par = result.best_values
#         res_fit = result.best_fit
#         ax[0].plot(f, x, 'o',ms = 3)
#         ax[1].plot(f, y, 'o', ms = 3)
#         ax[0].plot(f, np.real(res_fit), 'k--')
#         ax[1].plot(f, np.imag(res_fit), 'k--')
#         
#         a1s.append(par['a1'])
#         b1s.append(par['b1'])
#         a2s.append(par['a2'])
#         b2s.append(par['b2'])
#         drs.append(amp)
#         
# fig2,ax2 = plt.subplots(2,2)
# ax2[0,0].plot(drs,a1s,'o')
# ax2[1,0].plot(drs,b1s,'o')
# ax2[0,1].plot(drs,a2s,'o')
# ax2[1,1].plot(drs,b2s,'o')
# =============================================================================
        
        
    

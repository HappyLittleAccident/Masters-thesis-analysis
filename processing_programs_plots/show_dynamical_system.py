# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 13:16:52 2024

@author: Marek
"""

import numpy as np
import matplotlib.pyplot as plt
from format_figure import use_latex,polish


plt.close('all')



gs = np.linspace(0,4,1000)
Ls = np.linspace(0,4,1000)

GS,LS = np.meshgrid(gs,Ls)

def Ldot(L,d,g):
    return -L -d*(L/(L**2 + 1))**2 + g

colors = ['tab:blue','tab:orange']

figf,axf = plt.subplots(dpi = 200)
ds = [2,6]
for i,d in enumerate(ds):
    condition = np.logical_or(Ldot(LS-0.1,d,GS) < 0,Ldot(LS+0.1,d,GS) > 0)
    
    result = np.ma.MaskedArray(Ldot(LS,d,GS),condition)
    result_inverse=  np.ma.MaskedArray(Ldot(LS,d,GS),np.logical_not(condition))
    
    axf.contour(GS,LS,result,levels=[0],colors=[colors[i]],linewidths=[1])
    axf.contour(GS,LS,result_inverse,levels=[0],colors=[colors[i]],linestyles=['--'],linewidths=[1])

axf.axvline(2.35,lw=0.5,c='k')
axf.axvline(2.6,lw=0.5,c='k')

fig_ldot,ax_ldot = plt.subplots(dpi=200)






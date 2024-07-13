# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 11:15:41 2023

@author: Marek
"""

import numpy as np
import matplotlib.pyplot as plt
from format_figure import use_latex,polish


plt.close('all')
x = np.linspace(-1,1,50)
pos = np.array([-0.7,0.7])

# fig_F,ax_F = plt.subplots()

# ax_F.plot(x,x**2,'k-')
# ax_F.plot(pos,pos**2,'ro')
# ax_F.set_xticks([0],[0])
# ax_F.set_yticks([0],['$V_0$'])
# ax_F.set_xlabel('oscillator position x')
# ax_F.set_ylabel('potential energy V(x)')
# ax_F.xaxis.label.set_fontsize(15)
# ax_F.yaxis.label.set_fontsize(15)
# fig_F.tight_layout()
# fig_F.savefig('oscillator.png')

use_latex()
x0 = np.linspace(0,6,500)
figE, axE = plt.subplots(dpi=300)
axE.plot(x0,4.5-x0-15*(x0/(x0**2 +1))**2,'k-',label='$g = 4$')
axE.plot(x0,1-x0-15*(x0/(x0**2 +1))**2,'k--',label='$g = 1$')
axE.plot(0.25,0,'ko',markerfacecolor='black',ms=10,label='Stable fixed point')
axE.set_xticks([0,0.85,1.75,3.42,0.25],[0,'$x_1$','$x_2$','$x_3$',"$x_1'$"])
axE.set_yticks([0],[0])
axE.set_xlabel('Variable $x$')
axE.set_ylabel('Function $f(x)$')
axE.axhline(0,c='k')
axE.plot(0.85,0,'ko',markerfacecolor='black',ms=10)
axE.plot(1.75,0,'ko',markerfacecolor='white',ms=10,label='Unstable fixed point')
axE.plot(3.42,0,'ko',markerfacecolor='black',ms=10)
axE.set_xlim(0,6)
axE.annotate('g',(5,-3))
axE.annotate('',(4.9,-2),(4.9,-3.8),arrowprops=dict(arrowstyle='->'))
axE.legend()
polish(figE, 1, name='images//fixed_points', extension='.png', grid=True,tight_layout=True)        


# x0 = np.linspace(0.0001,0.01,10)
# figZ, axZ = plt.subplots()
# axZ.plot(x0,x0**2,'ko')
# axZ.set_xticks([0,*x0[4:6]],[0,'',''])
# axZ.set_yticks([0,*x0[4:6]**2],[0,'',''])
# axZ.set_xlabel('oscillator amplitude $x_0$')
# axZ.set_ylabel('oscillator energy E')
# axZ.xaxis.label.set_fontsize(15)
# axZ.yaxis.label.set_fontsize(15)
# figZ.tight_layout()
# figZ.savefig('energy_levels_zoomed.png')

# x0 = np.linspace(0,1,1000)
# figC, axC = plt.subplots()
# axC.plot(x0,np.sin(x0*2*np.pi),'b-')
# axC.set_xlabel('oscillator amplitude $x_0$')
# axC.set_ylabel('oscillator energy E')
# axC.set_xticks([],[])
# axC.set_yticks([],[])
# axC.xaxis.label.set_fontsize(15)
# axC.yaxis.label.set_fontsize(15)
# figC.tight_layout()
# figC.savefig('cavity.png')




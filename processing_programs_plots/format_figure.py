# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 16:55:10 2023

@author: Marek
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import time

n=1
thesis = 412.56496 #pt

golden_ratio = (5**.5 - 1) / 2
point_to_inch = 1/72.27


def use_latex(font_size = 12):
    tex_fonts = {
        # Use LaTeX to write all text    
        "text.usetex": True,
        "font.family": "serif",
        "text.latex.preamble":r'\usepackage{lmodern}',
        "font.size": font_size,
    }
    plt.rcParams.update(tex_fonts)
    return 0

def polish(fig,fraction,name='',pagewidth=thesis,width_to_height = golden_ratio,extension = '.pdf',grid=True):
    width_inch = point_to_inch*pagewidth
    width = fraction*width_inch
    height = width*width_to_height
    
    gs = fig.axes[0].get_gridspec()
    gs.nrows  # return  1
    gs.ncols  # returns 4
    

    if name == '':
        global n        
        name = str(n)
        n+=1
    

    fig.set_size_inches(width, height*gs.nrows/gs.ncols, forward=True)

    if grid:
        for ax in fig.axes:
            ax.grid('both')
    fig.savefig(name+extension)
    
    

    return 0
        



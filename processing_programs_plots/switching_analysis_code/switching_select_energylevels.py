# -*- coding: utf-8 -*-
"""
Created on Fri Apr 21 20:36:56 2023

@author: Marek
"""

import numpy as np
import matplotlib.pyplot as plt
import os.path as path
from glob import glob
from scipy import interpolate
import os

plt.close('all')
save = False
basefolder = r'D:\OneDrive_copy\OneDrive - Univerzita Karlova\DATA\2023_4_Helmholtz_resonators_recal_2'


class get_drive():
    """
    Selects energy levels for single temperature
    """
    def __init__(self,figure,axis,files_switching,temp):
        self.fig = figure
        self.ax = axis                      
        self.files = files_switching   
        
        self.v_points = []
        self.p_points = []
    
        self.fig.canvas.mpl_connect('button_press_event', self.mouse_click)
        
        self.rs = []
        self.amps = []

        amp_last = 0
        self.vs=[]
        self.ps=[]

        for f in self.files:
            data = np.load(f, allow_pickle=True).item()
            
            r = np.array(data['Velocity (m/s)'])
            P = np.array(data['Pressure (Pa/m)'])
            amp = np.array(data['App_gen (V)'])            
            
            self.rs.append(r)
            self.amps.append(amp)     
#%%            
            if amp_last == amp:
                self.vs[-1]=np.append(self.vs[-1],r)              
            else:
                self.ps.append(P)
                self.vs.append(r)
            amp_last = amp
                
            
        cmap = plt.get_cmap('inferno')
                
        vs_hists = []
        max_vs = np.max(np.concatenate(self.vs))
        min_vs = np.min(np.concatenate(self.vs))
        for vs in self.vs:
            hist,edges = np.histogram(vs,bins=np.linspace(min_vs,max_vs,400))
            vs_values=0.5*(edges[1:]+edges[:-1])
            vs_hists.append(hist)
            self.vs_edges=edges
        
        Ps,Vs = np.meshgrid(self.ps,vs_values)
        self.ax.pcolormesh(Ps,Vs,np.log(np.array(vs_hists).T),cmap = cmap,shading = 'nearest')
        self.ax.set_xlabel('Pressure (Pa/m)')
        self.ax.set_ylabel('Velocity (m/s)')
        self.ax.set_title(f'Temperature {temp}')
        self.fig.savefig(f'3Dhist_{temp}K.png')
        
        keyboard = False            
        while not keyboard:
            keyboard = plt.waitforbuttonpress()  
        
    def mouse_click(self,event):

        if event.button == 2:
            if event.ydata == None:
                print('Click the graph!')
            else:
                self.v_points.append(event.ydata)
                self.p_points.append(event.xdata)
                self.ax.plot(event.xdata,event.ydata,'r+')
                self.fig.canvas.draw()
                self.fig.canvas.flush_events()


        if event.button == 3:
            def interpolated(x):
                x_points = self.p_points
                y_points = self.v_points
            
                tck = interpolate.splrep(x_points, y_points,k=1)
                return interpolate.splev(x, tck)
            
            p_linspace =np.linspace(self.p_points[0],self.p_points[-1],100)
            self.ax.plot(p_linspace,interpolated(p_linspace))
            self.fig.canvas.draw()
            self.fig.canvas.flush_events()

switching_sweeps_names = {'1.35':'T1350','1.45':'T1450','1.65':'T1650','1.85':'T1850'}
temps = ['1.35','1.45','1.65','1.85']
letters = ['A','B','C','D']

for letter in letters:
    print(f'Resonator: {letter}')
    if letter in ['B']:
        for temp in temps[:]:
            params = {}
            print(f'Temperature: {temp}')
            date = '*'
            distance = 500
               
            folder = switching_sweeps_names[temp]
            basedir = fr'{basefolder}\Helmholtz_buffer_measurements\{folder}\buffer'
        
            resdir = path.join(basedir,r'resonator_{}{}'.format(distance,letter))
            files = glob(path.join(resdir, r'*.npy'))            
            files.sort(key=lambda f: np.load(f, allow_pickle=True).item()['App_gen (V)'])            
        
#%%               
            fig,ax = plt.subplots()
            obj = get_drive(fig,ax,files,temp)
            params['v points'] = obj.v_points
            params['p points'] = obj.p_points
            
            if save:
                folder = 'lifetimes_energylevels'
                os.makedirs(folder,exist_ok=True)
                np.save(os.path.join(folder,f'params_{letter}_T{temp}_new.npy'),params)












# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 09:10:47 2023

@author: Marek
"""


import numpy as np
import matplotlib.pyplot as plt
import os.path as path
from glob import glob
from matplotlib.lines import Line2D
import scipy.signal as signal
from scipy.optimize import curve_fit
from scipy import fft
from scipy import interpolate
import os

plt.close('all')
save = False
FFT = True
show = False


def gaussian(x,x0,w):
    return np.exp(-0.5*((x-x0)/w)**2)/(w*np.sqrt(2*np.pi))

def two_gaussian(x,x0,w,x01,w1,a):
    return gaussian(x,x0,w)*(a**2)/(a**2 + 1) + gaussian(x,x01,w1)/(a**2 + 1)


class get_drive():
    def __init__(self,figure,axis,files_switching,temp):
        self.fig = figure
        self.ax = axis
        
        self.figax = []
        self.figopen = False
        # self.fig_hist,self.ax_hist= plt.subplots()
        
        self.n_click = 0
        red_line = Line2D([0],[0],c='r',label = 'Start point')
        green_line = Line2D([0],[0],c='g',label = 'End point')
        blue_line = Line2D([0],[0],c='b',label = 'Shown measurement')
        # self.ax.legend(handles = [red_line,green_line,blue_line])
        if show:
            fig_fft,ax_fft = plt.subplots()
            self.fig_meas,self.ax_meas = plt.subplots()
        
        self.files = files_switching
        
        self.x1 = None
        self.x2 = None
        
        self.y1_list = []
        self.y2_list = []
        self.y3_list = []
        self.dy1_list = []
        self.dy2_list = []
        
        self.v_points = []
        self.p_points = []
        
        self.dy1_p = None
        self.dy1_m = None
        self.dy2_p = None
        self.dy2_m = None
        
        self.y1 = None
        self.y2 = None
        self.y3 = None
        self.dy1 = None
        self.dy2 = None
        
        self.pick_y1 = None
        self.pick_y2 = None
        self.pick_y3 = None
        
        self.red_line = None
        self.green_line = None
        
        self.data = None
        
        self.fig.canvas.mpl_connect('pick_event',self.pick_measurement)        
        self.fig.canvas.mpl_connect('button_press_event', self.pick_start)
        self.fig.canvas.mpl_connect('key_press_event', self.key_press)
        
        self.xs = []
        self.ys = []
        self.rs = []
        self.amps = []
        self.ts = []

        r_last = 0
        first = True
        p02_last=None
        amp_last = 0
        self.vs=[]
        self.ps=[]
        ffts = []
        for f in self.files:
            data = np.load(f, allow_pickle=True).item()
            
            upP = data['upPressure (Pa/m)']
            downP = data['downPressure (Pa/m)']
            jumpP = data['jumpPressure (Pa/m)']
            upV=data['upVelocity (m/s)']
            downV=data['downVelocity (m/s)']
            jumpV=data[ 'jumpVelocity (m/s)']
            
            # self.vs.append(upV)
            self.vs.append(downV)
            self.vs.append(jumpV)
            
            
            # self.ps.append(upP)
            self.ps.append(downP)
            self.ps.append(jumpP)
            

                
            
        cmap = plt.get_cmap('inferno')
        self.vs = np.concatenate(self.vs)
        self.ps = np.concatenate(self.ps)
        
        
        
        # hist,p_bins,v_bins = np.histogram2d(self.ps,self.vs,bins = [50,200])
        # vs_values=0.5*(v_bins[1:]+v_bins[:-1])
        # ps_values = 0.5*(p_bins[1:]+p_bins[:-1])
        
        
        # Ps,Vs = np.meshgrid(ps_values,vs_values)
        
        self.ax.hist2d(self.ps,self.vs,bins = 500,norm='log',cmap = cmap)
        
        x_min, x_max = self.ax.get_xlim()
        x_min -= 0.1*(x_max - x_min)
        x_max += 0.1*(x_max - x_min)
        self.ax.set_xlim(x_min,x_max)
        
        y_min, y_max = self.ax.get_ylim()
        y_min -= 0.1*(y_max - y_min)
        y_max += 0.1*(y_max - y_min)
        self.ax.set_ylim(y_min,y_max)
        
        self.ax.set_xlabel('Pressure (Pa/m)')
        self.ax.set_ylabel('Velocity (m/s)')
        self.ax.set_title(f'Temperature {temp}')
        self.fig.savefig(f'3Dhist_ampsweeps_{temp}K.png')

        while not plt.waitforbuttonpress():
            pass  



     

 
        
    def pick_start(self,event):
        
        previous= False
        if event.button == 2:
            if previous:
                self.p_points = []
                self.v_points = []
            self.v_points.append(event.ydata)
            self.p_points.append(event.xdata)
            self.ax.plot(event.xdata,event.ydata,'r+')
            self.fig.canvas.draw()
            self.fig.canvas.flush_events()
            # self.n_click += 1
            # if self.n_click%2 == 1:                
            #     # if not self.red_line == None:
            #     #     self.red_line.remove()
            #     # self.x1 = event.xdata
            #     # self.red_line = self.ax.axvline(self.x1,c='r')
            # elif self.n_click%2 == 0:
            #     if not self.green_line == None:
            #         self.green_line.remove()
            #     self.x2 = event.xdata
            #     self.green_line = self.ax.axvline(self.x2,c='g')
            
            # self.fig.canvas.draw()
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
            previous = True
            # if self.figopen:
            #     for fig_m in self.figax:
            #         plt.close(fig_m[0])
            #     self.figax = []
                
    def pick_level(self,event):
        if event.button == 1:
            self.y1 = event.ydata

            if not self.pick_y1 == None:
                self.pick_y1.remove()
            self.pick_y1 = self.ax_meas.axhline(self.y1,c='r')
        if event.button == 3:
            self.y2 = event.ydata

            if not self.pick_y2 == None:
                self.pick_y2.remove()
            self.pick_y2 = self.ax_meas.axhline(self.y2,c='b')
        if event.button == 2:
            self.y3 = event.ydata
            if not self.pick_y3 == None:
                self.pick_y3.remove()
            self.pick_y3 = self.ax_meas.axhline(self.y3,c='g')
            
            # if self.y2 == None:
            #     self.dy1 = np.abs(event.ydata-self.y1)
            # elif self.dy1==None:
            #     self.dy1 = np.abs(event.ydata-self.y1)
            # elif event.ydata > self.y1-self.dy1:
            #     self.dy1 = np.abs(event.ydata-self.y1)
            # else:
            #     self.dy2 = np.abs(event.ydata - self.y2)

        # arr_dy = [self.dy1_m,self.dy1_p,self.dy2_m,self.dy2_p]
        # for margin_line in arr_dy:
        #     if not margin_line == None:
        #         margin_line.remove()        

        # if not self.dy1 == None:
        #     self.dy1_m = self.ax_meas.axhline(-self.dy1 + self.y1,c='r',ls='--')
        #     self.dy1_p = self.ax_meas.axhline(self.dy1 + self.y1,c='r',ls='--')
        #     if not self.y2 == None and not self.dy2 == None:
        #         self.dy2_m = self.ax_meas.axhline(-self.dy2 + self.y2,c='b',ls='--')
        #         self.dy2_p = self.ax_meas.axhline(self.dy2 + self.y2,c='b',ls='--')
        
        if not self.y1 == None:
            self.ax_hist.axvline(self.y1,c='r')
        if not self.y2 == None:
            self.ax_hist.axvline(self.y2,c='b')
        if not self.y3 == None:
            self.ax_hist.axvline(self.y3,c='g')
        self.fig_meas.canvas.draw()
        self.fig_hist.canvas.draw()   
        
    
    def pick_measurement(self,event):
        if self.figopen:
            for fig_m in self.figax:
                plt.close(fig_m[0]) 
            self.figax = []
            self.marked.remove()
        self.artist = event.artist
        D_amp = self.artist.get_xdata()
        for i,amp in enumerate(self.amps):
            if D_amp == amp:
                fig_temp,ax_temp = plt.subplots(3,1,sharex = True)
                self.figax.append([fig_temp,ax_temp])

            
                self.figax[-1][1][0].plot(self.ts[i],self.xs[i],'.-',c= 'tab:blue',ms = 1,lw=1)
                self.figax[-1][1][1].plot(self.ts[i],self.ys[i],'.-',c = 'tab:orange',ms = 1,lw=1)
                self.figax[-1][1][2].plot(self.ts[i],np.abs(self.xs[i]+1j*self.ys[i]),'.-',c = 'tab:green',ms = 1,lw=1)
                self.figax[-1][1][2].set_title(f'Amp={D_amp}')
                self.figopen = True
        self.marked = self.ax.axvline(D_amp,c='b')
        self.fig.canvas.draw()
    def key_press(self,event):
        if self.figopen:
            for fig_m in self.figax:
                plt.close(fig_m[0])
            self.figax = []
    def key_press_select(self,event):
        self.data[0].remove()
        
        
    
        
    
def function(t,a,b,c,d):
    return d*t**3 + b*t**2 + c*t + d

switching_sweeps_names = {'1.35':'T1350_3','1.45':'T1450_3','1.65':'T1650_3','1.85':'T1850_2'}
switching_ampsweeps_names = {'1.35':'T1350','1.45':'T1450','1.65':'T1650','1.85':'T1850'}

temps = ['1.35','1.45','1.65','1.85']
letters = ['A','B','C','D']


ampsweeps_folders=glob(r'D:\OneDrive_copy\OneDrive - Univerzita Karlova\DATA\2023_4_Helmholtz_resonators_recal\Helmholtz_res_drive_dep\*')
ampsweeps_temps = [x.split('\\')[-1] for x in ampsweeps_folders if not '_' in x.split('\\')[-1]]


for letter in letters:
    print(f'Resonator: {letter}')
    if letter in ['B']:
        for temp in ampsweeps_temps:
            params = {}
            print(f'Temperature: {temp}')
            # date = '*'
            distance = 500
               
            # folder = switching_sweeps_names[temp]
            # folder_ampsweeps = switching_ampsweeps_names[temp]
            
            # basedir = r'D:\OneDrive_copy\OneDrive - Univerzita Karlova\DATA\2023_4_Helmholtz_resonators_recal\Helmholtz_res_drive_dep\{}\buffer'.format(folder)
            
            # ampsweeps_base = r'D:\OneDrive_copy\OneDrive - Univerzita Karlova\DATA\2023_4_Helmholtz_resonators_recal\Helmholtz_res_drive_dep\{}\ampsweeps_fast'.format(folder_ampsweeps)
            
            # resdir = path.join(basedir,r'resonator_{}{}'.format(distance,letter))
            # files = glob(path.join(resdir, r'*.npy'))
            
            # ampsweeps = path.join(ampsweeps_base,r'resonator_{}{}_updowndown'.format(distance,letter))
            # files_ampsweeps = glob(path.join(ampsweeps, r'*.npy'))
            
            # # files = files[:-53]
            folder_ampsweeps = fr'D:\OneDrive_copy\OneDrive - Univerzita Karlova\DATA\2023_4_Helmholtz_resonators_recal\Helmholtz_res_drive_dep\{temp}'
            files_ampsweeps = glob(folder_ampsweeps + fr'\ampsweeps_fast\resonator_{distance}{letter}_updowndown\*.npy')
            
            fig_ampsweeps, ax_ampsweeps = plt.subplots()
            
            cmap = plt.get_cmap('viridis')
            cmap_hist = plt.get_cmap('Reds')
        
            #%% plot ampsweeps
            for file in files_ampsweeps[:1]:#files[:20] and files[39:]:
                # print(file)
                d = np.load(file, allow_pickle=True).item()
                        
                upP = d['upPressure (Pa/m)']
                downP = d['downPressure (Pa/m)']
                jumpP = d['jumpPressure (Pa/m)']
                upV=d['upVelocity (m/s)']
                downV=d['downVelocity (m/s)']
                jumpV=d[ 'jumpVelocity (m/s)']
                
                ax_ampsweeps.plot(upP, upV, '-o',ms=0.5, lw=0.5, color='tab:blue',alpha=1)
                ax_ampsweeps.plot(downP,downV, '-o',ms=0.5, lw=0.5, color='tab:orange',alpha=1)
                ax_ampsweeps.plot(jumpP, jumpV, '-o',ms=0.5, lw=0.5, color='tab:green',alpha=1)
                ax_ampsweeps.set_xlabel('Pressure (Pa/m)')
                ax_ampsweeps.set_ylabel('Velocity (m/s)')
               
            fig,ax = plt.subplots()
            obj = get_drive(fig,ax,files_ampsweeps,temp)
            
            # keyboard = False
            # while not keyboard:
            #     keyboard = plt.waitforbuttonpress() 
    
            params['x1'] = obj.x1 #in normalised drive amplitude
            params['x2'] = obj.x2 
            params['y1'] = obj.y1_list
            params['y2'] = obj.y2_list
            params['dy1'] = obj.dy1_list
            params['dy2'] = obj.dy2_list
            params['v points'] = obj.v_points
            params['p points'] = obj.p_points
            # plt.close(fig_ampsweeps)
            if save:
                os.makedirs("lifetimes_params_ampsweeps_processed",exist_ok=True)
                np.save(f'lifetimes_params_ampsweeps_processed/params_{letter}_T{temp}_2.npy',params)












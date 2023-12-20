# -*- coding: utf-8 -*-
"""
Created on Fri Apr 21 20:36:56 2023

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

plt.close('all')
save = True
FFT = False
cutoff = (25,np.inf)

def gaussian(x,x0,w):
    return np.exp(-0.5*((x-x0)/w)**2)/(w*np.sqrt(2*np.pi))

def two_gaussian(x,x0,w,x01,w1,a):
    return gaussian(x,x0,w)*(a**2)/(a**2 + 1) + gaussian(x,x01,w1)/(a**2 + 1)


class get_drive():
    def __init__(self,figure,axis,files_switching):
        self.fig = figure
        self.ax = axis
        self.figax = []
        self.figopen = False
        self.fig_hist,self.ax_hist= plt.subplots()
        
        self.n_click = 0
        red_line = Line2D([0],[0],c='r',label = 'Start point')
        green_line = Line2D([0],[0],c='g',label = 'End point')
        blue_line = Line2D([0],[0],c='b',label = 'Shown measurement')
        self.ax.legend(handles = [red_line,green_line,blue_line])
        
        self.fig_meas,self.ax_meas = plt.subplots()
        
        self.files = files_switching
        
        self.x1 = None
        self.x2 = None
        
        self.y1_list = []
        self.y2_list = []
        self.y3_list = []
        self.dy1_list = []
        self.dy2_list = []
        
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
        for f in self.files:
            data = np.load(f, allow_pickle=True).item()
            x = np.array(data['x (A)'])
            y = np.array(data['y (A)'])
            t = data['time']
            amp = np.array(data['App_gen (V)'])/b
            self.sample_rate = data['sample rate (Hz)']
            r= np.abs(x + 1j*y)
            
            self.xs.append(x)
            self.ys.append(y)
            self.rs.append(r)
            self.ts.append(t)
            self.amps.append(amp)
            
            
            
            if r_last ==0:
                r_last = np.mean(r)
            shift = np.mean(r)-r_last
            if not self.y1 == None:
                self.y1 += shift
            if not self.y2 == None:
                self.y2 += shift
            if not self.y3 == None:    
                self.y3 += shift
            if not self.pick_y1 == None:
                self.pick_y1.remove()
                self.pick_y1 = self.ax_meas.axhline(self.y1,c='r')
            if not self.pick_y2 == None:
                self.pick_y2.remove()
                self.pick_y2 = self.ax_meas.axhline(self.y2,c='b')
            if not self.pick_y3 == None:
                self.pick_y3.remove()
                self.pick_y3 = self.ax_meas.axhline(self.y3,c='g')
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
        

            self.ax.plot(amp,np.mean(r),'ro',picker = 5,alpha = 0.5)
            self.fig.canvas.draw()
            

            meas_red_line = Line2D([0],[0], c='r',label = 'y1')
            meas_blue_line = Line2D([0],[0], c='b',label = 'y2')
            
            if FFT:
                sos = signal.butter(6, 2*30/self.sample_rate, output='sos')

                r = signal.sosfiltfilt(sos, r)
                # transformed = fft.fft(x+1j*y)
                # freq = fft.fftfreq(len(r),1/self.sample_rate)
                
                # for fq1,fq2 in [cutoff]:
                #     index = np.logical_and(np.abs(freq)>fq1,np.abs(freq)<fq2)
                #     transformed[index]=np.zeros_like(transformed)[index]
    
                # inverse = np.abs(fft.ifft(transformed))
                # r=inverse
            # w=250
            # r = np.convolve(r,np.ones(w),mode='same')[w:-w]/w
            # t = np.copy(t[w:-w])
            # b1 = [1.0 / n] * n
            # a = 1
            # zi = signal.lfilter_zi(b1, a)
            
            # r,_ = signal.lfilter(b1,a,r,zi = r[0]*zi)
            # r = signal.savgol_filter(r,200,4)
            
            self.data = self.ax_meas.plot(t,r,'.-',c = 'k',ms = 1,lw=0.1)
            
            self.ax_meas.legend(handles = [meas_red_line,meas_blue_line])
            self.fig_meas.canvas.draw()
            
            self.fig_meas.canvas.mpl_connect('button_press_event',self.pick_level)
            self.fig_meas.canvas.mpl_connect('key_press_event',self.key_press_select)
            self.fig_meas.canvas.draw()
            


            
            hist,edges = np.histogram(r,70,density=True)
            self.ax_hist.clear()
            qty_values=0.5*(edges[1:]+edges[:-1])
            self.ax_hist.plot(qty_values,hist,'.-k',ms=1)
            
            if not self.y1 == None:
                self.ax_hist.axvline(self.y1,c='r')
            if not self.y2 == None:
                self.ax_hist.axvline(self.y2,c='b')
            if not self.y3 == None:
                self.ax_hist.axvline(self.y3,c='g')
            skewness= np.mean((r-np.mean(r))**3)    
            """
            Fit gaussian on data
            """
            # try:
            #     p0s = [np.mean(r),np.std(r)]
            #     par,sig = curve_fit(gaussian,qty_values,hist,p0=p0s)
            #     chi = np.sum((gaussian(qty_values,*par)/(hist+1) -1)**2)/len(hist)
            #     # self.ax_hist.plot(qty_values,gaussian(qty_values,*par),'r-')

            # except:
            #     print('error gauss')
            
            
            # two_gauss = False
            try:

                
                if True or first or np.array_equal(p02_last,None):
                    p0s = [r[len(r)//4],(qty_values[-1]-qty_values[0])/8]
                    p0s1 = [r[3*len(r)//4],(qty_values[-1]-qty_values[0])/8,0.5]
                    init_par = [*p0s,*p0s1]
                    first = False
                else:
                    init_par = p02_last
                    init_par = [np.abs(parameter) for parameter in init_par]
                    init_par[0]+=shift
                    init_par[2]+=shift
                    
                par1,sig1 = curve_fit(two_gaussian,qty_values,hist,p0=init_par)
                self.ax_hist.plot(qty_values,two_gaussian(qty_values,*par1),'b-')
                self.ax_hist.plot(qty_values,par1[-1]**2 *gaussian(qty_values,*par1[:2])/(par1[-1]**2 +1),'b-.')
                self.ax_hist.plot(qty_values,gaussian(qty_values,*par1[2:4])/(par1[-1]**2 +1),'b--')
                
                # self.ax_hist.plot(qty_values,two_gaussian(qty_values,*init_par),'g-')
                # self.ax_hist.plot(qty_values,init_par[-1]**2 *gaussian(qty_values,*init_par[:2])/(init_par[-1]**2 +1),'g.-')
                # self.ax_hist.plot(qty_values,gaussian(qty_values,*init_par[2:4])/(init_par[-1]**2 +1),'g--')
                p02_last = np.copy(par1)
            except Exception as e:
                print(e)

                    
            self.ax_hist.set_title(#f'chi: {np.sum(chi) :.2e}'+
                                      # f'\n mean fit - hist: {par[0] - np.mean(qty):.2e}'+
                                      f'\n skew fit - hist {skewness:.2e}')

                
                
            # except Exception as e:
            #     print(e)
                
            self.fig_hist.canvas.draw()
            
            
        

            

            r_last = np.mean(r)
            # if 0: 
            #     if two_gauss:
            #         try:
            #             self.y1_list.append(par1[0])
            #             self.y2_list.append(par1[2])
            #             self.dy1_list.append(par1[1])
            #             self.dy2_list.append(par1[3])
            #         except:
            #             self.y1_list.append(par[0])
            #             self.y2_list.append(0)
            #             self.dy1_list.append(par[1])
            #             self.dy2_list.append(0)
            #     else:
            #         self.y1_list.append(par[0])
            #         self.y2_list.append(0)
            #         self.dy1_list.append(par[1])
            #         self.dy2_list.append(0)
            keyboard = False
            while not keyboard:
                keyboard = plt.waitforbuttonpress()            
            self.y1_list.append(self.y1)
            self.y2_list.append(self.y2)
            self.y3_list.append(self.y3)
            self.dy1_list.append(self.dy1)
            self.dy2_list.append(self.dy2)
            # plt.close('all')
            amp_last = amp

        
        plt.figure()
        plt.plot(self.amps,self.y1_list)        
        plt.close(self.fig_meas)     
        
 
        
    def pick_start(self,event):
        if event.button == 2:
            self.n_click += 1
            if self.n_click%2 == 1:                
                if not self.red_line == None:
                    self.red_line.remove()
                self.x1 = event.xdata
                self.red_line = self.ax.axvline(self.x1,c='r')
            elif self.n_click%2 == 0:
                if not self.green_line == None:
                    self.green_line.remove()
                self.x2 = event.xdata
                self.green_line = self.ax.axvline(self.x2,c='g')
            self.fig.canvas.draw()
        if event.button == 3:
            if self.figopen:
                for fig_m in self.figax:
                    plt.close(fig_m[0])
                self.figax = []
                
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
                
                
                if FFT:
                    transformed = fft.fft(self.xs[i]+1j*self.ys[i])
                    freq = fft.fftfreq(len(self.rs[i]),1/self.sample_rate)
                    
                    for fq1,fq2 in [cutoff]:
                        index = np.logical_and(np.abs(freq)>fq1,np.abs(freq)<fq2)
                        transformed[index]=np.zeros_like(transformed)[index]
        
                    inverse = fft.ifft(transformed)
                    
                    self.xs[i] = np.real(inverse)
                    self.ys[i]=np.imag(inverse)
            
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
switching_ampsweeps_names = {'1.35':'T1350_2','1.45':'T1450_2','1.65':'T1650_2','1.85':'T1850'}

temps = ['1.35','1.45','1.65','1.85']
letters = ['A','B','C','D']

for letter in letters:
    print(f'Resonator: {letter}')
    if letter in ['B']:
        for temp in temps[1:2]:
            params = {}
            print(f'Temperature: {temp}')
            date = '*'
            distance = 500
               
            folder = switching_sweeps_names[temp]
            folder_ampsweeps = switching_ampsweeps_names[temp]
            
            basedir = r'D:\OneDrive_copy\OneDrive - Univerzita Karlova\DATA\2023_4_Helmholtz_resonators_recal\Helmholtz_res_drive_dep\{}\buffer'.format(folder)
            
            ampsweeps_base = r'D:\OneDrive_copy\OneDrive - Univerzita Karlova\DATA\2023_4_Helmholtz_resonators_recal\Helmholtz_res_drive_dep\{}\ampsweeps_fast'.format(folder_ampsweeps)
            
            resdir = path.join(basedir,r'resonator_{}{}'.format(distance,letter))
            files = glob(path.join(resdir, r'*.npy'))
            
            ampsweeps = path.join(ampsweeps_base,r'resonator_{}{}_updowndown'.format(distance,letter))
            files_ampsweeps = glob(path.join(ampsweeps, r'*.npy'))
            
            # files = files[:-53]
            files.sort(key=lambda f: np.load(f, allow_pickle=True).item()['App_gen (V)'])
            
            fig_ampsweeps, ax_ampsweeps = plt.subplots()
            
            cmap = plt.get_cmap('viridis')
            cmap_hist = plt.get_cmap('Reds')
        
            #%% plot ampsweeps
            for file in files_ampsweeps[0:14]:#files[:20] and files[39:]:
                # print(file)
                d = np.load(file, allow_pickle=True).item()
                try:
                    if not isinstance(d['x (A)'][0],list):    
                        f = np.array(d['freq (Hz)'])
                        x = np.array(d['x (A)'])
                        y = np.array(d['y (A)'])
                        A = np.array(d['App (V)'])
                        # A+=5.0
                        # A /= 10
                        # base = d['App_base (V)']
                        N = len(x)
                        leg1 = np.arange(0, int(N/3))
                        leg2 = np.arange(int(N/3), int(2*N/3))
                        leg3 = np.arange(int(2*N/3), N)
                        
                        ax_ampsweeps.plot(A[leg1], np.sqrt(x**2 + y**2)[leg1], '-o',ms=0.5, lw=0.5, color='tab:blue',alpha=1)
                        ax_ampsweeps.plot(A[leg2], np.sqrt(x**2 + y**2)[leg2], '-o',ms=0.5, lw=0.5, color='tab:orange',alpha=1)
                        ax_ampsweeps.plot(A[leg3], np.sqrt(x**2 + y**2)[leg3], '-o',ms=0.5, lw=0.5, color='tab:green',alpha=1)
                        
                except:
                    pass

            obj = get_drive(fig_ampsweeps,ax_ampsweeps,files)
            
            # keyboard = False
            # while not keyboard:
            #     keyboard = plt.waitforbuttonpress() 
    
            params['x1'] = obj.x1 #in normalised drive amplitude
            params['x2'] = obj.x2 
            params['y1'] = obj.y1_list
            params['y2'] = obj.y2_list
            params['dy1'] = obj.dy1_list
            params['dy2'] = obj.dy2_list
            plt.close(fig_ampsweeps)
            if save:
                np.save(f'lifetimes_params/params_{letter}_T{temp}.npy',params)












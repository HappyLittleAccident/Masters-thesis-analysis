# -*- coding: utf-8 -*-
"""
Created on Mon May  6 11:49:24 2024

@author: Marek
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from glob import glob
from scipy.fft import fft2,fftfreq


#%%



def prob_activate(p,count):
    tot = 0
    if count <= 0:
        return tot
    else:
        for l in range(count):
            tot += p*(1-tot)
    return tot





all_s = []
N = 10000
Nx = 70

X,Y= np.meshgrid(np.linspace(0,1,Nx),np.linspace(0,1,Nx))
sin = (np.sin(8*np.pi*X)*np.sin(8*np.pi*Y)+1)/2

rand = np.random.rand(Nx,Nx)
s0 = np.full((Nx,Nx),0.5)

zeros = np.full_like(s0,0)
ones = np.ones_like(s0)

plt.close('all')  
s_point = (-1,0,1)


probs_difuse = [0.2,0.2,0.2,0.2]
probs_anihilate = [0.3,0.3,0.3,0.3]
probs_create = [0.0,0.01,0.02,0.03]
probs_fluctuate = [1e-2,1e-2,1e-2,1e-2]

periodic = False

for p_difuse,p_anihilate,p_create,p_fluctuate in zip(probs_difuse,probs_anihilate,probs_create,probs_fluctuate):
    print('start')
    # s0 = np.where(sin>0.9,ones,zeros) + np.where(sin<0.1,np.full_like(s0,-1) ,zeros)
    s0 = np.where(rand>0.95,ones,zeros) + np.where(rand<0.05,np.full_like(s0,-1) ,zeros)

    update_prob = 0.95
    
    
    indices = np.array([0,1,2,3,4])
    
    
    
    times = range(N)
    n = np.zeros_like(times)
    s = np.zeros((Nx,Nx,len(times)))
    for t in times:
        updates = np.random.rand(Nx,Nx)
        for i,rows in enumerate(s0):
            for j,sij in enumerate(rows):
                if updates[i,j]>update_prob:
                    try:
                        if not periodic and (j==Nx-1 or i==Nx-1 or j==0 or i==0):
                            s0[i,j] = 0
                        else:
                            array = np.array([sij,s0[(i+1)%Nx,j],s0[i,(j+1)%Nx],s0[i-1,j],s0[i,j-1]])
                            zeros_count = np.count_nonzero(array==0)
                            plus_count = np.count_nonzero(array==1)
                            minus_count = 5-zeros_count-plus_count
                            
                            prob= np.random.rand()
                            
                            if zeros_count == 5:
                                  if prob<p_fluctuate:
                                      if i>3*Nx//4 and j<Nx//5:
                                          array[0]=1
                                      elif i<Nx//4 and j>4*Nx//5:
                                          array[0]=-1
                                            
                            elif sij == 1:
                                tot_dif= prob_activate(p_difuse,zeros_count)+(plus_count>0)
                                tot_an= prob_activate(p_anihilate,minus_count)   
                                tot_create = prob_activate(p_create,zeros_count-1)
                                if prob<tot_an:
                                    move = np.random.rand()*minus_count
                                    array[0] = 0
        
                                    array[indices[array==-1][int(move//1)]]=0
                                elif prob<tot_create:
                                    move = np.random.rand()*zeros_count
                                    array[indices[array==0][int(move//1)]]=1
                                    move = np.random.rand()*(zeros_count -1)
                                    array[indices[array==0][int(move//1)]]=-1
                                
                                elif prob<tot_dif:
                                    move = np.random.rand()*zeros_count
        
                                    array[indices[array==0][int(move//1)]]=1 
        
                                    array[0]=0   
        
        
                            elif sij==-1:
                                tot_create = prob_activate(p_create,zeros_count-1)
                                tot_dif= prob_activate(p_difuse,zeros_count)+(minus_count>0)
                                tot_an= prob_activate(p_anihilate,plus_count)
        
                                
                                if prob<tot_an:
                                    move = np.random.rand()*plus_count
                                    array[0] = 0
                                    array[indices[array==1][int(move//1)]]=0 
                                
                                elif prob<tot_create:
                                    move = np.random.rand()*(zeros_count)
                                    array[indices[array==0][int(move//1)]]=1
                                    move = np.random.rand()*(zeros_count -1)
                                    array[indices[array==0][int(move//1)]]=-1
                                
                                
                                elif prob<tot_dif:
                                    move = np.random.rand()*zeros_count
        
                                    array[indices[array==0][int(move//1)]]=-1  
                                    
                                    array[0]=0
                                 
        
                   
                                                
                            
                            [s0[i,j],s0[(i+1)%Nx,j],s0[i,(j+1)%Nx],s0[i-1,j],s0[i,j-1]]=array 
                            
                    except Exception as e:
                        # print(e)
                        pass
        s[:,:,t]=s0
        # print(t)
        
        n[t] = np.sum(np.abs(s0))
    print(n[-1])
    all_s.append(s)
    # break

print('saving..')
try:
    np.save('all_s.npy',all_s,allow_pickle=True)
    all_s = None
except:
    print('saving failed')
print('saved')
plt.close('all')
#%%
fign,axn = plt.subplots(2,2)
axn = axn.flatten()
all_s_loaded = np.load('all_s.npy',allow_pickle=True)

for i,s in enumerate(all_s_loaded):
    axn[i].set_title(i)              
    axn[i].plot(times,np.sum(np.abs(s),axis=(0,1)))
all_s_loaded = None
#%%

fig,ax = plt.subplots(2,2)

anis = []
ims = []
all_s_loaded = np.load('all_s.npy',allow_pickle=True)

ax = ax.flatten()

for i,s in enumerate(all_s_loaded):
    ims.append([])
    cmap= plt.get_cmap('bwr')
    ax[i].imshow(s[:,:,0],cmap=cmap)
    ax[i].set_title(i)
    for t in times:
        im=ax[i].imshow(s[:,:,t],cmap = cmap ,animated= True)
        ims[i].append([im])
    anis.append(animation.ArtistAnimation(fig,ims[i],interval=10,blit=True,repeat_delay = 1000)   )
    # anis[i].save(filename=rf"anims/{i}.mp4", writer="ffmpeg") 
    # print('saved')
   
all_s_loaded = None    
#%%
plt.close('all')
anis = []

ims = []
norm_fft = 200
file = 'all_s.npy'
d = np.load(file,allow_pickle=True)
print('loaded')
fig,ax = plt.subplots(2,2)
# fign,axn = plt.subplots(2,2)


ax = ax.flatten()

for i,sf in enumerate(d):
    
    ims.append([])
    cmap= plt.get_cmap('inferno')
    fft = fft2(sf[:,:,0],(Nx,Nx))
    fft = np.abs(fft)[0:Nx//2,0:Nx//2]/norm_fft
    ax[i].imshow(fft,cmap=cmap)
    ax[i].set_title(i)
    for t in times:
        fft = fft2(sf[:,:,t],(Nx,Nx))[0:Nx//2,0:Nx//2]
        fft = np.abs(fft)/norm_fft
        im=ax[i].imshow(fft,cmap = cmap ,animated= True)
        ims[i].append([im])
    anis.append(animation.ArtistAnimation(fig,ims[i],interval=10,blit=True,repeat_delay = 1000)   )
d = None
#%%
plt.close('all')
anis = []

ims = []
norm_fft = 200
file = 'all_s.npy'
d = np.load(file,allow_pickle=True)
print('loaded')

k = fftfreq(Nx)

Kx,Ky = np.meshgrid(k,k)

k = np.sqrt(Kx**2 + Ky**2).flatten()
sorted_indices = k.argsort()
k = k[sorted_indices]
class PauseAnimation:
    def __init__(self):
        fig,ax = plt.subplots(2,2)
        self.ax = ax.flatten()

        # Start with a normal distribution
        self.p = [None]*4
        self.d = d
        for j,a in enumerate(self.ax):
            

            ims.append([])
            fft = fft2(self.d[j][:,:,0],(Nx,Nx)).flatten()[sorted_indices]
            fft = np.abs(fft)/norm_fft

            self.p[j] = a.plot(k,fft,'k-')[0]
            self.ax[j].set_title(j)
            
        plt.waitforbuttonpress()
        self.animation = animation.FuncAnimation(
            fig, self.update, frames=N, interval=50, blit=True)
        # self.animation.pause()
        self.paused = False
        
        fig.canvas.mpl_connect('button_press_event', self.toggle_pause)

    def toggle_pause(self, *args, **kwargs):
        if self.paused:
            self.animation.resume()
        else:
            self.animation.pause()
        self.paused = not self.paused

    def update(self, i):
        for j,a in enumerate(self.p):
            fft = fft2(self.d[j][:,:,i],(Nx,Nx)).flatten()[sorted_indices]
            fft = np.abs(fft)/norm_fft
            a.set_ydata(fft)
            self.ax[j].relim()
            self.ax[j].autoscale_view()
        return self.p


pa = PauseAnimation()
plt.show()

d = None
   
#%%
files = glob(r'D:\Github\Masters-thesis-analysis\simulations\simulations\*.npy')
anis = []

ims = []
norm_fft = 200
cmap = plt.get_cmap('viridis')
for file in files:
    d = np.load(file,allow_pickle=True)
    print('loaded')
    fig,ax = plt.subplots(2,2)
    # fign,axn = plt.subplots(2,2)


    ax = ax.flatten()

    for i,sf in enumerate(d):
        
        ims.append([])
        cmap= plt.get_cmap('inferno')
        fft = fft2(sf[:,:,0],(Nx,Nx))[0:len(fft)//2,0:len(fft)//2]
        fft = np.abs(fft)/norm_fft
        ax[i].imshow(fft,cmap=cmap)
        ax[i].set_title(i)
        for t in times:
            fft = fft2(sf[:,:,t],(Nx,Nx))
            fft = np.abs(fft)/norm_fft
            im=ax[i].imshow(fft,cmap = cmap ,animated= True)
            ims[i].append([im])
        anis.append(animation.ArtistAnimation(fig,ims[i],interval=10,blit=True,repeat_delay = 1000)   )
        # anis[i].save(filename=rf"anims/{i}.mp4", writer="ffmpeg") 
        # print('saved')
        # axn[i].set_title(i)              
        # axn[i].plot(times,np.sum(np.abs(sf),axis=(0,1)))
    
    
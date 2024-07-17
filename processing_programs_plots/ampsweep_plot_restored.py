# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 16:55:10 2023

@author: Marek
"""

import numpy as np
import matplotlib.pyplot as plt
import os
from glob import glob
from format_figure import use_latex,polish
import matplotlib as mpl

plt.close('all')
folder = r'Helmholtz_res_drive_dep'
temps = ['T1325', 'T1350', 'T1375', 'T1400', 'T1425', 'T1450', 'T1475',
         'T1500', 'T1550', 'T1650', 'T1700', 'T1750', 'T1850', 'T1950']
resonators = ['500A', '500B', '500C']
names = ['C', 'J', 'R']

use_latex()
# plt.rcdefaults()

cmap = plt.get_cmap('winter_r')
viridis = plt.get_cmap('viridis')

# plt.ioff()
plt.ion()

for i,resonator in enumerate(resonators[0:3]):
    units= 1 if resonator == '500C' else 1e-3
    fig,ax = plt.subplots(3,1,sharex=False,dpi=250,gridspec_kw={'height_ratios': [1, 1,0.05]})
    for temp in temps[0::]:
        if (temp in ['T1325','T1650','T1850'] and i<2) or (i==2 and temp in ['T1325','T1550','T1700']):
            scale = (float(temp[1:])*1e-3 - 1.325)/(1.95 - 1.325)
    
            filenames = glob(os.path.join(folder,f'{temp}','ampsweeps_fast',f'resonator_{resonator}_updowndown','*.npy'))
            # A_max = np.max([np.load(file_temp,allow_pickle=True).item()['App (V)'] for file_temp in filenames])
            for file in filenames[5:6]:
                data = np.load(file,allow_pickle=True).item()
                   
                upP = data['upPressure (Pa/m)']*units
                downP = data['downPressure (Pa/m)']*units
                jumpP = data['jumpPressure (Pa/m)']*units
                up = data['upVelocity (m/s)']*100
                down =  data['downVelocity (m/s)']*100
                jump = data['jumpVelocity (m/s)']*100

                if resonator == '500C':
                    upP/=25000
                    downP/=25000
                    up/=11
                    down/=11
                    

                ax[0].plot(upP,up,'o-',ms=0.5,lw=0.5,c='tab:blue',alpha=1)
                ax[0].plot(downP,down,'o-',ms=0.5,lw=0.2,c='tab:orange',alpha=0.2)
                # ax[0].arrow(maxx/2, maxy/4, maxx/10,maxy/10,width=maxy/30,head_length = maxx/30)
                
                
                if i==0:
                    r =0.15
                    angle_up = 60
                    angle_down = 180+60
                    
                    xy_up = (0.15, 0.30)
                    xy_end_up = (xy_up[0]+np.cos(angle_up*np.pi/180)*r, xy_up[1]+np.sin(angle_up*np.pi/180)*r)
                    
                    xy_down = (0.25,0.32)
                    xy_end_down = (xy_down[0]+np.cos(angle_down*np.pi/180)*r, xy_down[1]+np.sin(angle_down*np.pi/180)*r)
                    
                    ax[0].annotate("", xytext=xy_up,xycoords='axes fraction',
                           xy=xy_end_up,arrowprops=dict(arrowstyle="simple",color='tab:blue'))
                    ax[0].annotate("", xy=xy_end_down,xycoords='axes fraction',
                           xytext=xy_down,arrowprops=dict(arrowstyle="simple",color='tab:orange'))
                
                elif i==1:
                    r =0.15
                    angle_up = 68
                    angle_down = 180+68
                    
                    xy_up = (0.12, 0.33)
                    xy_end_up = (xy_up[0]+np.cos(angle_up*np.pi/180)*r, xy_up[1]+np.sin(angle_up*np.pi/180)*r)
                    
                    xy_down = (0.20,0.35)
                    xy_end_down = (xy_down[0]+np.cos(angle_down*np.pi/180)*r, xy_down[1]+np.sin(angle_down*np.pi/180)*r)
                    
                    ax[0].annotate("", xytext=xy_up,xycoords='axes fraction',
                           xy=xy_end_up,arrowprops=dict(arrowstyle="simple",color='tab:blue'))
                    ax[0].annotate("", xy=xy_end_down,xycoords='axes fraction',
                           xytext=xy_down,arrowprops=dict(arrowstyle="simple",color='tab:orange'))
                
                else:
                    r =0.15
                    angle_up = 60
                    angle_down = 180+60
                    
                    xy_up = (0.15, 0.28)
                    xy_end_up = (xy_up[0]+np.cos(angle_up*np.pi/180)*r, xy_up[1]+np.sin(angle_up*np.pi/180)*r)
                    
                    xy_down = (0.25,0.30)
                    xy_end_down = (xy_down[0]+np.cos(angle_down*np.pi/180)*r, xy_down[1]+np.sin(angle_down*np.pi/180)*r)
                    
                    ax[0].annotate("", xytext=xy_up,xycoords='axes fraction',
                           xy=xy_end_up,arrowprops=dict(arrowstyle="simple",color='tab:blue'))
                    ax[0].annotate("", xy=xy_end_down,xycoords='axes fraction',
                           xytext=xy_down,arrowprops=dict(arrowstyle="simple",color='tab:orange'))
                
                maxx = ax[0].get_xlim()[1]-ax[0].get_xlim()[0]
                maxy = ax[0].get_ylim()[1]-ax[0].get_ylim()[0]
                z = (upP[-1] - upP[upP > upP[-1]*0.8][0])/maxx + 1j*(up[-1] - up[upP > upP[-1]*0.8][0])/maxy               
                if i == 2 and temp=='T1700':
                    ax[0].text(upP[-1]*0.8,up[upP > upP[-1]*0.8][0]-maxy*0.07,temp[1:2]+'.'+temp[2:]+'K',rotation=np.angle(z,deg=True)-5)
                elif i == 2 and temp=='T1550':
                    ax[0].text(upP[-1]*0.8,up[upP > upP[-1]*0.8][0]+maxy*0.015,temp[1:2]+'.'+temp[2:]+'K',rotation=np.angle(z,deg=True)-5)
                    
                else:
                    ax[0].text(upP[-1]*0.65,up[upP > upP[-1]*0.65][0]+maxy*0.05,temp[1:2]+'.'+temp[2:]+'K',rotation=np.angle(z,deg=True)+
                               (temp=='T1650')*5 - (temp=='T1850')*5) 
                    
    for j,temp in enumerate(temps[0::]):
        scale = (float(temp[1:])*1e-3 - 1.325)/(1.95 - 1.325)

        filenames = glob(os.path.join(folder,f'{temp}','ampsweeps_fast',f'resonator_{resonator}_updowndown','*.npy'))
        # A_max = np.max([np.load(file_temp,allow_pickle=True).item()['App (V)'] for file_temp in filenames])
        vs = []
        ps = []    
        for file in filenames[:]:
            data = np.load(file,allow_pickle=True).item()
            
            indent= j*(700+ (i==1)*5000 + (i==2)*1500)
            upP = data['upPressure (Pa/m)'] + indent
            downP = data['downPressure (Pa/m)'] + indent
            jumpP = data['jumpPressure (Pa/m)'] + indent
            upV = data['upVelocity (m/s)']*100 
            downV =  data['downVelocity (m/s)']*100
            jumpV = data['jumpVelocity (m/s)']*100  
            
            
            vsend=5
            if resonator == '500C':
                upP/=25000
                downP/=25000
                jumpP/=25000
                indent/=25000
                vsend = vsend/20
                upV/=11
                downV/=11
                jumpV/=11
                if j>len(temps)-3:
                    continue
            
            

            vs.append(downV)
            vs.append(jumpV)
            # print(temp)

            ps.append(downP)
            ps.append(jumpP)
            ax[1].plot(upP*units,upV,'.',markeredgewidth=0.0,ms=0.5,c=cmap(scale),alpha=0.2)
            ax[1].plot(downP*units,downV,'.',markeredgewidth=0.0,ms=0.5,c=cmap(scale),alpha=0.2)
            ax[1].plot([indent*units,indent*units],[0,vsend],'--',lw=0.5,c=cmap(scale))
            ax[1].set_ylim(0,None)
            
    if i==2:
        ax[0].set_xlabel('Normalised pressure gradient (Arb.)')
        ax[1].set_xlabel('Normalised pressure gradient (Arb.)')
        ax[0].set_ylabel('Normalised superfluid velocity (Arb.)')
        ax[1].set_ylabel('Normalised superfluid velocity (Arb.)')
    else:
        ax[0].set_xlabel('Pressure gradient (Pa/mm)')
        ax[1].set_xlabel('Pressure gradient (Pa/mm)')
        ax[0].set_ylabel('Superfluid velocity (cm/s)')
        ax[1].set_ylabel('Superfluid velocity (cm/s)')
    
    
    fig.colorbar(plt.cm.ScalarMappable(norm=mpl.colors.Normalize(
    vmin=1.325, vmax=1.95, clip=False), cmap=cmap),orientation='horizontal',location="bottom", ax=None,cax=ax[2], label='Temperature (K)')

    # for j,temp in enumerate(temps[0::]):
    #     scale = (float(temp[1:])*1e-3 - 1.325)/(1.95 - 1.325)

    #     filenames = glob(os.path.join(folder+'_old_calib',f'{temp}','ampsweeps_fast',f'resonator_{resonator}_updowndown','*.npy'))
    #     # A_max = np.max([np.load(file_temp,allow_pickle=True).item()['App (V)'] for file_temp in filenames])
    #     vs = []
    #     ps = []    
    #     for file in filenames[:]:
    #         data = np.load(file,allow_pickle=True).item()
            
    #         indent= j*(700+ (i==1)*5000 + (i==2)*1500)
    #         upP = data['upPressure (Pa/m)'] + indent
    #         downP = data['downPressure (Pa/m)'] + indent
    #         jumpP = data['jumpPressure (Pa/m)'] + indent
    #         upV = data['upVelocity (m/s)']*100 
    #         downV =  data['downVelocity (m/s)']*100
    #         jumpV = data['jumpVelocity (m/s)']*100  
            
            
    #         if resonator == '500C':
    #             upP/=25000
    #             downP/=25000
    #             jumpP/=25000
    #             upV/=11
    #             downV/=11
    #             jumpV/=11
    #             if j>len(temps)-3:
    #                 continue
            
            

    #         vs.append(downV)
    #         vs.append(jumpV)
    #         # print(temp)

    #         ps.append(downP)
    #         ps.append(jumpP)
    #         ax[1].plot(upP*units,upV,'.',markeredgewidth=0.0,ms=0.5,c=viridis(scale),alpha=0.2)
    #         ax[1].plot(downP*units,downV,'.',markeredgewidth=0.0,ms=0.5,c=viridis(scale),alpha=0.2)
            
    # if i==2:
    #     ax[0].set_xlabel('Normalised pressure gradient (Arb.)')
    #     ax[1].set_xlabel('Normalised pressure gradient (Arb.)')
    #     ax[0].set_ylabel('Normalised superfluid velocity (Arb.)')
    #     ax[1].set_ylabel('Normalised superfluid velocity (Arb.)')
    # else:
    #     ax[0].set_xlabel('Pressure gradient (Pa/mm)')
    #     ax[1].set_xlabel('Pressure gradient (Pa/mm)')
    #     ax[0].set_ylabel('Superfluid velocity (cm/s)')
    #     ax[1].set_ylabel('Superfluid velocity (cm/s)')




    # # 

    # # polish(fig, 1, name=f'images//ampsweeps{names[i]}', extension='.png', grid=True,width_to_height = 1.3)        
            


    # # fig.colorbar(plt.cm.ScalarMappable(
    # # norm=mpl.colors.Normalize(vmin=1.325, vmax=1.85, clip=False), cmap=cmap),
    # # ax=ax, label='Temperature (K)')


    polish(fig, 1, name=f'images//ampsweeps{names[i]}', extension='.png', grid=True,width_to_height = 0.5,tight_layout = True)        



#%%
# =============================================================================
# plt.close('all')
# for i,resonator in enumerate(resonators[0:]):
#     fig,ax = plt.subplots(1,1,sharex=True,dpi=200)
#     for j,temp in enumerate(temps[0::]):
#         scale = (float(temp[1:])*1e-3 - 1.325)/(1.95 - 1.325)
# 
#         filenames = glob(os.path.join(folder,f'{temp}','ampsweeps_fast',f'resonator_{resonator}_updowndown','*.npy'))
#         # A_max = np.max([np.load(file_temp,allow_pickle=True).item()['App (V)'] for file_temp in filenames])
#         vs = []
#         ps = []    
#         for file in filenames[:]:
#             data = np.load(file,allow_pickle=True).item()
#             
#             indent= j*(700+ (i==1)*5000 + (i==2)*1500)
#             upP = data['upPressure (Pa/m)'] + indent
#             downP = data['downPressure (Pa/m)'] + indent
#             jumpP = data['jumpPressure (Pa/m)'] + indent
#             upV = data['upVelocity (m/s)']*100 
#             downV =  data['downVelocity (m/s)']*100
#             jumpV = data['jumpVelocity (m/s)']*100  
#             
#             if resonator == '500C':
#                 upP/=25000
#                 downP/=25000
#                 jumpP/=25000
#                 upV/=11
#                 downV/=11
#                 jumpV/=11
#                 if j>len(temp)-2:
#                     continue
#             
#             
# 
#             vs.append(downV)
#             vs.append(jumpV)
#             
# 
#             ps.append(downP)
#             ps.append(jumpP)
#             ax.plot(upP,upV,'.',markeredgewidth=0.0,ms=0.5,c=cmap(scale),alpha=0.2)
#             ax.plot(downP,downV,'.',markeredgewidth=0.0,ms=0.5,c=cmap(scale),alpha=0.2)
# 
# 
#     ax.set_xlabel('Pressure gradient (Pa/m)')
#     ax.set_ylabel('Superfluid velocity (cm/s)')
#     fig.tight_layout()
#     fig.colorbar(plt.cm.ScalarMappable(norm=mpl.colors.Normalize(
#     vmin=1.325, vmax=1.85, clip=False), cmap=cmap), ax=ax,orientation='horizontal', label='Temperature (K)')
# 
# 
#     polish(fig, 1, name=f'images//ampsweeps{names[i]}', extension='.png', grid=True,width_to_height = 1.3)        
# 
# 
# =============================================================================
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 16:45:42 2022

@author: TsapatsarisLeonidas
"""

import SAXS
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from time import process_time
plt.rcParams.update({'font.family':'serif'})
df = pd.DataFrame(np.array([[0,1,2],[3,4,5],[6,7,8]]),index = ['one','two','three'],columns = [0.1,0.2,0.3])
scdata = pd.read_csv('80C scattering data.csv',header=0,index_col=0,dtype=np.float64)
# print(scdata.columns[20])
array = np.empty(scdata.shape[1])
for i in range(0,scdata.shape[1]):
    array[i] = float(scdata.columns[i])
# print(type(array[0]))
# print(type(scdata.columns[-1]))
# print(type(scdata.index[-1]))

## Loop Example ###############################################################
files = ['RT scattering data.csv','80C scattering data.csv']
invariant_results = []
for i in range(0,2):
    dataEX = SAXS.ScatteringData(data=None,csv=True,csv_path = files[i])
    invariant=dataEX.invariant(lowlim=0,upplim=scdata.shape[1])
    invariant_results.append(invariant)
#################################################################################
print(invariant_results)
plt.plot(invariant_results[0][0],invariant_results[0][1],marker = '^',linestyle='',markeredgecolor='k',markerfacecolor = 'blue',label='RT')
plt.plot(invariant_results[1][0],invariant_results[1][1],marker = '^',linestyle='',markeredgecolor='k',markerfacecolor='red',label='80C')
plt.legend()
plt.xlabel('Aging Time [min]')
plt.ylabel('Q$^*$')
data80 = SAXS.ScatteringData(data=None,csv=True,csv_path = '80C scattering data.csv')
intensity = data80.intensity()
invariant = data80.invariant(lowlim = 0,upplim=scdata.shape[1])
t1_start = process_time() 
IQmax,cl = data80.IQmax(lowlim=0,upplim=500,char_length=True)
fit,params = data80.curve_fit(lowlim=100,upplim=300,fitfunc = 'Guinier',parameters=True)
q,iq,tage = data80.get_data(aging_time=True)
t1_stop = process_time()
print("Elapsed time:", t1_stop-t1_start,' seconds') 
fig2,ax2 = plt.subplots(1,2,figsize=(8,6))
ax2[0].plot(invariant[:,0],invariant[:,1],marker='o',markersize=12,markeredgecolor='k',markerfacecolor='green',linestyle='')
ax2[0].set_title('Invariant',fontsize=14)
ax2[1].plot(cl[:,0],cl[:,1],marker='s',markersize=12,markeredgecolor='k',markerfacecolor = 'green',linestyle='')
ax2[0].grid('on')
ax2[1].grid('on')
ax2[0].set_xlabel('Time [min]',fontsize=12)
ax2[1].set_xlabel('Time [min]',fontsize=12)
ax2[0].set_ylabel('Q$^*$',fontsize=12)
ax2[1].set_ylabel('[$\AA$]',fontsize=12)
ax2[1].set_title('Characteristic Length',fontsize=14)
fig2.tight_layout()
fig2.savefig('invariant and cl.png',facecolor=fig2.get_facecolor(),edgecolor=None,dpi=300)
fig,ax = plt.subplots(1,2,figsize = (8,6))
for i in range(tage.shape[0]-4,tage.shape[0]):
    if i < tage.shape[0]-3:
    # if i == 0:
        ax[0].loglog(q,iq[i,:]*2,marker='^',markersize = 2)
        ax[1].loglog(intensity.columns,intensity.iloc[i,:]*2,marker='^',markersize=2)
        ax[0].loglog(q[100:300], fit[i,:]*2,'k--',label = 'Guinier Fit, R$_g$ = {:.2f}'.format(params[i,2]))
    # else:
    #     ax[0].loglog(q,iq[i,:]*i,marker='^',markersize = 2)
    #     ax[1].loglog(intensity.columns,intensity.iloc[i,:]*i,marker='^',markersize=2)
    # elif i < tage.shape[0]-2:
    #     ax[0].loglog(q,iq[i,:]*2,marker='^',markersize = 2)
    #     ax[1].loglog(intensity.columns,intensity.iloc[i,:],marker='^',markersize=2)
    elif i < tage.shape[0]-2:
        ax[0].loglog(q,iq[i,:]*3,marker='^',markersize=2)
        ax[1].loglog(intensity.columns,intensity.iloc[i,:]*3,marker='^',markersize=2)
        ax[0].loglog(q[100:300], fit[i,:]*3,'k--',label = 'Guinier Fit, R$_g$ = {:.2f}'.format(params[i,2]))
    else:
        ax[0].loglog(q,iq[i,:]*4,marker='^',markersize=2)
        ax[1].loglog(intensity.columns,intensity.iloc[i,:]*4,marker='^',markersize=2)
        ax[0].loglog(q[100:300], fit[i,:]*4,'k--',label = 'Guinier Fit, R$_g$ = {:.2f}'.format(params[i,2]))
#plt.legend()
ax[0].grid('on')
ax[1].grid('on')
ax[0].set_xlabel('Q [$\AA^{-1}$]',fontsize=12)
ax[1].set_xlabel('Q [$\AA^{-1}$]',fontsize=12)
ax[0].set_ylabel('I(Q) [a.u.]',fontsize=12)
ax[1].set_ylabel('I(Q)*Q$^2$ [a.u.]',fontsize=12)
ax[0].set_title('Scattering Data',fontsize=14)
ax[1].set_title('Intensity',fontsize=14)
fig.tight_layout()
#fig.savefig('Scattering and Intensity Subplot w fits.png', facecolor=fig.get_facecolor(),edgecolor='None',dpi=300)
#plt.plot(IQmax[:,0],IQmax[:,1])
#fig,ax = plt.subplots(figsize=(8,6))
#ax.plot(cl[:,0],cl[:,1],'ko',markersize = 12)
# fig1,ax1 = plt.subplots(figsize=(8,6))
# for i in range(intensity.shape[0]):
#       ax1.plot(intensity.columns,intensity.iloc[i])
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Set up dataframe with even samples
Rolling windowed xcorr

Created on Tue Sep  7 20:25:30 2021

@author: lizz
"""
import numpy as np
import pandas as pd
import seaborn as sns

point_to_plot=5
t_start = max(min(d_interp), min(smb_d_interp), min(tm_d_interp))
t_end = min(max(d_interp), max(smb_d_interp), max(tm_d_interp))

df = pd.DataFrame(index=t_grid[t_grid<=t_end], columns=['Velocity','SMB','Runoff','Terminus'])
for i, d in enumerate(t_grid[t_grid<=t_end]):
    df.loc[d] = [preds[point_to_plot]['full'][i],
        smb_func(d),
        runoff_func(d),
        termini_func(d)]
df = df.diff()

## difference and normalize the dataframe
for c in df.columns:
    if c=='Velocity':
        df[c] = (df[c] - df[c].mean())/df[c].std()
    else:
        df[c] = (df[c] - df[c].mean())/(df[c].std()*len(df))
df = df.dropna(axis=0, how='all')


## rolling cross-corr
spacing = 3.5
half_window_size = 52 #samples
t_start = 0
t_end = t_start + window_size
step_size = 5
rss_t=[]
rss_s=[]
rss_r=[]
while t_end < len(df):
    d1 = df['Velocity'].iloc[t_start:t_end]
    d2 = df['Terminus'].iloc[t_start:t_end]
    d3 = df['SMB'].iloc[t_start:t_end]
    d4 = df['Runoff'].iloc[t_start:t_end]
    rs_t = np.correlate(d2, d1, mode='full')
    rs_s = np.correlate(d3, d1, mode='full')
    rs_r = np.correlate(d4, d1, mode='full')    
    rss_t.append(rs_t)
    rss_s.append(rs_s)
    rss_r.append(rs_r)
    t_start = t_start + step_size
    t_end = t_end + step_size
rss_t = pd.DataFrame(rss_t)
rss_s = pd.DataFrame(rss_s)
rss_r = pd.DataFrame(rss_r)

## plot with seaborn heatmap
f,(ax1,ax2,ax3) = plt.subplots(1,3,figsize=(10,10), sharey=True)
sns.heatmap(rss_s,cmap='RdBu_r',ax=ax1,
            vmin=-0.1, vmax=0.1, cbar=False)
ax1.set(title='SMB', ylabel='Epochs')
sns.heatmap(rss_r,cmap='RdBu_r',ax=ax2,
            vmin=-0.1, vmax=0.1, cbar=False)
ax2.set(title='Runoff')
sns.heatmap(rss_t,cmap='RdBu_r',ax=ax3,
            vmin=-0.1, vmax=0.1, cbar=True)
ax3.set(title='Terminus')

for ax in (ax1,ax2,ax3):
    ax.axvline(x=52, ls=':', lw=1.0, color='k')
    ax.set(xlim=[0,104], xlabel='Offset')
    ax.set_xticks([0, 26, 52, 78, 104])
    ax.set_xticklabels([-180, -90, 0, 90, 180]);
    

## limit to positive lags (potentially causal)
f,(ax1,ax2,ax3) = plt.subplots(1,3,figsize=(10,10), sharey=True)
sns.heatmap(rss_s,cmap='RdBu_r',ax=ax1,
            vmin=-0.1, vmax=0.1, cbar=False)
ax1.set(title='SMB', ylabel='Epochs')
sns.heatmap(rss_r,cmap='RdBu_r',ax=ax2,
            vmin=-0.1, vmax=0.1, cbar=False)
ax2.set(title='Runoff')
sns.heatmap(rss_t,cmap='RdBu_r',ax=ax3,
            vmin=-0.1, vmax=0.1, cbar=True)
ax3.set(title='Terminus')

for ax in (ax1,ax2,ax3):
    # ax.axvline(x=52, ls=':', lw=1.0, color='k')
    ax.set(xlim=[52,104], xlabel='Offset')
    ax.set_xticks([52, 78, 104])
    ax.set_xticklabels([0, 90, 180]);

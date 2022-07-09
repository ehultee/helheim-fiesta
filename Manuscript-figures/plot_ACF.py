#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot ACFs of all variables
Created on Fri Jun 25 14:51:18 2021

@author: lizz
"""

import statsmodels.api as sm
from statsmodels.graphics.tsaplots import plot_acf

a_vel_raw = sm.tsa.stattools.acf(preds[0]['full'])[1]
a_vel = sm.tsa.stattools.acf(np.diff(preds[0]['full']))[1]
a_vel_lt = sm.tsa.stattools.acf(preds[0]['secular']+preds[0]['transient'])[1]
b_runoff_raw = sm.tsa.stattools.acf(runoff)[1]
b_runoff = sm.tsa.stattools.acf(np.diff(runoff))[1]
b_runoff_lt = sm.tsa.stattools.acf(rf_lowfreq(t_grid))[1]
F_runoff_raw = np.sqrt((1+(a_vel_raw*b_runoff_raw))/(1-(a_vel_raw*b_runoff_raw)))
F_runoff = np.sqrt((1+(a_vel*b_runoff))/(1-(a_vel*b_runoff)))
F_runoff_lt = np.sqrt((1+(a_vel_lt*b_runoff_lt))/(1-(a_vel_lt*b_runoff_lt)))

## modify CCSL for terminus
b_terminus_raw = sm.tsa.stattools.acf(td)[1]
b_terminus = sm.tsa.stattools.acf(np.diff(td))[1]
b_terminus_lt = sm.tsa.stattools.acf(tf_lowfreq(t_grid))[1]
F_terminus_raw = np.sqrt((1+(a_vel_raw*b_terminus_raw))/(1-(a_vel_raw*b_terminus_raw)))
F_terminus = np.sqrt((1+(a_vel*b_terminus))/(1-(a_vel*b_terminus)))
F_terminus_lt = np.sqrt((1+(a_vel_lt*b_terminus_lt))/(1-(a_vel_lt*b_terminus_lt)))


## modify CCSL for smb
b_smb_raw = sm.tsa.stattools.acf(smb)[1]
b_smb = sm.tsa.stattools.acf(np.diff(smb))[1]
b_smb_lt = sm.tsa.stattools.acf(smb_lowfreq(t_grid))[1]
F_smb_raw = np.sqrt((1+(a_vel_raw*b_smb_raw))/(1-(a_vel_raw*b_smb_raw)))
F_smb = np.sqrt((1+(a_vel*b_smb))/(1-(a_vel*b_smb)))
F_smb_lt = np.sqrt((1+(a_vel_lt*b_smb_lt))/(1-(a_vel_lt*b_smb_lt)))

## Plot ACFs
dt = np.mean(np.diff(t_grid))*365.26 ## convert lags to physical units of time

fig, axs = plt.subplots(4,3, sharey=True, sharex=True, figsize=(8,4))
vars_toplot = (preds[0]['full'], np.diff(preds[0]['full']), preds[0]['full']-preds[0]['seasonal'],
               smb, np.diff(smb), smb_lowfreq(t_grid),
               runoff, np.diff(runoff), rf_lowfreq(t_grid),
               td, np.diff(td), tf_lowfreq(t_grid))
names = ('', '', 'Velocity',
         '', '', 'SMB',
         '', '', 'Runoff',
         '', '', 'Terminus')
factors = (F_smb_raw, F_smb, F_smb_lt,
           F_runoff_raw, F_runoff, F_runoff_lt,
           F_terminus_raw, F_terminus, F_terminus_lt)
for var,ax, n in zip(vars_toplot, axs.ravel(), names):
    plot_acf(var, ax, lags=30)
    ax.text(1.1, 0.5, str(n), transform=ax.transAxes)
    xt = ax.get_xticks()
    ax.set_xticklabels(['{:1.0f}'.format(x*dt) for x in xt])
    ax.set(title='')
for ax in axs[:,0]:
    ax.set(ylabel='Acorr')
for i, ax in enumerate(axs[1::, :].ravel()):
    ax.text(0.73, 0.85, 'F={:5.2f}'.format(factors[i]), transform=ax.transAxes)
for ax in axs[3,:]:
    ax.set(xlabel='Lags [d]')
axs[0,0].set(title='Full')
axs[0,1].set(title='Single-diff')
axs[0,2].set(title='Long-term')
plt.tight_layout()
plt.show()
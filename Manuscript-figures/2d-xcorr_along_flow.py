#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot xcorr along flow
Created on Thu Feb  4 21:10:31 2021

@author: lizz
"""
import nifl_helper


## Extract bed topo along flowline
xyvals = np.array([(xh[i], yh[i]) for i in range(len(xh))])
bed_vals = [float(B_helheim(xh[i], yh[i])) for i in range(len(xh))]
surface_vals = [float(S_helheim(xh[i], yh[i])) for i in range(len(xh))]
xvals = (0.001*np.array(nifl_helper.ArcArray(xyvals)))

## Create arclength array for selected points
s = 0.001*np.array(nifl_helper.ArcArray(np.array(xys)))+2 # trim the 2km that was removed

## Plot

fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, constrained_layout=True)

ax1.axhline(0)
ax1.plot(s, smb_corr_amax, marker='.', color='k', ls=':', label='SMB')
ax1.plot(s, runoff_corr_amax, marker='*', color='k', ls='--', label='Runoff')
ax1.plot(s, terminus_corr_amax, marker='d', color='k', label='Terminus')
ax1.legend(loc='lower left')
ax1.set(ylabel='AMax cross-correlation')

## add long-term?
ax1.plot(s, smb_lt_corr_amax, marker='.', color='DarkGrey', ls=':', label='SMB long-term')
ax1.plot(s, rf_lt_corr_amax, marker='*', color='DarkGrey', ls='--', label='Runoff long-term')
ax1.plot(s, term_lt_corr_amax, marker='d', color='DarkGrey', label='Terminus long-term')

ax2.plot(xvals, bed_vals, color='saddlebrown')
ax2.plot(xvals, surface_vals, color='darkgrey')
plt.fill_between(xvals, surface_vals, bed_vals, color='darkgrey', alpha=0.5)
plt.fill_between(xvals, bed_vals, y2=-1300, color='saddlebrown', alpha=0.5, hatch='/')
ax2.set(xlim=(25, 0), ylim=(-1300, 1500), aspect=0.002, 
        xlabel='Upstream distance [km]', ylabel='Elevation [m a.s.l.]')

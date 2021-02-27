#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot xcorr along flow
Created on Thu Feb  4 21:10:31 2021

@author: lizz
"""
import nifl_helper
from matplotlib import cm


## Extract bed topo along flowline
xyvals = np.array([(xh[i], yh[i]) for i in range(len(xh))])
bed_vals = [float(B_helheim(xh[i], yh[i])) for i in range(len(xh))]
surface_vals = [float(S_helheim(xh[i], yh[i])) for i in range(len(xh))]
xvals = (0.001*np.array(nifl_helper.ArcArray(xyvals)))

## Create arclength array for selected points
s = 0.001*np.array(nifl_helper.ArcArray(np.array(xys)))+2 # trim the 2km that was removed


### Plot all on one set of axes
smb_offset=4000 # how many 'meters' above topo to plot this line
runoff_offset=5000 # move runoff to top to accom. vertical indicator of sign switch
terminus_offset=3000
scaling=800  # vertical multiplicative factor to display xcorr on same axes as topo
qual_colors = cm.get_cmap('tab20b')

fig, ax = plt.subplots(1, constrained_layout=True)

# ax.axvline(10, color='k', lw=1.0, ls='--')
ax.axvline(14, color='k', lw=1.0, ls='--')
ax.plot(10*np.ones_like(s), np.linspace(-1300, smb_offset, len(s)), color='k', lw=1.0, ls='--')
ax.plot(s, smb_offset*np.ones_like(s), color='k', lw=0.5, ls=':')
ax.plot(s, runoff_offset*np.ones_like(s), color='k', lw=0.5, ls=':')
ax.plot(s, terminus_offset*np.ones_like(s), color='k', lw=0.5, ls=':')
ax.plot(s, smb_offset+ scaling*np.array(smb_corr_amax), color=qual_colors(0), lw=2.0, label='SMB')
ax.plot(s, runoff_offset+ scaling*np.array(runoff_corr_amax), color=qual_colors(2), lw=2.0, label='Runoff')
ax.plot(s, terminus_offset+ scaling*np.array(terminus_corr_amax), color=qual_colors(4), lw=2.0, label='Terminus')
# ax.plot(s, smb_offset+ scaling*np.array(smb_corr_amax), color='k', lw=2.0, label='SMB')
# ax.plot(s, runoff_offset+ scaling*np.array(runoff_corr_amax), color='k', lw=2.0, label='Runoff')
# ax.plot(s, terminus_offset+ scaling*np.array(terminus_corr_amax), color='k', lw=2.0, label='Terminus')

## add long-term
ax.plot(s, smb_offset+ scaling*np.array(smb_lt_corr_amax), color=qual_colors(1), ls=':', lw=2.0, label='SMB long-term')
ax.plot(s, runoff_offset+ scaling*np.array(rf_lt_corr_amax), color=qual_colors(3), ls=':', lw=2.0, label='Runoff long-term')
ax.plot(s, terminus_offset+ scaling*np.array(term_lt_corr_amax), color=qual_colors(5), ls=':', lw=2.0, label='Terminus long-term')

## fill to show sign
plt.fill_between(s, smb_offset+ scaling*np.array(smb_corr_amax), y2=smb_offset, 
                 color=qual_colors(0), alpha=0.4, hatch='|')
plt.fill_between(s, runoff_offset+ scaling*np.array(runoff_corr_amax), y2=runoff_offset, 
                 color=qual_colors(2), alpha=0.4, hatch='|')
plt.fill_between(s, terminus_offset+ scaling*np.array(terminus_corr_amax), y2=terminus_offset, 
                 color=qual_colors(4), alpha=0.4, hatch='|')

ax.plot(xvals, bed_vals, color='saddlebrown')
ax.plot(xvals, surface_vals, color='darkgrey')
plt.fill_between(xvals, surface_vals, bed_vals, color='darkgrey', alpha=0.5)
plt.fill_between(xvals, bed_vals, y2=-1300, color='saddlebrown', alpha=0.5, hatch='/')
ax.set(xlim=(25, 0), ylim=(-1300, 5500), aspect=0.002, 
        xlabel='Upstream distance [km]',
        yticks=(-1000, 0, 1000, terminus_offset, runoff_offset, smb_offset),
        yticklabels=('-1000', '0', '1000 m a.s.l.', 'Terminus', 'Runoff', 'SMB'))
ax.get_yticklabels()[3].set_color(qual_colors(4)) # match the ticks to the plotted color
ax.get_yticklabels()[4].set_color(qual_colors(2))
ax.get_yticklabels()[5].set_color(qual_colors(0))


ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot xcorr along flow
Created on Thu Feb  4 21:10:31 2021

@author: lizz
"""
import nifl_helper
from matplotlib import cm
from matplotlib.ticker import FixedLocator

## set matplotlib font size defaults
SMALL_SIZE = 10
MEDIUM_SIZE = 11
BIGGER_SIZE = 14

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title



## Extract bed topo along flowline
xyvals = np.array([(xh[i], yh[i]) for i in range(len(xh))])
bed_vals = [float(B_helheim(xh[i], yh[i])) for i in range(len(xh))]
surface_vals = [float(S_helheim(xh[i], yh[i])) for i in range(len(xh))]
xvals = (0.001*np.array(nifl_helper.ArcArray(xyvals)))


## Compute mean ice speed for each point
mean_speed = [np.mean(pred['full']) for pred in preds]
scaled_speed = [m/np.mean(mean_speed) for m in mean_speed]
scaled_ticks = (np.array([4.5, 6, 7.5])/np.mean(mean_speed)) - np.mean(scaled_speed) # for display


## Create arclength array for selected points
s = 0.001*np.array(nifl_helper.ArcArray(np.array(xys)))+2 # trim the 2km that was removed

### Plot all on one set of axes
smb_offset=9000 # how many 'meters' above topo to plot this line
runoff_offset=7000 
terminus_offset=5000
velocity_offset=2500 
scaling=1000  # vertical multiplicative factor to display xcorr on same axes as topo
qual_colors = cm.get_cmap('tab20b')

fig, ax = plt.subplots(1, constrained_layout=True)

# ax.axvline(10, color='k', lw=1.0, ls='--')
ax.axvline(14, color='k', lw=0.5, ls='--')
# ax.plot(10*np.ones_like(s), np.linspace(-1300, smb_offset, len(s)), color='k', lw=1.0, ls='--')
ax.plot(s, smb_offset*np.ones_like(s), color='k', lw=0.5, ls=':')
ax.plot(s, runoff_offset*np.ones_like(s), color='k', lw=0.5, ls=':')
ax.plot(s, terminus_offset*np.ones_like(s), color='k', lw=0.5, ls=':')
ax.plot(s, smb_offset+ scaling*np.array(smb_corr_amax), color=qual_colors(0), lw=2.5, label='SMB')
ax.plot(s, runoff_offset+ scaling*np.array(runoff_corr_amax), color=qual_colors(2), lw=2.5, label='Runoff')
ax.plot(s, terminus_offset+ scaling*np.array(terminus_corr_amax), color=qual_colors(4), lw=2.5, label='Terminus')
# ax.plot(s, smb_offset+ scaling*np.array(smb_corr_amax), color='k', lw=2.0, label='SMB')
# ax.plot(s, runoff_offset+ scaling*np.array(runoff_corr_amax), color='k', lw=2.0, label='Runoff')
# ax.plot(s, terminus_offset+ scaling*np.array(terminus_corr_amax), color='k', lw=2.0, label='Terminus')

## add long-term
ax.plot(s, smb_offset+ scaling*np.array(smb_lt_corr_amax), color=qual_colors(1), ls=':', lw=2.5, label='SMB long-term')
ax.plot(s, runoff_offset+ scaling*np.array(rf_lt_corr_amax), color=qual_colors(3), ls=':', lw=2.5, label='Runoff long-term')
ax.plot(s, terminus_offset+ scaling*np.array(term_lt_corr_amax), color=qual_colors(5), ls=':', lw=2.5, label='Terminus long-term')

## add velocity
ax.plot(s, velocity_offset+scaling*np.array(scaled_speed - np.mean(scaled_speed)), color=qual_colors(17), label='Mean speed', lw=2.5)


ax.plot(xvals, bed_vals, color='saddlebrown', lw=2.0)
ax.plot(xvals, surface_vals, color='darkgrey', lw=2.0)
plt.fill_between(xvals, surface_vals, bed_vals, color='darkgrey', alpha=0.5)
plt.fill_between(xvals, bed_vals, y2=-1300, color='saddlebrown', alpha=0.5, hatch='/')
ax.set(xlim=(25, 0), ylim=(-1300, 10000), aspect=0.002, 
        xlabel='Upstream distance [km]',
        yticks=(-1000, 0, 1000, 
                velocity_offset+scaling*scaled_ticks[1],
                terminus_offset,  
                runoff_offset, 
                smb_offset),
        yticklabels=('-1000', '0', '1000 m a.s.l.', 
                     'Surface speed',
                     'Terminus',
                     'Runoff',
                     'SMB'))
minor_locator = FixedLocator([velocity_offset+scaling*scaled_ticks[0], velocity_offset+scaling*scaled_ticks[2],
                              terminus_offset-0.5*scaling, terminus_offset+0.5*scaling, 
                              runoff_offset-0.5*scaling, runoff_offset+0.5*scaling,
                              smb_offset-0.5*scaling, smb_offset+0.5*scaling])
ax.yaxis.set_minor_locator(minor_locator)
ax.tick_params(axis='both', which='major', length=4)
ax.tick_params(axis='y', which='minor', length=2)
ax.get_yticklabels()[3].set_color(qual_colors(17)) # match the ticks to the plotted color
ax.get_yticklabels()[4].set_color(qual_colors(4)) 
ax.get_yticklabels()[5].set_color(qual_colors(2))
ax.get_yticklabels()[6].set_color(qual_colors(0))


ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
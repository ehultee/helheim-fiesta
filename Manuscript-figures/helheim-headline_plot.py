#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Headline fig components
Created on Wed Jan  6 12:08:12 2021

@author: lizz
"""
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import pandas as pd

## set matplotlib font size defaults
SMALL_SIZE = 10
MEDIUM_SIZE = 11
BIGGER_SIZE = 14

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

# In[]:
## Import calving fronts to display
## hand-picked from MEaSUREs tiles
fronts_fpath = '/Users/lizz/Documents/GitHub/Data_unsynced/Helheim-processed/Helheim-selected_fronts.csv' 
fronts = pd.read_csv(fronts_fpath)
front09 = fronts['Front_points'][0].split()
f09 = np.array([(float(front09[i].strip('()[],')), float(front09[i+1].strip('()[],'))) 
       for i in range(0, len(front09), 2)]) # create an array of points from the list of strings

# In[]:
## Import velocity composite to show mean flow
vcomp_fpath ='/Users/lizz/Documents/GitHub/Data_unsynced/Sentinel-velocity/greenland_iv_500m_s1_20161223_20170227_v1_0.nc'
fh = Dataset(vcomp_fpath, mode='r')
xx = fh.variables['x'][:].copy() #x-coord (polar stereo (70, 45))
yy = fh.variables['y'][:].copy() #y-coord
v_raw = fh.variables['land_ice_surface_velocity_magnitude'][:].copy() #surface elevation
fh.close()

## Re-scale and restrict to Helheim region
xv_l, xv_r = 1700, 1950
yv_t, yv_b = 3640, 3950
sentinel_v = 0.36526*np.ma.asarray(v_raw)[yv_t:yv_b, xv_l:xv_r] ## convert from m/d to km/a


# In[]:
## COMPONENT PLOTS

## Flowline to plot
# xy_plot = [(xh[i], yh[i]) for i in range(0, upstream_max)]
xy_plot = np.array(xys)

clrs = plt.get_cmap('plasma')(np.array(range(len(xys)))/len(xys))
t_grid_trimmed = t_grid[np.argwhere(t_grid<2017)] ## valid range for interpolated funcs

highlight_idx = (5, 15) ## indices of points to plot for data/fit panel

rgb = ls.shade(np.asarray(b_hel), cmap=plt.get_cmap('gist_earth'), blend_mode='overlay',
                dx=np.mean(np.diff(x_hel)), dy=np.mean(np.diff(y_hel)), vert_exag=5.)
rgb3 = ls.shade(np.asarray(b_hel), cmap=plt.get_cmap('gray'), blend_mode='overlay',
                dx=np.mean(np.diff(x_hel)), dy=np.mean(np.diff(y_hel)), vert_exag=5.)

# ## Plot velocity stack and forcing series
# fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(6, 8), sharex=True)
# for i in range(len(xys)):
#     # ax1.plot(hel_stack.tdec, series[i], '.')
#     ax1.plot(t_grid, preds[i]['full'], label='Point {}'.format(i), color=clrs[i], lw=2.0)
#     # ax1.plot(t_grid, preds[i]['secular']+preds[i]['transient'], color=clrs[i], lw=2.0, alpha=0.3)
# ax1.set(ylabel='Surf. speed [km/a]',
#         yticks=(4, 6, 8), xlim=(2009,2017))
# ax2.scatter(smb_d_interp, np.array(smb), color='k', alpha=0.5) # raw SMB data
# ax2.plot(t_grid_trimmed, np.array(smb_func(t_grid_trimmed)), color='k', lw=2.0)
# ax2.plot(t_grid_trimmed, 1E9*np.array(smb_lowfreq(t_grid_trimmed)), color='k', lw=2.0, alpha=0.3)
# ax2.set(ylabel='Int. SMB [m$^3$ w.e.]')
# ax3.scatter(d_interp, np.array(runoff), color='k', alpha=0.5) # raw runoff data
# ax3.plot(t_grid_trimmed, np.array(runoff_func(t_grid_trimmed)), color='k', lw=2.0)
# ax3.plot(t_grid_trimmed, np.array(rf_lowfreq(t_grid_trimmed)), color='k', lw=2.0, alpha=0.3)
# ax3.set(ylabel='Int. runoff [m$^3$ w.e.]')
# ax3.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
# ax4.scatter(tm_d_interp, tm, color='k', alpha=0.5) # raw terminus data
# ax4.plot(t_grid_trimmed, termini_func(t_grid_trimmed), color='k', lw=2.0)
# ax4.plot(t_grid_trimmed, tf_lowfreq(t_grid_trimmed), color='k', lw=2.0, alpha=0.3)
# ax4.set(ylabel='Term. pos. [km]')
# ax4.set(xlim=(2009,2017), xlabel='Year')
# for ax in (ax1, ax2, ax3, ax4):
#     ax.grid(True, which='major', axis='x', ls=':', color='k', alpha=0.5)
# plt.tight_layout()


# ## plot splines with data from two points
# fig, ax = plt.subplots(1, figsize=(6,2))
# for j in highlight_idx:
#     pt = xys[j]
#     series = hel_stack.timeseries(xy=pt, key=data_key)
#     ax.plot(hel_stack.tdec, series, '.', color=clrs[j], alpha=0.5)
#     ax.plot(t_grid, preds[j]['full'], label='Point {}'.format(j), color=clrs[j], lw=2.0)
#     ax.plot(t_grid, preds[j]['secular']+preds[j]['transient'], color=clrs[j], lw=1.0, alpha=0.5)
# ax.set(xlabel='Year', ylabel='Surface speed [km/a]',
#     yticks=(4, 6, 8), xlim=(2009,2017))
# plt.tight_layout()
# plt.show()
    

# ## Plot map of flowline
# ls = LightSource(azdeg=225, altdeg=80)

# fig, ax = plt.subplots(1, figsize=(6,6))
# ax.imshow(rgb, origin='lower', extent=(x_hel[0], x_hel[-1], y_hel[0], y_hel[-1]))
# sc = ax.scatter(xy_plot[:,0], xy_plot[:,1], c=clrs)
# ax.plot(f09[:,0], f09[:,1], color='k')
# for j in highlight_idx:
#     ax.plot(xy_plot[j,0], xy_plot[j,1], marker='*', ms=15., mec='w', c=clrs[j], lw=0.)
# ax.set(xlim=(270000, 320000), xticks=(280000, 300000, 320000), 
#       ylim=(-2590000, -2550000), yticks=(-2590000, -2570000, -2550000), 
#         xticklabels=('280', '300', '320'), yticklabels=('-2590', '-2570', '-2550'),
#       xlabel='Easting [km]', ylabel='Northing [km]')
# im = ax.imshow(b_hel, cmap=plt.get_cmap('gist_earth'))
# im.remove()
# cbaxes = inset_axes(ax, width="33%", height="4%", loc='lower left', borderpad=1.5)
# cb = plt.colorbar(im, cax=cbaxes, orientation='horizontal')
# cb.set_ticks([-1000,0,2000])
# cb.ax.set_xticklabels(['-1000', '0', '2000 m a.s.l.'])
# cb.ax.xaxis.set_ticks_position('top')
# cb.ax.xaxis.set_label_position('top')
# plt.tight_layout()
# plt.show()
# # plt.savefig('/Users/lizz/Desktop/{}-helheim_map_with_colorbar'.format(
# #    datetime.date.today().strftime('%Y%m%d')))


# In[]:
## COMPOSITE FIGURE

## GridSpec version
fig3 = plt.figure(figsize=(10,8), constrained_layout=True)
gs = fig3.add_gridspec(nrows=10, ncols=2)
f3_ax1 = fig3.add_subplot(gs[0:5,0])
f3_ax2 = fig3.add_subplot(gs[0:2,1])
f3_ax3 = fig3.add_subplot(gs[2:4,1])
f3_ax4 = fig3.add_subplot(gs[4:6,1])
f3_ax5 = fig3.add_subplot(gs[6:8,1])
f3_ax6 = fig3.add_subplot(gs[8:10,1])
f3_ax7 = fig3.add_subplot(gs[5:10,0], sharex=f3_ax1, sharey=f3_ax1)

## map
f3_ax1.imshow(rgb, origin='lower', extent=(x_hel[0], x_hel[-1], y_hel[0], y_hel[-1]))
sc = f3_ax1.scatter(xy_plot[:,0], xy_plot[:,1], c=clrs)
f3_ax1.plot(f09[:,0], f09[:,1], color='w')
for j in highlight_idx:
    f3_ax1.plot(xy_plot[j,0], xy_plot[j,1], marker='*', ms=15., mec='w', c=clrs[j], lw=0.)
f3_ax1.set(xlim=(278000, 320000), xticks=(280000, 300000, 320000), 
      ylim=(-2590000, -2550000), yticks=(-2590000, -2570000, -2550000), 
      # xticklabels=('280', '300', '320'), xlabel='Easting [km]', 
      yticklabels=('-2590', '-2570', '-2550'), ylabel='Northing [km]')
plt.setp(f3_ax1.get_xticklabels(), visible=False) #hide these because the shared axis will show them
im = f3_ax1.imshow(b_hel, cmap=plt.get_cmap('gist_earth'))
im.remove()
cbaxes = inset_axes(f3_ax1, width="33%", height="4%", loc='lower left', borderpad=1.5)
cb = plt.colorbar(im, cax=cbaxes, orientation='horizontal')
cb.set_ticks([-1000,0,2000])
cb.ax.set_xticklabels(['-1000', '0', '2000 m a.s.l.'])
cb.ax.xaxis.set_ticks_position('top')
cb.ax.xaxis.set_label_position('top')

## mean velocity map
## black-white hillshade topo 
f3_ax7.imshow(rgb3, origin='lower', extent=(x_hel[0], x_hel[-1], y_hel[0], y_hel[-1]))
vcomp = f3_ax7.contourf(xx[xv_l:xv_r], yy[yv_t:yv_b], sentinel_v, levels=30, alpha=0.5)
sc = f3_ax7.scatter(xy_plot[:,0], xy_plot[:,1], facecolors='none', edgecolors='k')
f3_ax7.plot(f09[:,0], f09[:,1], color='w')
for j in highlight_idx:
    f3_ax7.plot(xy_plot[j,0], xy_plot[j,1], marker='*', ms=15., mec='w', color='k', fillstyle='none', lw=0.)
f3_ax7.set(xlim=(278000, 320000), xticks=(280000, 300000, 320000), 
      ylim=(-2590000, -2550000), yticks=(-2590000, -2570000, -2550000), 
        xticklabels=('280', '300', '320'), yticklabels=('-2590', '-2570', '-2550'),
      xlabel='Easting [km]', ylabel='Northing [km]')
cbaxes = inset_axes(f3_ax7, width="33%", height="4%", loc='lower left', borderpad=1.5)
cb = plt.colorbar(vcomp, cax=cbaxes, orientation='horizontal')
cb.set_ticks([0,3,6,9])
cb.ax.set_xticklabels(['0', '3', '6', '9 km/a'], ha='left', color='w')
cb.ax.xaxis.set_ticks_position('top')
cb.ax.xaxis.set_label_position('top')


## splines with data
for j in highlight_idx:
    pt = xys[j]
    series = hel_stack.timeseries(xy=pt, key=data_key)
    f3_ax2.plot(hel_stack.tdec, series, '.', color=clrs[j], alpha=0.5)
    f3_ax2.plot(t_grid, preds[j]['full'], label='Point {}'.format(j), color=clrs[j], lw=2.0)
    f3_ax2.plot(t_grid, preds[j]['secular']+preds[j]['transient'], color=clrs[j], lw=1.0, alpha=0.5)
f3_ax2.set(ylabel='Surf. speed [km/a]', yticks=(4, 6, 8), 
           xlim=(2009,2017), xticklabels=())
## Plot velocity stack and forcing series
for i in range(len(xys)):
    # ax1.plot(hel_stack.tdec, series[i], '.')
    f3_ax3.plot(t_grid, preds[i]['full'], label='Point {}'.format(i), color=clrs[i], lw=2.0)
    # f3_ax3.plot(t_grid, preds[i]['secular']+preds[i]['transient'], color=clrs[i], lw=2.0, alpha=0.3)
f3_ax3.set(ylabel='Surf. speed [km/a]', yticks=(4, 6, 8), 
           xticklabels=())
f3_ax4.scatter(smb_d_interp, np.array(smb), color='k', alpha=0.5) # raw SMB data
f3_ax4.plot(t_grid_trimmed, np.array(smb_func(t_grid_trimmed)), color='k', lw=2.0)
f3_ax4.plot(t_grid_trimmed, 1E9*np.array(smb_lowfreq(t_grid_trimmed)), color='k', lw=2.0, alpha=0.3)
f3_ax4.set(ylabel='Int. SMB [m$^3$ w.e.]',
           xticklabels=())
f3_ax5.scatter(d_interp, np.array(runoff), color='k', alpha=0.5) # raw runoff data
f3_ax5.plot(t_grid_trimmed, np.array(runoff_func(t_grid_trimmed)), color='k', lw=2.0)
f3_ax5.plot(t_grid_trimmed, np.array(rf_lowfreq(t_grid_trimmed)), color='k', lw=2.0, alpha=0.3)
f3_ax5.set(ylabel='Int. runoff [m$^3$ w.e.]', yticks=(0, 1E9),
           xticklabels=())
f3_ax6.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
f3_ax6.scatter(tm_d_interp, tm, color='k', alpha=0.5) # raw terminus data
f3_ax6.plot(t_grid_trimmed, termini_func(t_grid_trimmed), color='k', lw=2.0)
f3_ax6.plot(t_grid_trimmed, tf_lowfreq(t_grid_trimmed), color='k', lw=2.0, alpha=0.3)
f3_ax6.set(ylabel='Term. pos. [km]',
           xlim=(2009,2017), xlabel='Year',
           xticks=np.arange(2009,2018), xticklabels=np.arange(2009,2018))
for ax in (f3_ax2, f3_ax3, f3_ax4, f3_ax5, f3_ax6):
    ax.grid(True, which='major', axis='x', ls=':', color='k', alpha=0.5)
    ax.set(xlim=(2009,2017))
# plt.tight_layout()
plt.show()
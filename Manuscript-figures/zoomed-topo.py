#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Surface plot of bump

Created on Tue Mar 23 11:45:22 2021

@author: lizz
"""

from matplotlib import cm
from matplotlib.colors import LightSource
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import numpy as np
from matplotlib.patches import Rectangle 


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

## black-white hillshade topo 
rgb2 = ls.shade(np.asarray(b_hel)[200:300, 200:300], cmap=plt.get_cmap('gray'), blend_mode='overlay',
                dx=np.mean(np.diff(x_hel)), dy=np.mean(np.diff(y_hel)), vert_exag=5.)

fig, ax = plt.subplots()
im_dx = (x_hel[201]-x_hel[200])/2.
im_dy = (y_hel[201]-y_hel[200])/2.
im_extent = [x_hel[200]-im_dx, x_hel[300]+im_dx, y_hel[300]+im_dy, y_hel[200]-im_dy]
ax.imshow(rgb2, cmap='gray', extent=im_extent)
cf = ax.contour(x_hel[200:300], y_hel[200:300], b_hel[200:300, 200:300], 
           levels=30, 
           cmap='gist_earth', vmin=np.min(b_hel), vmax=np.max(b_hel)) #base on full map cbar
ax.plot(np.asarray(xys)[:,0], np.asarray(xys)[:,1], color='k', marker='o', path_effects=[pe.Stroke(linewidth=4, foreground='w'), pe.Normal()])
## set up a colorbar that aligns with the earlier maps displayed
cbaxes = inset_axes(ax, width="33%", height="4%", loc='lower left', borderpad=1.5)
cb = plt.colorbar(im, cax=cbaxes, orientation='horizontal')
cb.set_ticks([-1000, 0, 2000])
# cb.ax.set_xlim([-1100, 700])
cb.ax.set_xticklabels(['-1000', '0', '2000 m a.s.l.'], fontsize=11)
cb.ax.xaxis.set_ticks_position('top')
cb.ax.xaxis.set_label_position('top')
ax.set(xlim=(292000, 307000), xticks=(292000, 296000, 300000, 304000), 
       ylim=(-2582500, -2568000), yticks=(-2580000, -2570000), 
      aspect=1.,
        xticklabels=('292', '296', '300', '304'), yticklabels=('-2580', '-2570',),
      xlabel='Easting [km]', ylabel='Northing [km]')
plt.show()

## Make inset to show location
fig1, ax1 = plt.subplots()
ax1.imshow(rgb2, origin='lower', extent=(x_hel[0], x_hel[-1], y_hel[0], y_hel[-1]))
ax1.plot(np.asarray(xys)[:,0], np.asarray(xys)[:,1], color='k', marker='o', path_effects=[pe.Stroke(linewidth=2, foreground='w'), pe.Normal()])
rect=Rectangle(xy=(292000, -2582500), width=15000, height=14500, fill=False, 
                edgecolor='b', lw=4.0)
ax1.add_patch(rect)
ax1.set(xlim=(278000, 320000), xticks=(280000, 300000, 320000), 
      ylim=(-2590000, -2550000), yticks=(-2590000, -2570000, -2550000), 
       xticklabels=('280', '300', '320'), yticklabels=('-2590', '-2570', '-2550'))
plt.show()
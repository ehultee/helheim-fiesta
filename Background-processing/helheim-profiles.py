#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  7 14:26:34 2020
Flowline profiles for Helheim analysis

@author: lizz
"""
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from scipy import interpolate
import sys
sys.path.insert(0, '/Users/lizz/Documents/GitHub/nifl')
import nifl_helper


## Read in flowline for analysis
flowline_fpath = '/Users/lizz/Documents/GitHub/Data_unsynced/Felikson-flowlines/netcdfs/glaciera199.nc'
ncfile = Dataset(flowline_fpath, 'r')
xh = ncfile['flowline05'].variables['x'][:]
yh = ncfile['flowline05'].variables['y'][:]
ncfile.close()


## Read in and interpolate BedMachine topography
gl_bed_path ='/Users/lizz/Documents/GitHub/Data_unsynced/BedMachine-Greenland/BedMachineGreenland-2017-09-20.nc'
fh = Dataset(gl_bed_path, mode='r')
xx = fh.variables['x'][:].copy() #x-coord (polar stereo (70, 45))
yy = fh.variables['y'][:].copy() #y-coord
s_raw = fh.variables['surface'][:].copy() #surface elevation
h_raw=fh.variables['thickness'][:].copy() # Gridded thickness
b_raw = fh.variables['bed'][:].copy() # bed topo
thick_mask = fh.variables['mask'][:].copy()
ss = np.ma.masked_where(thick_mask !=2, s_raw)#mask values: 0=ocean, 1=ice-free land, 2=grounded ice, 3=floating ice, 4=non-Greenland land
hh = np.ma.masked_where(thick_mask !=2, h_raw) 
bb = b_raw #don't mask, to allow bed sampling from modern bathymetry (was subglacial in ~2006)
fh.close()

## Interpolate in area of Helheim
xl, xr = 6100, 6600
yt, yb = 12700, 13100
x_hel = xx[xl:xr]
y_hel = yy[yt:yb]
s_hel = ss[yt:yb, xl:xr]
b_hel = bb[yt:yb, xl:xr]
S_helheim = interpolate.RectBivariateSpline(x_hel, y_hel[::-1], s_hel.T[::,::-1]) #interpolating surface elevation provided
B_helheim = interpolate.RectBivariateSpline(x_hel, y_hel[::-1], b_hel.T[::,::-1]) #interpolating surface elevation provided

## Extract along flowlines
xyvals = np.array([(xh[i], yh[i]) for i in range(len(xh))])
bed_vals = [float(B_helheim(xh[i], yh[i])) for i in range(len(xh))]
surface_vals = [float(S_helheim(xh[i], yh[i])) for i in range(len(xh))]
xvals = (0.001*np.array(nifl_helper.ArcArray(xyvals)))-2 # align with flowline that has lower 2km trimmed

## Plot
fig, ax = plt.subplots(1, figsize=(14, 2))
ax.plot(xvals, bed_vals, color='saddlebrown')
ax.plot(xvals, surface_vals, color='darkgrey')
plt.fill_between(xvals, surface_vals, bed_vals, color='darkgrey', alpha=0.5)
plt.fill_between(xvals, bed_vals, y2=-1300, color='saddlebrown', alpha=0.5, hatch='/')
ax.set(xlim=(30, 0), ylim=(-1300, 1500), aspect=0.001, xlabel='Upstream distance [km]', ylabel='Elevation [m a.s.l.]')
plt.tight_layout()
plt.show()


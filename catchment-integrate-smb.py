#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Catchment-integrate SMB over Mankoff catchment
Uses packages in nifl-process environment -- shapely and pyshp (shapefile)

Created on Wed Nov 11 17:33:09 2020

@author: lizz
"""
from shapely.geometry import MultiPoint
from shapely.ops import triangulate
import shapefile
from netCDF4 import Dataset
import numpy as np
import pandas as pd
import pyproj as pyproj
from scipy import interpolate
import csv


## Read in Mankoff catchment from shapefile
catchment_fn = '/Users/lizz/Documents/GitHub/Data_unsynced/Helheim-processed/catchment/helheim_ice_catchment_mankoff'
sf = shapefile.Reader(catchment_fn) 
outline = sf.shapes()[0]

## Aid reprojection with BedMachine background
gl_bed_path ='/Users/lizz/Documents/GitHub/Data_unsynced/BedMachine-Greenland/BedMachineGreenland-2017-09-20.nc'
fh = Dataset(gl_bed_path, mode='r')
xx = fh.variables['x'][:].copy() #x-coord (polar stereo (70, 45))
yy = fh.variables['y'][:].copy() #y-coord
s_raw = fh.variables['surface'][:].copy() #surface elevation
fh.close()

## Read in HIRHAM SMB
gl_smb_path = '/Users/lizz/Documents/GitHub/Data_unsynced/HIRHAM5-SMB/DMI-HIRHAM5_GL2_ERAI_1980_2016_SMB_MM.nc'
fh2 = Dataset(gl_smb_path, mode='r')
x_lon = fh2.variables['lon'][:].copy() #x-coord (latlon)
y_lat = fh2.variables['lat'][:].copy() #y-coord (latlon)
ts = fh2.variables['time'][:].copy()
smb_raw = fh2.variables['smb'][:].copy()
fh2.close()

## Select topo in area of Helheim
xl1, xr1 = 5000, 7000
yt1, yb1 = 10000, 14000
x_hel = xx[xl1:xr1:20] # sample at ~3 km resolution
y_hel = yy[yt1:yb1:20]

## Select large area of SMB around Helheim
xl2, xr2 = 150, 300
yt2, yb2 = 305, 410
x_lon_h = x_lon[yt2:yb2, xl2:xr2] 
y_lat_h = y_lat[yt2:yb2, xl2:xr2] # resolution ~0.015 deg or abut 2.5 km

wgs84 = pyproj.Proj("+init=EPSG:4326") # LatLon with WGS84 datum used by HIRHAM
psn_gl = pyproj.Proj("+init=epsg:3413") # Polar Stereographic North used by BedMachine (as stated in NetCDF header)
xs, ys = pyproj.transform(wgs84, psn_gl, x_lon_h, y_lat_h)
Xmat, Ymat = np.meshgrid(x_hel, y_hel) # BedMachine coords from helheim-profiles

SMB_dict = {} #set up a dictionary of surface mass balance fields indexed by year
time_indices = range(311, 444) # go from Jan 2006 to Dec 2016 in monthly series
smb_dates = pd.date_range(start='2006-01-01', end='2016-12-31', periods=len(time_indices))
# smb_d = [d.utctimetuple() for d in smb_dates]
# dates_interp = [ice.timeutils.datestr2tdec(d[0], d[1], d[2]) for d in smb_d]
for t,d in zip(time_indices, smb_dates):
    smb_t = smb_raw[t][0][::-1, ::][yt2:yb2, xl2:xr2]
    regridded_smb_t = interpolate.griddata((xs.ravel(), ys.ravel()), smb_t.ravel(), (Xmat, Ymat), method='nearest')
    SMB_dict[d] = regridded_smb_t

## Perform Delaunay triangulation over catchment region
catchment = MultiPoint(outline.points)
triangles = triangulate(catchment)

## Sample SMB at each triangle and sum, for each time step
catchment_sum = np.zeros(len(time_indices))
for tri in triangles:
    rep_x, rep_y = tri.representative_point().x, tri.representative_point().y
    area_m2 = tri.area
    smb_x = (np.abs(x_hel - rep_x)).argmin()
    smb_y = (np.abs(y_hel - rep_y)).argmin()
    local_series = [SMB_dict[d][smb_y, smb_x]*area_m2 for d in smb_dates]
    catchment_sum = catchment_sum + np.asarray(local_series)
    
## Write out results
fn_out = '/Users/lizz/Documents/GitHub/Data_unsynced/Helheim-processed/HIRHAM_integrated_SMB.csv'
with open(fn_out, 'a', newline='') as csvfile:
    fieldnames = ['Date', 'SMB_int']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    
    for d, s in zip(smb_dates, catchment_sum):
        writer.writerow({'Date': d, 'SMB_int': s})
    

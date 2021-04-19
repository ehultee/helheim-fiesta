#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 15:39:25 2020
Test load-in of Felikson flowlines

@author: lizz
"""

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import iceutils as ice
import nifl_helper as nifl

## Pull single flowline from Denis Felikson
fp1 = '/Users/lizz/Documents/GitHub/Data_unsynced/Felikson-flowlines/netcdfs/glaciera199.nc'
ncfile = Dataset(fp1, 'r')
xh = ncfile['flowline05'].variables['x'][:]
yh = ncfile['flowline05'].variables['y'][:]
s = ncfile['flowline05']['geometry']['surface']['GIMP']['nominal'].variables['h'][:] # GIMP DEM
b = ncfile['flowline05']['geometry']['bed']['BedMachine']['nominal'].variables['h'][:] # BedMachine v3
dh = ncfile['flowline05']['dh']['GIMP-Arctic']['nominal']['dh'][:]
# ncfile.close()

fp2 = '/Users/lizz/Documents/GitHub/Data_unsynced/Felikson-flowlines/netcdfs/glacierb199.nc'
ncfile2 = Dataset(fp2, 'r')
xh2 = ncfile2['flowline05'].variables['x'][:]
yh2 = ncfile2['flowline05'].variables['y'][:]
s2 = ncfile2['flowline05']['geometry']['surface']['GIMP']['nominal'].variables['h'][:] # GIMP DEM
b2 = ncfile2['flowline05']['geometry']['bed']['BedMachine']['nominal'].variables['h'][:] # BedMachine v3
dh2 = ncfile2['flowline05']['dh']['GIMP-Arctic']['nominal']['dh'][:]
# ncfile2.close()

## Set up combined hdf5 stack
fpath='/Users/lizz/Documents/GitHub/Data_unsynced/Gld-Stack/'
hel_stack = ice.MagStack(files=[fpath+'vx.h5', fpath+'vy.h5'])
data_key = 'igram' # B. Riel convention for access to datasets in hdf5 stack

## Extract time series at selected points
upstream_max = 500 # index of last xh,yh within given distance of terminus--pts roughly 50m apart
xys = [(xh[i], yh[i]) for i in range(0, upstream_max, 25)]

# Create an evenly spaced time array for time series predictions
t_grid = np.linspace(hel_stack.tdec[0], hel_stack.tdec[-1], 1000)

# First convert the time vectors to a list of datetime
dates = ice.tdec2datestr(hel_stack.tdec, returndate=True)
dates_grid = ice.tdec2datestr(t_grid, returndate=True)

# Build the collection
collection = nifl.build_collection(dates)

# Construct a priori covariance
Cm = nifl.computeCm(collection)
iCm = np.linalg.inv(Cm)

# Instantiate a model for inversion
model = ice.tseries.Model(dates, collection=collection)

# Instantiate a model for prediction
model_pred = ice.tseries.Model(dates_grid, collection=collection)

## Access the design matrix for plotting
G = model.G

# Create lasso regression solver that does the following:
# i) Uses an a priori covariance matrix for damping out the B-splines
# ii) Uses sparsity-enforcing regularization (lasso) on the integrated B-splines
solver = ice.tseries.select_solver('lasso', reg_indices=model.itransient, penalty=0.05,
                                   rw_iter=1, regMat=iCm)

corr_max = []
lag_max = []
preds = []
for xy in xys:
    pred, st, lt = nifl.VSeriesAtPoint(xy, vel_stack=hel_stack, collection=collection, 
                                  model=model, model_pred=model_pred, solver=solver, 
                                  t_grid=t_grid, sigma=1.5, data_key='igram')
    # corr, lags, ci = nifl.SmbXcorr(xy, smb_dictionary=SMB_dict, smb_dates=smb_dates, 
    #                            velocity_pred=pred, t_grid=t_grid, diff=1)
    # corr, lags, ci = nifl.RunoffXcorr(xy, runoff_func=runoff_func, runoff_dates=d_interp, 
    #                           velocity_pred=pred, t_grid=t_grid, diff=1)
    # corr_max.append(max(corr))
    # lag_max.append(lags[np.argmax(corr)])
    preds.append(pred)

# fig, ax = plt.subplots(1)
# ax.contourf(x_hel, y_hel, b_hel, cmap='gist_earth', alpha=0.5)
# sc = ax.scatter(np.asarray(xys)[:,0], np.asarray(xys)[:,1], c=corr_max)
# fig.colorbar(sc, ax=ax)
# plt.show()

# fig1, ax1 = plt.subplots(1)
# ax1.contourf(x_hel, y_hel, b_hel, cmap='gist_earth', alpha=0.5)
# sc1 = ax1.scatter(np.asarray(xys)[:,0], np.asarray(xys)[:,1], c=lag_max, cmap='plasma')
# fig1.colorbar(sc1, ax=ax1)
# plt.show()
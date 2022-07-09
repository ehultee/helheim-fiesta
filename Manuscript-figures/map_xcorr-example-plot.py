#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pub-quality figures: maps of xcorr and lag for different variables
See Jupyter notebook for full guided analysis

Created on Fri Dec 18 13:59:51 2020
Updated Mon Feb 22 2021

@author: lizz
"""

from netCDF4 import Dataset
from scipy import interpolate
import pyproj as pyproj
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource
from mpl_toolkits.axes_grid1 import make_axes_locatable 
import iceutils as ice
import nifl_helper as nifl


# ### Define where the necessary data lives

# In[ ]:


flowline_fpath = '/Users/lizz/Documents/GitHub/Data_unsynced/Felikson-flowlines/netcdfs/glaciera199.nc'
velocity_fpath='/Volumes/Backup Plus/Large data moved 20220331/Gld-Stack/'
gl_bed_fpath ='/Users/lizz/Documents/GitHub/Data_unsynced/BedMachine-Greenland/BedMachineGreenland-2017-09-20.nc'
catchment_smb_fpath = '/Users/lizz/Documents/GitHub/Data_unsynced/Helheim-processed/smb_rec._.BN_RACMO2.3p2_ERA5_3h_FGRN055.1km.MM.csv'
runoff_fpath = '/Users/lizz/Documents/GitHub/Data_unsynced/Helheim-processed/runoff._.BN_RACMO2.3p2_ERA5_3h_FGRN055.1km.MM.csv'
termini_fpath = '/Users/lizz/Documents/GitHub/Data_unsynced/Helheim-processed/HLM_terminus_widthAVE.csv'

# ### Define the domain of analysis

# We will analyse along flowlines defined by Denis Felikson in his previous work, saved and shared as NetCDF files.  The flowlines are numbered 01-10 across the terminus; flowline 05 is close to the middle.  Note that Helheim Glacier has two large branches.  For now we'll study the main trunk, `glaciera199.nc`.  The more southerly trunk is `glacierb199.nc`.

# In[ ]:


ncfile = Dataset(flowline_fpath, 'r')
xh = ncfile['flowline05'].variables['x'][:]
yh = ncfile['flowline05'].variables['y'][:]
ncfile.close()


# In[ ]:


## Define points at which to extract
upstream_max = 500 # index of last xh,yh within given distance of terminus--pts roughly 50m apart
xys = [(xh[i], yh[i]) for i in range(0, upstream_max, 20)][2::]

# ## Import and invert velocity observations

# In[ ]:

## Set up combined hdf5 stack
hel_stack = ice.MagStack(files=[velocity_fpath+'vx.h5', velocity_fpath+'vy.h5'])
data_key = 'igram' # B. Riel convention for access to datasets in hdf5 stack


# In[ ]:


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


# Now that we are set up with our data and machinery, we'll ask the inversion to make us a continuous time series of velocity at each point we wish to study.



# In[ ]:


preds = []
for j, xy in enumerate(xys):
    try:
        pred, st, lt = nifl.VSeriesAtPoint(xy, vel_stack=hel_stack, collection=collection, 
                                  model=model, model_pred=model_pred, solver=solver, 
                                  t_grid=t_grid, sigma=2.5, data_key='igram')
        preds.append(pred)
    except AssertionError: # catches failed inversion
        print('Insufficient data for point {}. Removing'.format(j))
        xys.remove(xy)
        continue

    
# ## Comparison data sets


# ### Catchment-integrated SMB

# We load in a 1D timeseries of surface mass balance integrated over the whole Helheim catchment.  This data is monthly surface mass balance from the HIRHAM5 model, integrated over the Helheim catchment defined by K. Mankoff, with processing steps (coordinate reprojection, Delaunay triangulation, nearest-neighbor search and area summing) in `catchment-integrate-smb.py`.

# In[ ]:


## Read in RACMO monthly int from Denis
smb_racmo = pd.read_csv(catchment_smb_fpath, index_col=0, parse_dates=True)
smb_tr = smb_racmo.loc[smb_racmo.index.year >= 2006]
smb = smb_tr.loc[smb_tr.index.year <2018].squeeze()

smb_d = [d.utctimetuple() for d in smb.index]
smb_d_interp = [ice.timeutils.datestr2tdec(d[0], d[1], d[2]) for d in smb_d]
smb_func = interpolate.interp1d(smb_d_interp, smb)

# Now, we compute the normalized cross-correlation between catchment-integrated SMB and surface velocity at each point along the flowline.  We will draw on the inverted velocity series saved in `preds` above.  We save the value of the maximum normalized cross-correlation, and the value in days of the lag where it occurs, to compare with other variables later.

# In[ ]:

smb_corr_amax = []
smb_lag_amax = []
smb_significance = []

for xy, pred in zip(xys, preds):
    corr, lags, ci = nifl.Xcorr1D(xy, series_func=smb_func, series_dates=smb_d_interp, 
                              velocity_pred=pred, t_grid=t_grid, t_limits=(2009,2017), 
                              diff=1, normalize=True, pos_only=True)
    smb_corr_amax.append(corr[abs(corr).argmax()])
    smb_lag_amax.append(lags[abs(corr).argmax()])
    smb_significance.append(abs(corr[abs(corr).argmax()]) > ci[abs(corr).argmax()])


# ### Runoff

# We import monthly runoff from the RACMO model, integrated over the Helheim catchment and shared as a CSV by Denis Felikson.  Because this data is catchment-integrated, we interpolate a single 1D time series that will be used at all points.

# In[ ]:

## Read in RACMO monthly int from Denis
runoff_racmo = pd.read_csv(runoff_fpath, index_col=0, parse_dates=True)
runoff_tr = runoff_racmo.loc[runoff_racmo.index.year >= 2006]
runoff = runoff_tr.loc[runoff_tr.index.year <2018].squeeze()

runoff_d = [d.utctimetuple() for d in runoff.index]
d_interp = [ice.timeutils.datestr2tdec(d[0], d[1], d[2]) for d in runoff_d]
runoff_func = interpolate.interp1d(d_interp, runoff)


# We compute the normalized cross-correlation between catchment-integrated runoff and surface velocity at each same point.  Again we save the value of the maximum normalized cross-correlation, and the value in days of the lag where it occurs, to compare with other variables.

# In[ ]:


runoff_corr_amax = []
runoff_lag_amax = []
runoff_significance = []

for xy, pred in zip(xys, preds):
    corr, lags, ci = nifl.Xcorr1D(xy, series_func=runoff_func, series_dates=d_interp, 
                              velocity_pred=pred, t_grid=t_grid, t_limits=(2009,2017), 
                              diff=1, normalize=True, pos_only=True)
    runoff_corr_amax.append(corr[abs(corr).argmax()])
    runoff_lag_amax.append(lags[abs(corr).argmax()])
    runoff_significance.append(abs(corr[abs(corr).argmax()]) > ci[abs(corr).argmax()])


# ### Terminus position change

# We import width-averaged terminus position change processed by Leigh Stearns.  These data give terminus position in km from a baseline, so they do not need to be processed into a coordinate system.

# In[ ]:

termini = pd.read_csv(termini_fpath, index_col=0, parse_dates=True, usecols=[0,1])
trmn = termini.loc[termini.index.year >= 2006]
tm = trmn.loc[trmn.index.year <2017].squeeze()

## smooth a little to make more comparable with SMB and runoff
td = tm.rolling('10D').mean() # approximately 3 measurements per window

termini_d = [d.utctimetuple() for d in td.index]
tm_d_interp = [ice.timeutils.datestr2tdec(d[0], d[1], d[2]) for d in termini_d]
termini_func = interpolate.interp1d(tm_d_interp, td)


# In[ ]:


terminus_corr_amax = []
terminus_lag_amax = []
terminus_significance = []

for xy, pred in zip(xys, preds):
    corr, lags, ci = nifl.Xcorr1D(xy, series_func=termini_func, series_dates=tm_d_interp, 
                              velocity_pred=pred, t_grid=t_grid, t_limits=(2009,2017), 
                              diff=1, normalize=True, pos_only=True)
    terminus_corr_amax.append(corr[abs(corr).argmax()])
    terminus_lag_amax.append(lags[abs(corr).argmax()])
    terminus_significance.append(abs(corr[abs(corr).argmax()]) > ci[abs(corr).argmax()])
    

# ### Bed topography

# Mostly we will use this for plotting and for defining a standard coordinate system.  However, future analyses could combine bed topography with calving position or other variables to analyse effect on surface velocity.

# In[ ]:


## Read in and interpolate BedMachine topography
fh = Dataset(gl_bed_fpath, mode='r')
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


# In[ ]:


## Interpolate in area of Helheim
xl, xr = 6100, 6600
yt, yb = 12700, 13100
x_hel = xx[xl:xr]
y_hel = yy[yt:yb]
s_hel = ss[yt:yb, xl:xr]
b_hel = bb[yt:yb, xl:xr]
S_helheim = interpolate.RectBivariateSpline(x_hel, y_hel[::-1], s_hel.T[::,::-1]) #interpolating surface elevation provided
B_helheim = interpolate.RectBivariateSpline(x_hel, y_hel[::-1], b_hel.T[::,::-1]) #interpolating surface elevation provided

# ## Plotting




# First, we plot the max correlation at each point for a single variable.

# In[ ]:
ls = LightSource(azdeg=225, altdeg=80)

fig, ax = plt.subplots(1)
# ax.contourf(x_hel, y_hel, b_hel, cmap='gist_earth', alpha=0.5)
rgb = ls.shade(np.asarray(b_hel), cmap=plt.get_cmap('gist_earth'), blend_mode='overlay',
               dx=np.mean(np.diff(x_hel)), dy=np.mean(np.diff(y_hel)), vert_exag=5.)
ax.imshow(rgb, origin='lower', extent=(x_hel[0], x_hel[-1], y_hel[0], y_hel[-1]))
sc = ax.scatter(np.asarray(xys)[:,0], np.asarray(xys)[:,1], c=terminus_corr_amax, cmap='RdBu', vmin=-0.5, vmax=0.5)
cb = fig.colorbar(sc, ax=ax)
cb.ax.set_title('Max. xcorr')
ax.set(xlim=(270000, 320000), xticks=(280000, 300000, 320000), 
      ylim=(-2590000, -2550000), yticks=(-2590000, -2570000, -2550000), 
       xticklabels=('280', '300', '320'), yticklabels=('-2590', '-2570', '-2550'),
      xlabel='Easting [km]', ylabel='Northing [km]')
plt.show()


# Now, let's compare the patterns of correlation and lag for each variable.

# In[ ]:


div_colors = 'RdBu' # choose divergent colormap
corrnorm_min, corrnorm_max = -0.3, 0.3

fig, (ax1, ax2, ax3) = plt.subplots(nrows=1,ncols=3, figsize=(14, 4))
# ax1.contourf(x_hel, y_hel, b_hel, cmap='gist_earth', alpha=0.5)
ax1.imshow(rgb, origin='lower', extent=(x_hel[0], x_hel[-1], y_hel[0], y_hel[-1]))
sc1 = ax1.scatter(np.asarray(xys)[:,0], np.asarray(xys)[:,1], c=smb_corr_amax, cmap=div_colors,
                 vmin=corrnorm_min, vmax=corrnorm_max)
## set up correctly scaled colorbar
div1 = make_axes_locatable(ax1)
cax1 = div1.append_axes("right", size="5%", pad=0.1)
plt.colorbar(sc1, cax=cax1)
# cb1.ax.set_title('AMax. xcorr')
ax1.set(xlim=(270000, 320000), xticks=(280000, 300000, 320000), 
      ylim=(-2590000, -2550000), yticks=(-2590000, -2570000, -2550000), 
       xticklabels=('280', '300', '320'), yticklabels=('-2590', '-2570', '-2550'),
      xlabel='Easting [km]', ylabel='Northing [km]', title='Catchment SMB')

# ax2.contourf(x_hel, y_hel, b_hel, cmap='gist_earth', alpha=0.5)
ax2.imshow(rgb, origin='lower', extent=(x_hel[0], x_hel[-1], y_hel[0], y_hel[-1]))
sc2 = ax2.scatter(np.asarray(xys)[:,0], np.asarray(xys)[:,1], c=runoff_corr_amax, cmap=div_colors,
                 vmin=corrnorm_min, vmax=corrnorm_max)
## set up correctly scaled colorbar
div2 = make_axes_locatable(ax2)
cax2 = div2.append_axes("right", size="5%", pad=0.1)
fig.colorbar(sc2, cax=cax2)
# cb2.ax.set_title('AMax. xcorr')
ax2.set(xlim=(270000, 320000), xticks=(280000, 300000, 320000), 
      ylim=(-2590000, -2550000), yticks=(-2590000, -2570000, -2550000), 
       xticklabels=('280', '300', '320'), yticklabels=('-2590', '-2570', '-2550'),
      xlabel='Easting [km]', ylabel='Northing [km]', title='Catchment-integrated runoff')

# ax3.contourf(x_hel, y_hel, b_hel, cmap='gist_earth', alpha=0.5)
ax3.imshow(rgb, origin='lower', extent=(x_hel[0], x_hel[-1], y_hel[0], y_hel[-1]))
sc3 = ax3.scatter(np.asarray(xys)[:,0], np.asarray(xys)[:,1], c=terminus_corr_amax, cmap=div_colors,
                 vmin=corrnorm_min, vmax=corrnorm_max)
## set up correctly scaled colorbar
div3 = make_axes_locatable(ax3)
cax3 = div3.append_axes("right", size="5%", pad=0.1)
cb3 = fig.colorbar(sc3, cax=cax3)
cb3.ax.set_title('AMax. xcorr')
ax3.set(xlim=(270000, 320000), xticks=(280000, 300000, 320000), 
      ylim=(-2590000, -2550000), yticks=(-2590000, -2570000, -2550000), 
       xticklabels=('280', '300', '320'), yticklabels=('-2590', '-2570', '-2550'),
      xlabel='Easting [km]', ylabel='Northing [km]', title='Terminus position')
plt.tight_layout()
# plt.savefig('/Users/lizz/Desktop/20210105-map_xcorr_amax.png')

# In[ ]:
## Plot spatial pattern of lag at the absolute max xcorr

div_colors = 'RdBu' # choose divergent colormap
lagnorm_min, lagnorm_max = -365, 365

fig, (ax1, ax2, ax3) = plt.subplots(nrows=1,ncols=3, figsize=(14, 4))
# ax1.contourf(x_hel, y_hel, b_hel, cmap='gist_earth', alpha=0.5)
ax1.imshow(rgb, origin='lower', extent=(x_hel[0], x_hel[-1], y_hel[0], y_hel[-1]))
sc1 = ax1.scatter(np.asarray(xys)[:,0], np.asarray(xys)[:,1], c=smb_lag_amax, cmap=div_colors,
                  vmin=lagnorm_min, vmax=lagnorm_max)
## set up correctly scaled colorbar
div1 = make_axes_locatable(ax1)
cax1 = div1.append_axes("right", size="5%", pad=0.1)
plt.colorbar(sc1, cax=cax1)
# cb1.ax.set_title('Lag [d] at peak xcorr')
ax1.set(xlim=(270000, 320000), xticks=(280000, 300000, 320000), 
      ylim=(-2590000, -2550000), yticks=(-2590000, -2570000, -2550000), 
       xticklabels=('280', '300', '320'), yticklabels=('-2590', '-2570', '-2550'),
      xlabel='Easting [km]', ylabel='Northing [km]', title='Catchment SMB')

# ax2.contourf(x_hel, y_hel, b_hel, cmap='gist_earth', alpha=0.5)
ax2.imshow(rgb, origin='lower', extent=(x_hel[0], x_hel[-1], y_hel[0], y_hel[-1]))
sc2 = ax2.scatter(np.asarray(xys)[:,0], np.asarray(xys)[:,1], c=runoff_lag_amax, cmap=div_colors,
                  vmin=lagnorm_min, vmax=lagnorm_max)
## set up correctly scaled colorbar
div2 = make_axes_locatable(ax2)
cax2 = div2.append_axes("right", size="5%", pad=0.1)
fig.colorbar(sc2, cax=cax2)
# cb2.ax.set_title('Lag [d] at peak xcorr')
ax2.set(xlim=(270000, 320000), xticks=(280000, 300000, 320000), 
      ylim=(-2590000, -2550000), yticks=(-2590000, -2570000, -2550000), 
       xticklabels=('280', '300', '320'), yticklabels=('-2590', '-2570', '-2550'),
      xlabel='Easting [km]', ylabel='Northing [km]', title='Catchment runoff')

# ax3.contourf(x_hel, y_hel, b_hel, cmap='gist_earth', alpha=0.5)
ax3.imshow(rgb, origin='lower', extent=(x_hel[0], x_hel[-1], y_hel[0], y_hel[-1]))
sc3 = ax3.scatter(np.asarray(xys)[:,0], np.asarray(xys)[:,1], c=terminus_lag_amax, cmap=div_colors,
                  vmin=lagnorm_min, vmax=lagnorm_max)
## set up correctly scaled colorbar
div3 = make_axes_locatable(ax3)
cax3 = div3.append_axes("right", size="5%", pad=0.1)
cb3 = fig.colorbar(sc3, cax=cax3)
cb3.ax.set_title('Lag [d] at peak xcorr')
ax3.set(xlim=(270000, 320000), xticks=(280000, 300000, 320000), 
      ylim=(-2590000, -2550000), yticks=(-2590000, -2570000, -2550000), 
       xticklabels=('280', '300', '320'), yticklabels=('-2590', '-2570', '-2550'),
      xlabel='Easting [km]', ylabel='Northing [km]', title='Terminus position')
plt.tight_layout()
plt.show()
# plt.savefig('/Users/lizz/Desktop/20210105-map_lag_amax.png')

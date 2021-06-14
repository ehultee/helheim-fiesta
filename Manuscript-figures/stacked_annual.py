#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Stack annual chunks along flow
Created on Mon Jun 14 17:53:46 2021

@author: lizz
"""

from netCDF4 import Dataset
from scipy import interpolate
import pyproj as pyproj
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import iceutils as ice
import nifl_helper as nifl
import datetime # for figure labelling
from matplotlib import ticker

# ### Define where the necessary data lives

# In[ ]:


flowline_fpath = '/Users/lizz/Documents/GitHub/Data_unsynced/Felikson-flowlines/netcdfs/glaciera199.nc'
velocity_fpath='/Users/lizz/Documents/GitHub/Data_unsynced/Gld-Stack/'
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
                                  t_grid=t_grid, sigma=1.5, data_key='igram')
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

# ### Runoff

# We import monthly runoff from the RACMO model, integrated over the Helheim catchment and shared as a CSV by Denis Felikson.  Because this data is catchment-integrated, we interpolate a single 1D time series that will be used at all points.

# In[ ]:


runoff_racmo = pd.read_csv(runoff_fpath, index_col=0, parse_dates=True)
runoff_tr = runoff_racmo.loc[runoff_racmo.index.year >= 2006]
runoff = runoff_tr.loc[runoff_tr.index.year <2018].squeeze()

runoff_d = [d.utctimetuple() for d in runoff.index]
d_interp = [ice.timeutils.datestr2tdec(d[0], d[1], d[2]) for d in runoff_d]
runoff_func = interpolate.interp1d(d_interp, runoff)

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

# In[ ]:
## Now we analyse each point along flowline and store in a dictionary
annual_vals = {i: {'smb': [], 'runoff': [], 'terminus': []} 
               for i in range(len(xys))}

for j in range(len(xys)):
    point_to_plot = j
        
    # In[ ]:
    # ## Annual chunks to compare changing seasonal cycle
    
    # We break signals into annual subsets and compute the cross-correlation signal for each single year of data.
    
    
    smb_annual_corrs = []
    smb_annual_lags = []
    smb_annual_ci = []
    
    date_chks = range(2009, 2018)
    for i in range(len(date_chks)-1):
        corr, lags, ci = nifl.Xcorr1D(xys[point_to_plot], series_func=smb_func, series_dates=smb_d_interp, 
                                  velocity_pred=preds[point_to_plot], t_grid=t_grid, t_limits=(date_chks[i], date_chks[i+1]),
                                      diff=1, normalize=True, pos_only=False)
        smb_annual_corrs.append(corr)
        smb_annual_lags.append(lags)
        smb_annual_ci.append(ci)
    annual_vals[j]['smb'] = smb_annual_corrs
    
    
    # In[ ]:
    
    
    rf_annual_corrs = []
    rf_annual_lags = []
    rf_annual_ci = []
    
    for i in range(len(date_chks)-1):
        corr, lags, ci = nifl.Xcorr1D(xys[point_to_plot], series_func=runoff_func, series_dates=d_interp, 
                                  velocity_pred=preds[point_to_plot], t_grid=t_grid, t_limits=(date_chks[i], date_chks[i+1]),
                                      diff=1, normalize=True, pos_only=False)
        rf_annual_corrs.append(corr)
        rf_annual_lags.append(lags)
        rf_annual_ci.append(ci)
    annual_vals[j]['runoff'] = rf_annual_corrs
    # In[ ]:
    
    
    tm_annual_corrs = []
    tm_annual_lags = []
    tm_annual_ci = []
    
    for i in range(len(date_chks)-1):
        corr, lags, ci = nifl.Xcorr1D(xys[point_to_plot], series_func=termini_func, series_dates=tm_d_interp, 
                                  velocity_pred=preds[point_to_plot], t_grid=t_grid, t_limits=(date_chks[i], date_chks[i+1]),
                                      diff=1, normalize=True)
        tm_annual_corrs.append(corr)
        tm_annual_lags.append(lags)
        tm_annual_ci.append(ci)
    annual_vals[j]['terminus'] = tm_annual_corrs




# In[ ]:

clrs = plt.get_cmap('plasma')(np.array(range(len(xys)))/len(xys))

fig, axs = plt.subplots(len(rf_annual_corrs), ncols=3, figsize=(6, 14), sharex=True, sharey=True)
for i in range(len(smb_annual_ci)):
    ax = axs[i][0]
    ax.axvline(x=0, color='k', alpha=0.5)
    ax.axhline(y=0, color='k', alpha=0.5)
    ax.plot(smb_annual_lags[i], smb_annual_ci[i], ls=':', color='k')
    ax.plot(smb_annual_lags[i], -1*np.array(smb_annual_ci[i]), ls=':', color='k')
    for n in range(len(xys)):
        corr_color = clrs[n]
        ax.plot(smb_annual_lags[i], annual_vals[n]['smb'][i], color=corr_color, alpha=0.3)
        ax.fill_between(smb_annual_lags[i], y1=annual_vals[n]['smb'][i], y2=0, 
                        where=abs(annual_vals[n]['smb'][i])>smb_annual_ci[i], color=corr_color, alpha=0.3)
    if i==0:
        ax.set(title='SMB', ylabel='xcorr')
    elif i==len(axs)-1:
        ax.set(xlabel='Lag [d]', ylabel='xcorr')
    else:
        ax.set(ylabel='xcorr')
for j in range(len(rf_annual_ci)):
    ax = axs[j][1]
    ax.axvline(x=0, color='k', alpha=0.5)
    ax.axhline(y=0, color='k', alpha=0.5)
    ax.plot(rf_annual_lags[j], rf_annual_ci[j], ls=':', color='k')
    ax.plot(rf_annual_lags[j], -1*np.array(rf_annual_ci[j]), ls=':', color='k')
    for n in range(len(xys)):
        corr_color = clrs[n]
        ax.plot(rf_annual_lags[j], annual_vals[n]['runoff'][j], color=corr_color, alpha=0.3)
        ax.fill_between(rf_annual_lags[j], y1=annual_vals[n]['runoff'][j], y2=0, 
                        where=abs(annual_vals[n]['runoff'][j])>rf_annual_ci[j], color=corr_color, alpha=0.3)
    if j==0:
        ax.set(title='Runoff')
    elif j==len(axs)-1:
        ax.set(xlabel='Lag [d]')
    else:
        continue
for k in range(len(tm_annual_ci)):
    ax = axs[k][2]
    ax.axvline(x=0, color='k', alpha=0.5)
    ax.axhline(y=0, color='k', alpha=0.5)
    ax.plot(tm_annual_lags[k], tm_annual_ci[k], ls=':', color='k')
    ax.plot(tm_annual_lags[k], -1*np.array(tm_annual_ci[k]), ls=':', color='k')
    for n in range(len(xys)):
        corr_color = clrs[n]
        ax.plot(tm_annual_lags[k], annual_vals[n]['terminus'][k], color=corr_color, alpha=0.3)
        ax.fill_between(tm_annual_lags[k], y1=annual_vals[n]['terminus'][k], y2=0, 
                        where=abs(annual_vals[n]['terminus'][k])>tm_annual_ci[k], color=corr_color, alpha=0.3)
    ax.text(150, 0.65, str(date_chks[k]), ha='center', size=10, weight=500, color='k')
    if k==0:
        ax.set(title='Terminus pos.')
    elif k==len(axs)-1:
        ax.set(xlabel='Lag [d]')
    else:
        continue
for ax in axs.ravel():
    ax.set(ylim=(-1,1), xticks=(-300, 0, 300))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(60))
plt.tight_layout()
# plt.savefig('/Users/lizz/Desktop/{}-annual_chunk-allvars.png'.format(datetime.date.today().strftime('%Y%m%d')))






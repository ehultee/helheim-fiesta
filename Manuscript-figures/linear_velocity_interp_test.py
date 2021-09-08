#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Linear interpolation of velocity data

Created on Fri Jun 25 18:19:45 2021

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
import statsmodels.api as sm


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


# In[ ]:
# Redefine XCorr1D to use linear interp function instead
def Xcorr1D_linearV(pt, series_func, series_dates, velocity_func, t_grid, t_limits, diff=1, normalize=True, pos_only=False):
    """
    Compute cross-correlation on coincident series of a 1D time series
    (e.g. catchment-integrated runoff or SMB) versus velocity at a point.

    Parameters
    ----------
    pt : tuple
        Position (x,y) at which to pull velocity series.
    series_func : interpolate.interp1d
        1D-interpolated function with values of data series over time.
    series_dates : list
        Decimal dates of data points
    velocity_series : dict
        Output of iceutils prediction.
    t_grid : ndarray
        Evenly spaced decimal times at which spline-fit velocity is sampled
    t_limits : tuple
        Start and end dates (decimal) of the time period to study
    diff : int, optional
        Number of discrete differences to apply to data. Default is 1.
        Setting diff=0 will process the input data as-is.
    normalize : bool, optional
        Whether to normalize for a cross-correlation in [-1,1]. Default is True.
        This makes the output inter-comparable with normalized output for other
        variables.  If set to False, the signal amplitude will be larger but
        the correlation values may exceed 1.
    pos_only : bool, optional
    	Whether to analyse only xcorrs with positive lag values.  Default is False.
    	This allows a bidirectional causal relationship.  For a causal relationship 
    	hypothesised to be single-directional, choose True to display only positive
    	lag values.

    Returns
    -------
    corr : array
        Cross-correlation coefficients between SMB, velocity
    lags : array
        Time lag for each correlation value
    ci : array
        Confidence intervals for evaluation

    """
    t_min = max(min(series_dates), t_limits[0])
    t_max = min(max(series_dates), t_limits[1])
    coincident_dates = np.asarray([t for t in t_grid if (t>=t_min and t<t_max)])
    coincident_series = series_func(coincident_dates) # sample at same dates as velocity series
    series_diff = np.diff(coincident_series, n=diff)
    
    tg_0 = t_grid[np.where(t_grid>=t_min)]
    tg = tg_0[np.where(t_grid[np.where(t_grid>=t_min)]<t_max)] # trim dates to match t_limits
    vel_series = velocity_func(tg)
    vel_diff = np.diff(vel_series, n=diff)
    if normalize:
        series_diff = (series_diff-np.mean(series_diff)) / (np.std(series_diff)*len(series_diff))
        vel_diff = (vel_diff-np.mean(vel_diff)) / (np.std(vel_diff))
    corr = np.correlate(series_diff, vel_diff, mode='full')
    lags = range(int(-0.5*len(corr)), int(0.5*len(corr)+1))
    ci = [1.96/np.sqrt(len(coincident_series)-abs(k)) for k in lags]
	
    ## convert lags to physical units
    lags = np.mean(np.diff(t_grid))*365.26*np.asarray(lags)
    
    if pos_only:
    	corr = corr[np.argwhere(lags>=0)].squeeze()
    	ci = np.asarray(ci)[np.argwhere(lags>=0)].squeeze()
    	lags = lags[np.argwhere(lags>=0)].squeeze()
    	
    return corr, lags, ci

# In[ ]:


    
# ## Comparison data sets

### Linear interp velocity
linear_v = []
for j, xy in enumerate(xys):
    pt = xys[j]
    series = hel_stack.timeseries(xy=pt, key=data_key)
    v_func = interpolate.interp1d(hel_stack.tdec, series, kind='linear', fill_value='extrapolate')
    linear_v.append(v_func)

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

## Read in RACMO monthly int from Denis
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


# In[]:
## Modify the CCSL for autocorrelation
a_vel = sm.tsa.stattools.acf(np.diff(linear_v[0](t_grid)))[1]
b_runoff = sm.tsa.stattools.acf(np.diff(runoff))[1]
F_runoff = np.sqrt((1+(a_vel*b_runoff))/(1-(a_vel*b_runoff)))

## modify CCSL for terminus
b_terminus = sm.tsa.stattools.acf(np.diff(td))[1]
F_terminus = np.sqrt((1+(a_vel*b_terminus))/(1-(a_vel*b_terminus)))


## modify CCSL for smb
b_smb = sm.tsa.stattools.acf(np.diff(smb))[1]
F_smb = np.sqrt((1+(a_vel*b_smb))/(1-(a_vel*b_smb)))


# Now, we compute the normalized cross-correlation between catchment-integrated SMB and surface velocity at each point along the flowline.  We will draw on the inverted velocity series saved in `preds` above.  We save the value of the maximum normalized cross-correlation, and the value in days of the lag where it occurs, to compare with other variables later.

# In[ ]:

smb_corr_amax = []
smb_lag_amax = []
smb_significance = []

for xy, pred in zip(xys, linear_v):
    corr, lags, ci = Xcorr1D_linearV(xy, series_func=smb_func, series_dates=smb_d_interp, 
                              velocity_func=pred, t_grid=t_grid, t_limits=(2009,2017), 
                              diff=1, normalize=True, pos_only=True)
    ci_mod = F_smb*np.asarray(ci)
    smb_corr_amax.append(corr[abs(corr).argmax()])
    smb_lag_amax.append(lags[abs(corr).argmax()])
    smb_significance.append(abs(corr[abs(corr).argmax()]) > ci_mod[abs(corr).argmax()])




# We compute the normalized cross-correlation between catchment-integrated runoff and surface velocity at each same point.  Again we save the value of the maximum normalized cross-correlation, and the value in days of the lag where it occurs, to compare with other variables.

# In[ ]:


runoff_corr_amax = []
runoff_lag_amax = []
runoff_significance = []

for xy, pred in zip(xys, linear_v):
    corr, lags, ci = Xcorr1D_linearV(xy, series_func=runoff_func, series_dates=d_interp, 
                              velocity_func=pred, t_grid=t_grid, t_limits=(2009,2017), 
                              diff=1, normalize=True, pos_only=True)
    ci_mod = F_runoff*np.asarray(ci)
    runoff_corr_amax.append(corr[abs(corr).argmax()])
    runoff_lag_amax.append(lags[abs(corr).argmax()])
    runoff_significance.append(abs(corr[abs(corr).argmax()]) > ci_mod[abs(corr).argmax()])



# In[ ]:


terminus_corr_amax = []
terminus_lag_amax = []
terminus_significance = []

for xy, pred in zip(xys, linear_v):
    corr, lags, ci = Xcorr1D_linearV(xy, series_func=termini_func, series_dates=tm_d_interp, 
                              velocity_func=pred, t_grid=t_grid, t_limits=(2009,2017), 
                              diff=1, normalize=True)
    ci_mod = F_terminus *np.asarray(ci)
    terminus_corr_amax.append(corr[abs(corr).argmax()])
    terminus_lag_amax.append(lags[abs(corr).argmax()])
    terminus_significance.append(abs(corr[abs(corr).argmax()]) > ci_mod[abs(corr).argmax()])
    

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




# # First, we plot the max correlation at each point for a single variable.

# # In[ ]:
# ls = LightSource(azdeg=225, altdeg=80)

# fig, ax = plt.subplots(1)
# # ax.contourf(x_hel, y_hel, b_hel, cmap='gist_earth', alpha=0.5)
# rgb = ls.shade(np.asarray(b_hel), cmap=plt.get_cmap('gist_earth'), blend_mode='overlay',
#                dx=np.mean(np.diff(x_hel)), dy=np.mean(np.diff(y_hel)), vert_exag=5.)
# ax.imshow(rgb, origin='lower', extent=(x_hel[0], x_hel[-1], y_hel[0], y_hel[-1]))
# sc = ax.scatter(np.asarray(xys)[:,0], np.asarray(xys)[:,1], c=terminus_corr_amax, cmap='RdBu', vmin=-0.5, vmax=0.5)
# cb = fig.colorbar(sc, ax=ax)
# cb.ax.set_title('Max. xcorr')
# ax.set(xlim=(270000, 320000), xticks=(280000, 300000, 320000), 
#       ylim=(-2590000, -2550000), yticks=(-2590000, -2570000, -2550000), 
#        xticklabels=('280', '300', '320'), yticklabels=('-2590', '-2570', '-2550'),
#       xlabel='Easting [km]', ylabel='Northing [km]')
# plt.show()


# # Now, let's compare the patterns of correlation and lag for each variable.

# # In[ ]:


# div_colors = 'RdBu' # choose divergent colormap
# corrnorm_min, corrnorm_max = -0.3, 0.3

# fig, (ax1, ax2, ax3) = plt.subplots(nrows=1,ncols=3, figsize=(14, 4))
# # ax1.contourf(x_hel, y_hel, b_hel, cmap='gist_earth', alpha=0.5)
# ax1.imshow(rgb, origin='lower', extent=(x_hel[0], x_hel[-1], y_hel[0], y_hel[-1]))
# sc1 = ax1.scatter(np.asarray(xys)[:,0], np.asarray(xys)[:,1], c=smb_corr_amax, cmap=div_colors,
#                  vmin=corrnorm_min, vmax=corrnorm_max)
# ## set up correctly scaled colorbar
# div1 = make_axes_locatable(ax1)
# cax1 = div1.append_axes("right", size="5%", pad=0.1)
# plt.colorbar(sc1, cax=cax1)
# # cb1.ax.set_title('AMax. xcorr')
# ax1.set(xlim=(270000, 320000), xticks=(280000, 300000, 320000), 
#       ylim=(-2590000, -2550000), yticks=(-2590000, -2570000, -2550000), 
#        xticklabels=('280', '300', '320'), yticklabels=('-2590', '-2570', '-2550'),
#       xlabel='Easting [km]', ylabel='Northing [km]', title='Catchment SMB')

# # ax2.contourf(x_hel, y_hel, b_hel, cmap='gist_earth', alpha=0.5)
# ax2.imshow(rgb, origin='lower', extent=(x_hel[0], x_hel[-1], y_hel[0], y_hel[-1]))
# sc2 = ax2.scatter(np.asarray(xys)[:,0], np.asarray(xys)[:,1], c=runoff_corr_amax, cmap=div_colors,
#                  vmin=corrnorm_min, vmax=corrnorm_max)
# ## set up correctly scaled colorbar
# div2 = make_axes_locatable(ax2)
# cax2 = div2.append_axes("right", size="5%", pad=0.1)
# fig.colorbar(sc2, cax=cax2)
# # cb2.ax.set_title('AMax. xcorr')
# ax2.set(xlim=(270000, 320000), xticks=(280000, 300000, 320000), 
#       ylim=(-2590000, -2550000), yticks=(-2590000, -2570000, -2550000), 
#        xticklabels=('280', '300', '320'), yticklabels=('-2590', '-2570', '-2550'),
#       xlabel='Easting [km]', ylabel='Northing [km]', title='Catchment-integrated runoff')

# # ax3.contourf(x_hel, y_hel, b_hel, cmap='gist_earth', alpha=0.5)
# ax3.imshow(rgb, origin='lower', extent=(x_hel[0], x_hel[-1], y_hel[0], y_hel[-1]))
# sc3 = ax3.scatter(np.asarray(xys)[:,0], np.asarray(xys)[:,1], c=terminus_corr_amax, cmap=div_colors,
#                  vmin=corrnorm_min, vmax=corrnorm_max)
# ## set up correctly scaled colorbar
# div3 = make_axes_locatable(ax3)
# cax3 = div3.append_axes("right", size="5%", pad=0.1)
# cb3 = fig.colorbar(sc3, cax=cax3)
# cb3.ax.set_title('AMax. xcorr')
# ax3.set(xlim=(270000, 320000), xticks=(280000, 300000, 320000), 
#       ylim=(-2590000, -2550000), yticks=(-2590000, -2570000, -2550000), 
#        xticklabels=('280', '300', '320'), yticklabels=('-2590', '-2570', '-2550'),
#       xlabel='Easting [km]', ylabel='Northing [km]', title='Terminus position')
# plt.tight_layout()
# # plt.savefig('/Users/lizz/Desktop/20210105-map_xcorr_amax.png')

# # In[ ]:
# ## Plot spatial pattern of lag at the absolute max xcorr

# div_colors = 'RdBu' # choose divergent colormap
# lagnorm_min, lagnorm_max = -365, 365

# fig, (ax1, ax2, ax3) = plt.subplots(nrows=1,ncols=3, figsize=(14, 4))
# # ax1.contourf(x_hel, y_hel, b_hel, cmap='gist_earth', alpha=0.5)
# ax1.imshow(rgb, origin='lower', extent=(x_hel[0], x_hel[-1], y_hel[0], y_hel[-1]))
# sc1 = ax1.scatter(np.asarray(xys)[:,0], np.asarray(xys)[:,1], c=smb_lag_amax, cmap=div_colors,
#                   vmin=lagnorm_min, vmax=lagnorm_max)
# ## set up correctly scaled colorbar
# div1 = make_axes_locatable(ax1)
# cax1 = div1.append_axes("right", size="5%", pad=0.1)
# plt.colorbar(sc1, cax=cax1)
# # cb1.ax.set_title('Lag [d] at peak xcorr')
# ax1.set(xlim=(270000, 320000), xticks=(280000, 300000, 320000), 
#       ylim=(-2590000, -2550000), yticks=(-2590000, -2570000, -2550000), 
#        xticklabels=('280', '300', '320'), yticklabels=('-2590', '-2570', '-2550'),
#       xlabel='Easting [km]', ylabel='Northing [km]', title='Catchment SMB')

# # ax2.contourf(x_hel, y_hel, b_hel, cmap='gist_earth', alpha=0.5)
# ax2.imshow(rgb, origin='lower', extent=(x_hel[0], x_hel[-1], y_hel[0], y_hel[-1]))
# sc2 = ax2.scatter(np.asarray(xys)[:,0], np.asarray(xys)[:,1], c=runoff_lag_amax, cmap=div_colors,
#                   vmin=lagnorm_min, vmax=lagnorm_max)
# ## set up correctly scaled colorbar
# div2 = make_axes_locatable(ax2)
# cax2 = div2.append_axes("right", size="5%", pad=0.1)
# fig.colorbar(sc2, cax=cax2)
# # cb2.ax.set_title('Lag [d] at peak xcorr')
# ax2.set(xlim=(270000, 320000), xticks=(280000, 300000, 320000), 
#       ylim=(-2590000, -2550000), yticks=(-2590000, -2570000, -2550000), 
#        xticklabels=('280', '300', '320'), yticklabels=('-2590', '-2570', '-2550'),
#       xlabel='Easting [km]', ylabel='Northing [km]', title='Catchment runoff')

# # ax3.contourf(x_hel, y_hel, b_hel, cmap='gist_earth', alpha=0.5)
# ax3.imshow(rgb, origin='lower', extent=(x_hel[0], x_hel[-1], y_hel[0], y_hel[-1]))
# sc3 = ax3.scatter(np.asarray(xys)[:,0], np.asarray(xys)[:,1], c=terminus_lag_amax, cmap=div_colors,
#                   vmin=lagnorm_min, vmax=lagnorm_max)
# ## set up correctly scaled colorbar
# div3 = make_axes_locatable(ax3)
# cax3 = div3.append_axes("right", size="5%", pad=0.1)
# cb3 = fig.colorbar(sc3, cax=cax3)
# cb3.ax.set_title('Lag [d] at peak xcorr')
# ax3.set(xlim=(270000, 320000), xticks=(280000, 300000, 320000), 
#       ylim=(-2590000, -2550000), yticks=(-2590000, -2570000, -2550000), 
#        xticklabels=('280', '300', '320'), yticklabels=('-2590', '-2570', '-2550'),
#       xlabel='Easting [km]', ylabel='Northing [km]', title='Terminus position')
# plt.tight_layout()
# plt.show()
# # plt.savefig('/Users/lizz/Desktop/20210105-map_lag_amax.png')

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 10:34:57 2021
Cross-correlation of isolated low frequency or high frequency variability
Stack plot of these signals

@author: lizz
"""

import numpy as np
import pandas as pd
import csv
import iceutils as ice
from scipy import interpolate
from scipy import ndimage


## Test low-frequency variability
def Xcorr1D_lt(pt, series_func, series_dates, velocity_pred, t_grid, t_limits, 
               diff=1, normalize=True, pos_only=False):
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
    
    vel_longterm = velocity_pred['secular'] + velocity_pred['transient']
    series_diff = np.diff(coincident_series, n=diff)
    vel_series_0 = vel_longterm[np.where(t_grid>=t_min)]
    vel_series = vel_series_0[np.where(t_grid[np.where(t_grid>=t_min)]<t_max)] # trim dates to match t_limits
    vel_diff = np.diff(vel_series, n=diff)
    if normalize:
        series_diff = (series_diff-np.mean(series_diff)) / (np.std(series_diff)*len(series_diff))
        vel_diff = (vel_diff-np.mean(vel_diff)) / (np.std(vel_diff))
    corr = np.correlate(series_diff, vel_diff, mode='full')
    lags = range(int(-0.5*len(corr)), int(0.5*len(corr)+1))
    ci = [2/np.sqrt(len(coincident_series)-abs(k)) for k in lags]

    ## convert lags to physical units
    lags = np.mean(np.diff(t_grid))*365.26*np.asarray(lags)
    
    if pos_only:
    	corr = corr[np.argwhere(lags>=0)].squeeze()
    	ci = np.asarray(ci)[np.argwhere(lags>=0)].squeeze()
    	lags = lags[np.argwhere(lags>=0)].squeeze()
    
    return corr, lags, ci

# In[ ]:
## Low-frequency terminus variability
t_grid_trimmed = t_grid[np.argwhere(t_grid<2017)] ## valid range for interpolated funcs

tm_evensampled = termini_func(t_grid_trimmed).squeeze()
window = np.int(2./np.mean(np.diff(t_grid_trimmed.squeeze()))) # set window size to the number of time steps in 2 years
tm_filtered = ndimage.uniform_filter1d(tm_evensampled, size=window)
tf_lowfreq = interpolate.UnivariateSpline(t_grid_trimmed, tm_filtered, s=0)


term_lt_corr_amax = []
term_lt_lag_amax = []
for xy, pred in zip(xys, preds):
    corr, lags, ci = Xcorr1D_lt(xy, series_func=tf_lowfreq, series_dates=tm_d_interp, 
                              velocity_pred=pred, t_grid=t_grid, t_limits=(2009,2017), 
                              diff=0, normalize=True)
    term_lt_corr_amax.append(corr[abs(corr).argmax()])
    term_lt_lag_amax.append(lags[abs(corr).argmax()])

# In[ ]:
## Low-frequency runoff variability
rnf_evensampled = runoff_func(t_grid_trimmed).squeeze()
rnf_filtered = ndimage.uniform_filter1d(rnf_evensampled, size=window)
rf_lowfreq = interpolate.UnivariateSpline(t_grid_trimmed, rnf_filtered, s=0)

rf_lt_corr_amax = []
rf_lt_lag_amax = []
for xy, pred in zip(xys, preds):
    corr, lags, ci = Xcorr1D_lt(xy, series_func=rf_lowfreq, series_dates=d_interp, 
                              velocity_pred=pred, t_grid=t_grid, t_limits=(2009,2017), 
                              diff=0, normalize=True, pos_only=True)
    rf_lt_corr_amax.append(corr[abs(corr).argmax()])
    rf_lt_lag_amax.append(lags[abs(corr).argmax()])

# In[ ]:
## Low-frequency SMB variability
smb_evensampled = 1E-9*np.array(smb_func(t_grid_trimmed).squeeze())
smb_filtered = ndimage.uniform_filter1d(smb_evensampled, size=window)
smb_lowfreq = interpolate.UnivariateSpline(t_grid_trimmed, smb_filtered, s=0)

smb_lt_corr_amax = []
smb_lt_lag_amax = []
for xy, pred in zip(xys, preds):
    corr, lags, ci = Xcorr1D_lt(xy, series_func=smb_lowfreq, series_dates=smb_d_interp, 
                              velocity_pred=pred, t_grid=t_grid, t_limits=(2009,2017), 
                              diff=0, normalize=True, pos_only=True)
    smb_lt_corr_amax.append(corr[abs(corr).argmax()])
    smb_lt_lag_amax.append(lags[abs(corr).argmax()])

# # In[ ]:
# ## Plot the low-frequency signals in stack
# fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(8, 8), sharex=True)
# for i in range(len(xys)):
#     # ax1.plot(hel_stack.tdec, series[i], '.')
#     ax1.plot(t_grid, preds[i]['full'], label='Point {}'.format(i), color=clrs[i], lw=1.0, alpha=0.3)
#     ax1.plot(t_grid, preds[i]['secular']+preds[i]['transient'], color=clrs[i], lw=2.0)
# ax1.set(ylabel='Surf. speed [km/a]',
#         yticks=(4, 6, 8), xlim=(2009,2017))

# ax2.scatter(smb_d_interp, 0.001*np.array(smb['SMB_int']), color='k', alpha=0.3) # raw SMB data
# ax2.plot(t_grid_trimmed, 0.001*np.array(smb_func(t_grid_trimmed)), color='k', alpha=0.3)
# ax2.plot(t_grid_trimmed, 1E9*np.array(smb_lowfreq(t_grid_trimmed)), color='k', alpha=0.7)
# ax2.set(ylabel='Int. SMB [m3 w.e.]')

# ax3.scatter(d_interp, 1000*np.array(rf[:,2]), color='k', alpha=0.3) # raw runoff data
# ax3.plot(t_grid_trimmed, 1000*np.array(runoff_func(t_grid_trimmed)), color='k', alpha=0.3)
# ax3.plot(t_grid_trimmed, 1000*np.array(rf_lowfreq(t_grid_trimmed)), color='k', alpha=0.7)
# ax3.set(ylabel='Int. runoff [m3 w.e.]')
# ax3.ticklabel_format(axis='y', style='sci', scilimits=(0,0))

# ax4.scatter(tm_d_interp, tm['term_km'], color='k', alpha=0.3) # raw terminus data
# ax4.plot(t_grid_trimmed, termini_func(t_grid_trimmed), color='k', alpha=0.3)
# ax4.plot(t_grid_trimmed, tf_lowfreq(t_grid_trimmed), color='k', alpha=0.7)
# ax4.set(ylabel='Term. pos. [km]')
# ax4.set(xlim=(2009,2017), xlabel='Year')
# for ax in (ax1, ax2, ax3, ax4):
#     ax.grid(True, which='major', axis='x', ls=':', color='k', alpha=0.5)
# plt.tight_layout()


# In[ ]:
## Set up six-panel composite
div_colors = 'RdBu' # choose divergent colormap for xcorr
lag_colors = 'PiYG' # choose divergent colormap for lag
corrnorm_min, corrnorm_max = -0.3, 0.3
lagnorm_min, lagnorm_max = -365, 365

## set matplotlib font size defaults
SMALL_SIZE = 10
MEDIUM_SIZE = 12
BIGGER_SIZE = 14

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

## black-white hillshade topo underneath
rgb2 = ls.shade(np.asarray(b_hel), cmap=plt.get_cmap('gray'), blend_mode='overlay',
                dx=np.mean(np.diff(x_hel)), dy=np.mean(np.diff(y_hel)), vert_exag=5.)

fig, ((ax1, ax2, ax3), (ax4,ax5,ax6)) = plt.subplots(nrows=2,ncols=3, figsize=(12, 8), 
                                                      # constrained_layout=True, 
                                                      sharex=True, sharey=True,
                                                      gridspec_kw={'wspace':0.01})
    
ax1.imshow(rgb2, origin='lower', extent=(x_hel[0], x_hel[-1], y_hel[0], y_hel[-1]))
sc1 = ax1.scatter(np.asarray(xys)[:,0], np.asarray(xys)[:,1], c=smb_lt_corr_amax, cmap=div_colors,
                  vmin=corrnorm_min, vmax=corrnorm_max)
# ## set up correctly scaled colorbar
# div1 = make_axes_locatable(ax1)
# cax1 = div1.append_axes("right", size="5%", pad=0.1)
# plt.colorbar(sc1, cax=cax1)
# cb1.ax.set_title('AMax. xcorr')
ax1.set(xlim=(278000, 320000), xticks=(280000, 300000, 320000), 
        ylim=(-2590000, -2550000), yticks=(-2590000, -2570000, -2550000), 
        xticklabels=('280', '300', '320'), yticklabels=('-2590', '-2570', '-2550'),
        ylabel='Northing [km]', title='Catchment SMB')

ax2.imshow(rgb2, origin='lower', extent=(x_hel[0], x_hel[-1], y_hel[0], y_hel[-1]))
sc2 = ax2.scatter(np.asarray(xys)[:,0], np.asarray(xys)[:,1], c=rf_lt_corr_amax, cmap=div_colors,
                  vmin=corrnorm_min, vmax=corrnorm_max)
# ## set up correctly scaled colorbar
# div2 = make_axes_locatable(ax2)
# cax2 = div2.append_axes("right", size="5%", pad=0.1)
# fig.colorbar(sc2, cax=cax2)
# cb2.ax.set_title('AMax. xcorr')
ax2.set(xlim=(278000, 320000), xticks=(280000, 300000, 320000), 
      ylim=(-2590000, -2550000), yticks=(-2590000, -2570000, -2550000), 
        xticklabels=('280', '300', '320'), yticklabels=('-2590', '-2570', '-2550'),
        title='Catchment runoff')

ax3.imshow(rgb2, origin='lower', extent=(x_hel[0], x_hel[-1], y_hel[0], y_hel[-1]))
sc3 = ax3.scatter(np.asarray(xys)[:,0], np.asarray(xys)[:,1], c=term_lt_corr_amax, cmap=div_colors,
                  vmin=corrnorm_min, vmax=corrnorm_max)
## set up correctly scaled colorbar - one for all xcorr plots
div3 = make_axes_locatable(ax3)
cax3 = div3.append_axes("right", size="5%", pad=0.1)
cb3 = fig.colorbar(sc3, cax=cax3, extend='both')
cb3.ax.set_ylabel('AMax. xcorr')
ax3.set(xlim=(278000, 320000), xticks=(280000, 300000, 320000), 
      ylim=(-2590000, -2550000), yticks=(-2590000, -2570000, -2550000), 
        xticklabels=('280', '300', '320'), yticklabels=('-2590', '-2570', '-2550'),
        title='Terminus position', aspect=1.)

## SECOND ROW
ax4.imshow(rgb2, origin='lower', extent=(x_hel[0], x_hel[-1], y_hel[0], y_hel[-1]))
sc4 = ax4.scatter(np.asarray(xys)[:,0], np.asarray(xys)[:,1], c=smb_lt_lag_amax, cmap=lag_colors,
                  vmin=lagnorm_min, vmax=lagnorm_max)
# ## set up correctly scaled colorbar
# div4 = make_axes_locatable(ax4)
# cax4 = div4.append_axes("right", size="5%", pad=0.1)
# plt.colorbar(sc4, cax=cax4)
# cb1.ax.set_title('Lag [d] at peak xcorr')
ax4.set(xlim=(278000, 320000), xticks=(280000, 300000, 320000), 
      ylim=(-2590000, -2550000), yticks=(-2590000, -2570000, -2550000), 
        xticklabels=('280', '300', '320'), yticklabels=('-2590', '-2570', '-2550'),
      xlabel='Easting [km]', ylabel='Northing [km]')

ax5.imshow(rgb2, origin='lower', extent=(x_hel[0], x_hel[-1], y_hel[0], y_hel[-1]))
sc5 = ax5.scatter(np.asarray(xys)[:,0], np.asarray(xys)[:,1], c=rf_lt_lag_amax, cmap=lag_colors,
                  vmin=lagnorm_min, vmax=lagnorm_max)
# ## set up correctly scaled colorbar
# div5 = make_axes_locatable(ax5)
# cax5 = div5.append_axes("right", size="5%", pad=0.1)
# fig.colorbar(sc5, cax=cax5)
# cb2.ax.set_title('Lag [d] at peak xcorr')
ax5.set(xlim=(278000, 320000), xticks=(280000, 300000, 320000), 
      ylim=(-2590000, -2550000), yticks=(-2590000, -2570000, -2550000), 
        xticklabels=('280', '300', '320'), yticklabels=('-2590', '-2570', '-2550'),
      xlabel='Easting [km]')

ax6.imshow(rgb2, origin='lower', extent=(x_hel[0], x_hel[-1], y_hel[0], y_hel[-1]))
sc6 = ax6.scatter(np.asarray(xys)[:,0], np.asarray(xys)[:,1], c=term_lt_lag_amax, cmap=lag_colors,
                  vmin=lagnorm_min, vmax=lagnorm_max)
## set up correctly scaled colorbar
div6 = make_axes_locatable(ax6)
cax6 = div6.append_axes("right", size="5%", pad=0.1)
cb6 = fig.colorbar(sc6, cax=cax6, extend='min')
cb6.ax.set_ylabel('Lag [d] at peak xcorr')
ax6.set(xlim=(278000, 320000), xticks=(280000, 300000, 320000), 
      ylim=(-2590000, -2550000), yticks=(-2590000, -2570000, -2550000), 
        xticklabels=('280', '300', '320'), yticklabels=('-2590', '-2570', '-2550'),
      xlabel='Easting [km]', aspect=1.)
# plt.tight_layout()
plt.show()
# plt.savefig('/Users/lizz/Desktop/20210204-helheim-longterm_xcorr_lag_composite')
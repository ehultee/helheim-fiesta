#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot correlations of smoothed UnivariateSpline in familiar grid 
for inclusion  in supplement?
Created on Thu Jun 24 17:43:47 2021

@author: lizz
"""

## base these on previously plotted correlograms
spline_runoff = interpolate.interp1d(d_interp, runoff, kind='cubic')
spline_terminus = interpolate.interp1d(tm_d_interp, td, kind='cubic')


corr, lags, ci = nifl.Xcorr1D(xys[0], series_func=spline_runoff, series_dates=d_interp, 
                          velocity_pred=preds[0], t_grid=t_grid, t_limits=(2009,2017), 
                          diff=1, normalize=True, pos_only=False)
ci_mod = F_runoff*np.asarray(ci)


clrs = plt.get_cmap('plasma')(np.array(range(len(xys)))/len(xys))

fig, ax = plt.subplots(1)
ax.axhline(0, color='k', alpha=0.5)
ax.axvline(0, color='k', alpha=0.5)
ax.plot(lags, corr, color=clrs[0])
ax.plot(lags, ci_mod, ls=':', color='k')
ax.plot(lags, -1*ci_mod, ls=':', color='k')
ax.fill_between(lags, y1=corr, y2=0, 
                    where=abs(corr)>ci_mod, color=clrs[0])
ax.set(title='Runoff-velocity correlogram', ylim=(-0.5, 0.5), xlim=(-2500, 2500),
       xlabel='Lag [d]', ylabel='Xcorr')



corrT, lagsT, ciT = nifl.Xcorr1D(xys[0], series_func=spline_terminus, series_dates=tm_d_interp, 
                          velocity_pred=preds[0], t_grid=t_grid, t_limits=(2009,2017), 
                          diff=1, normalize=True)
ci_modT = F_terminus *np.asarray(ci)


clrs = plt.get_cmap('plasma')(np.array(range(len(xys)))/len(xys))

fig, ax = plt.subplots(1)
ax.axhline(0, color='k', alpha=0.5)
ax.axvline(0, color='k', alpha=0.5)
ax.plot(lagsT, corrT, color=clrs[0])
ax.plot(lagsT, ci_modT, ls=':', color='k')
ax.plot(lagsT, -1*ci_modT, ls=':', color='k')
ax.fill_between(lagsT, y1=corrT, y2=0, 
                    where=abs(corrT)>ci_modT, color=clrs[0])
ax.set(title='Terminus-velocity correlogram', ylim=(-0.5, 0.5), xlim=(-2500, 2500),
       xlabel='Lag [d]', ylabel='Xcorr')

fig, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=True)
for ax in (ax1, ax2):
    ax.axhline(0, color='k', alpha=0.5)
    ax.axvline(0, color='k', alpha=0.5)
ax1.plot(lags, corr, color=clrs[0])
ax1.plot(lags, ci_mod, ls=':', color='k')
ax1.plot(lags, -1*ci_mod, ls=':', color='k')
ax1.fill_between(lags, y1=corr, y2=0, 
                    where=abs(corr)>ci_mod, color=clrs[0])
ax1.set(ylim=(-0.5, 0.5), xlim=(-2500, 2500), ylabel='Runoff Xcorr')
ax2.plot(lagsT, corrT, color=clrs[0])
ax2.plot(lagsT, ci_modT, ls=':', color='k')
ax2.plot(lagsT, -1*ci_modT, ls=':', color='k')
ax2.fill_between(lagsT, y1=corrT, y2=0, 
                    where=abs(corrT)>ci_modT, color=clrs[0])
ax2.set(ylim=(-0.5, 0.5), xlim=(-2500, 2500),
       xlabel='Lag [d]', ylabel='Terminus Xcorr')
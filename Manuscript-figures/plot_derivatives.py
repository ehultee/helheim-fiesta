#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot derivatives of signals
Created on Fri Jun 25 10:39:09 2021

@author: lizz
"""

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

for yr in (2013,):
    
    ## GridSpec version
    fig, ((ax1, ax3), (ax2, ax4)) = plt.subplots(2,2, sharex=True, figsize=(8,6))
    
    
    ## splines with data - full signal
    for j in highlight_idx:
        pt = xys[j]
        series = hel_stack.timeseries(xy=pt, key=data_key)
        ax1.plot(hel_stack.tdec, series, '.', color=clrs[j], alpha=0.5)
        ax1.plot(t_grid, preds[j]['full'], label='Point {}'.format(j), color=clrs[j], lw=2.0)
        ax3.plot(t_grid[1::], np.diff(preds[j]['full']), color=clrs[j], lw=2.0)
    ax1.set(ylabel='Surf. speed [km/a]', yticks=(4, 6, 8), 
               xlim=(yr,yr+1), xticklabels=())
    ax3.set(ylabel=r'$\Delta$ surf. speed')
    ax3.axhline(0, ls=':', color='k')
    
    
    
    # ## Plot velocity stack and forcing series
    # for i in range(len(xys)):
    #     # ax1.plot(hel_stack.tdec, series[i], '.')
    #     f3_ax3.plot(t_grid, preds[i]['full'], label='Point {}'.format(i), color=clrs[i], lw=2.0)
    #     # f3_ax3.plot(t_grid, preds[i]['secular']+preds[i]['transient'], color=clrs[i], lw=2.0, alpha=0.3)
    # f3_ax3.set(ylabel='Surf. speed [km/a]', yticks=(4, 6, 8), 
    #            xticklabels=())
    
    # f3_ax4.scatter(smb_d_interp, np.array(smb), color='k', alpha=0.5) # raw SMB data
    # f3_ax4.plot(t_grid_trimmed, np.array(smb_func(t_grid_trimmed)), color='k', lw=2.0)
    # f3_ax4.plot(t_grid_trimmed, 1E9*np.array(smb_lowfreq(t_grid_trimmed)), color='k', lw=2.0, alpha=0.3)
    # f3_ax4.set(ylabel='Int. SMB [m$^3$ w.e.]',
    #            xticklabels=())
    
    ax2.scatter(d_interp, np.array(runoff), color='k', alpha=0.5) # raw runoff data
    ax2.plot(t_grid_trimmed, np.array(runoff_func(t_grid_trimmed)), color='k', lw=2.0)
    # ax2.plot(t_grid_trimmed, np.array(rf_lowfreq(t_grid_trimmed)), color='k', lw=2.0, alpha=0.3)
    ax2.set(ylabel='Int. runoff [m$^3$ w.e.]', yticks=(0, 1E9),
               xticklabels=())
    
    ax4.plot(t_grid_trimmed[1::].ravel(), np.diff(runoff_func(t_grid_trimmed), axis=0), color='k', lw=2.0)
    # ax3.plot(t_grid_trimmed, np.array(rf_lowfreq(t_grid_trimmed)), color='k', lw=2.0, alpha=0.3)
    ax4.axhline(0, ls=':', color='k')
    ax4.set(ylabel=r'$\Delta$ runoff',
               xticks=np.arange(yr, yr+1.1, 0.25),
               xticklabels=np.arange(yr, yr+1.1, 0.25))
    
    # f3_ax6.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    # f3_ax6.scatter(tm_d_interp, tm, color='k', alpha=0.5) # raw terminus data
    # f3_ax6.plot(t_grid_trimmed, termini_func(t_grid_trimmed), color='k', lw=2.0)
    # f3_ax6.plot(t_grid_trimmed, tf_lowfreq(t_grid_trimmed), color='k', lw=2.0, alpha=0.3)
    # f3_ax6.set(ylabel='Term. pos. [km]',
    #            xlim=(2009,2017), xlabel='Year',
    #            xticks=np.arange(2009,2018), xticklabels=np.arange(2009,2018))
    for ax in (ax1, ax2, ax3, ax4):
        ax.grid(True, which='both', axis='x', ls=':', color='k', alpha=0.5)
        ax.set(xlim=(yr,yr+1))
    plt.tight_layout()
    plt.show()
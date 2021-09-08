#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Six-panel xcorr *and* lag
Created on Thu Feb  4 18:21:08 2021

@author: lizz
"""

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

sig_markers = ['o', 'x']

## black-white hillshade topo underneath
rgb2 = ls.shade(np.asarray(b_hel), cmap=plt.get_cmap('gray'), blend_mode='overlay',
               dx=np.mean(np.diff(x_hel)), dy=np.mean(np.diff(y_hel)), vert_exag=5.)

fig, ((ax1, ax2, ax3), (ax4,ax5,ax6)) = plt.subplots(nrows=2,ncols=3, figsize=(12, 8), 
                                                     # constrained_layout=True, 
                                                     sharex=True, sharey=True,
                                                     gridspec_kw={'wspace':0.01})
    
ax1.imshow(rgb2, origin='lower', extent=(x_hel[0], x_hel[-1], y_hel[0], y_hel[-1]))
sc1 = ax1.scatter(np.asarray(xys)[smb_significance,0], np.asarray(xys)[smb_significance,1], 
                  c=np.asarray(smb_corr_amax)[smb_significance], cmap=div_colors, marker=sig_markers[0], 
                  vmin=corrnorm_min, vmax=corrnorm_max)
ax1.scatter(np.asarray(xys)[np.invert(smb_significance),0], np.asarray(xys)[np.invert(smb_significance),1], 
                  c=np.asarray(smb_corr_amax)[np.invert(smb_significance)], cmap=div_colors, marker=sig_markers[1], 
                  vmin=corrnorm_min, vmax=corrnorm_max) #different marker for insig values
# sc1 = ax1.scatter(np.asarray(xys)[:,0], np.asarray(xys)[:,1], c=smb_corr_amax, cmap=div_colors,
#                  vmin=corrnorm_min, vmax=corrnorm_max)
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
sc2 = ax2.scatter(np.asarray(xys)[runoff_significance,0], np.asarray(xys)[runoff_significance,1], 
                  c=np.asarray(runoff_corr_amax)[runoff_significance], cmap=div_colors, marker=sig_markers[0],
                  vmin=corrnorm_min, vmax=corrnorm_max)
ax2.scatter(np.asarray(xys)[np.invert(runoff_significance),0], np.asarray(xys)[np.invert(runoff_significance),1],
            c=np.asarray(runoff_corr_amax)[np.invert(runoff_significance)], cmap=div_colors, marker=sig_markers[1],
            vmin=corrnorm_min, vmax=corrnorm_max) # distinguish insig values
# sc2 = ax2.scatter(np.asarray(xys)[:,0], np.asarray(xys)[:,1], c=runoff_corr_amax, cmap=div_colors,
#                  vmin=corrnorm_min, vmax=corrnorm_max)
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
sc3 = ax3.scatter(np.asarray(xys)[terminus_significance,0], np.asarray(xys)[terminus_significance,1], 
                  c=np.asarray(terminus_corr_amax)[terminus_significance], cmap=div_colors, marker=sig_markers[0],
                  vmin=corrnorm_min, vmax=corrnorm_max)
ax3.scatter(np.asarray(xys)[np.invert(terminus_significance),0], np.asarray(xys)[np.invert(terminus_significance),1],
            c=np.asarray(terminus_corr_amax)[np.invert(terminus_significance)], cmap=div_colors, marker=sig_markers[1],
            vmin=corrnorm_min, vmax=corrnorm_max)
# sc3 = ax3.scatter(np.asarray(xys)[:,0], np.asarray(xys)[:,1], c=terminus_corr_amax, cmap=div_colors,
#                  vmin=corrnorm_min, vmax=corrnorm_max)
## set up correctly scaled colorbar - one for all xcorr plots
div3 = make_axes_locatable(ax3)
cax3 = div3.append_axes("right", size="5%", pad=0.1)
cb3 = fig.colorbar(sc3, cax=cax3)
cb3.ax.set_ylabel('AMax. xcorr')
ax3.set(xlim=(278000, 320000), xticks=(280000, 300000, 320000), 
      ylim=(-2590000, -2550000), yticks=(-2590000, -2570000, -2550000), 
       xticklabels=('280', '300', '320'), yticklabels=('-2590', '-2570', '-2550'),
       title='Terminus position', aspect=1.)

## SECOND ROW
ax4.imshow(rgb2, origin='lower', extent=(x_hel[0], x_hel[-1], y_hel[0], y_hel[-1]))
sc4 = ax4.scatter(np.asarray(xys)[smb_significance,0], np.asarray(xys)[smb_significance,1], 
                  c=np.asarray(smb_lag_amax)[smb_significance], cmap=lag_colors, marker=sig_markers[0],
                  vmin=lagnorm_min, vmax=lagnorm_max)
ax4.scatter(np.asarray(xys)[np.invert(smb_significance),0], np.asarray(xys)[np.invert(smb_significance),1], 
            c=np.asarray(smb_lag_amax)[np.invert(smb_significance)], cmap=lag_colors, marker=sig_markers[1],
                  vmin=lagnorm_min, vmax=lagnorm_max)
# sc4 = ax4.scatter(np.asarray(xys)[:,0], np.asarray(xys)[:,1], c=smb_lag_amax, cmap=lag_colors,
#                   vmin=lagnorm_min, vmax=lagnorm_max)
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
sc5 = ax5.scatter(np.asarray(xys)[runoff_significance,0], np.asarray(xys)[runoff_significance,1], 
                  c=np.asarray(runoff_lag_amax)[runoff_significance], cmap=lag_colors, marker=sig_markers[0],
                  vmin=lagnorm_min, vmax=lagnorm_max)
ax5.scatter(np.asarray(xys)[np.invert(runoff_significance),0], np.asarray(xys)[np.invert(runoff_significance),1], 
            c=np.asarray(runoff_lag_amax)[np.invert(runoff_significance)], cmap=lag_colors, marker=sig_markers[1],
                  vmin=lagnorm_min, vmax=lagnorm_max)
# sc5 = ax5.scatter(np.asarray(xys)[:,0], np.asarray(xys)[:,1], c=runoff_lag_amax, cmap=lag_colors,
#                   vmin=lagnorm_min, vmax=lagnorm_max)
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
sc6 = ax6.scatter(np.asarray(xys)[terminus_significance,0], np.asarray(xys)[terminus_significance,1], 
                  c=np.asarray(terminus_lag_amax)[terminus_significance], cmap=lag_colors, marker=sig_markers[0],
                  vmin=lagnorm_min, vmax=lagnorm_max)
ax6.scatter(np.asarray(xys)[np.invert(terminus_significance),0], np.asarray(xys)[np.invert(terminus_significance),1], 
            c=np.asarray(terminus_lag_amax)[np.invert(terminus_significance)], cmap=lag_colors, marker=sig_markers[1],
                  vmin=lagnorm_min, vmax=lagnorm_max)
# sc6 = ax6.scatter(np.asarray(xys)[:,0], np.asarray(xys)[:,1], c=terminus_lag_amax, cmap=lag_colors,
#                   vmin=lagnorm_min, vmax=lagnorm_max)
## set up correctly scaled colorbar
div6 = make_axes_locatable(ax6)
cax6 = div6.append_axes("right", size="5%", pad=0.1)
cb6 = fig.colorbar(sc6, cax=cax6)
cb6.ax.set_ylabel('Lag [d] at peak xcorr')
ax6.set(xlim=(278000, 320000), xticks=(280000, 300000, 320000), 
      ylim=(-2590000, -2550000), yticks=(-2590000, -2570000, -2550000), 
       xticklabels=('280', '300', '320'), yticklabels=('-2590', '-2570', '-2550'),
      xlabel='Easting [km]', aspect=1.)
# plt.tight_layout()
# plt.show()
# plt.savefig('/Users/lizz/Desktop/20210204-helheim-xcorr_lag_composite')
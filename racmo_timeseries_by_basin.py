#!/usr/bin/env python

import os, sys, argparse
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta
#from subprocess import call
import numpy as np
from netCDF4 import Dataset
#from pyproj import Proj, transform

import glob

#from scipy.interpolate import interp2d
from scipy import ndimage
from scipy.io import savemat, loadmat

import pickle

sys.path.append(os.environ['PYTHON_modules_DIR'])
import raster
import write_netcdf

import copy

import multiprocessing as mp
from progress.bar import Bar
import time

# ----------------
# Nested functions
# ----------------
def read_runoff(ncfilename,xVar,yVar,RACMOvar,runoffCalculation):
#{{{
   # Read file
   ncfile = Dataset(ncfilename, 'r')

   x = ncfile.variables[xVar][:]
   y = ncfile.variables[yVar][:]
   lat = ncfile.variables['LAT'][:,:]
   lon = ncfile.variables['LON'][:,:]
   runoff = ncfile.variables[RACMOvar]
   
   # Extract year
   year = int(os.path.basename(ncfilename).split('-')[-1].split('.')[0])
   #year = datetime.strptime(os.path.basename(ncfilename), netcdf_filename_datetime_template).year

   if runoffCalculation == 'monthly':
      # Monthly means
      runoff_mean = runoff[:,:,:]
      for month in range(1,13):
         dt = (datetime.strptime('{:02d}/1/{:4.0f}'.format(month,year), '%m/%d/%Y') - timeEpoch).days
         t0 = datetime.strptime('{:02d}/01/{:4.0f}'.format(month,year), '%m/%d/%Y')
         t1 = t0 + relativedelta(months=1)
         time_bounds.append( np.array( [[(t0 - timeEpoch).days, \
                                         (t1 - timeEpoch).days]] ).T)

   if runoffCalculation == 'yearly-sum':
      # Annual mean
      runoff_mean = np.expand_dims(np.sum(runoff[:,:,:], axis=0), axis=0) # runoff / year
      dt = (datetime.strptime('07/01/{:4.0f}'.format(year), '%m/%d/%Y') - timeEpoch).days
      time_bounds.append( np.array( [[(datetime.strptime('01/01/{:4.0f}'.format(year), '%m/%d/%Y') - timeEpoch).days, \
                                      (datetime.strptime('12/31/{:4.0f}'.format(year), '%m/%d/%Y') - timeEpoch).days]] ).T)

   if runoffCalculation == 'JJA-mean':
      # JJA mean
      runoff_mean = np.expand_dims(np.mean(runoff[5:8,:,:], axis=0), axis=0)
      dt = (datetime.strptime('07/01/{:4.0f}'.format(year), '%m/%d/%Y') - timeEpoch).days
      time_bounds.append( np.array( [[(datetime.strptime('06/01/{:4.0f}'.format(year), '%m/%d/%Y') - timeEpoch).days, \
                                      (datetime.strptime('08/31/{:4.0f}'.format(year), '%m/%d/%Y') - timeEpoch).days]] ).T)

   return dt, time_bounds, x, y, runoff_mean, year
#}}}

def calculate_bias(basinNum, basinArrayBaseline, basinArray, runoffBaseline, runoff):
#{{{
   #runoff_bl_basin = np.where(basinArrayBaseline==basinNum, np.flip(runoffBaseline, 1), 0.)
   #runoff_bl_basin_timeseries = np.sum(runoff_bl_basin, axis=(1,2))
   #runoff_pd_basin = np.where(basinArray==basinNum, np.flip(runoff, 1), 0.)
   #runoff_pd_basin_timeseries = np.sum(runoff_pd_basin, axis=(1,2))

   mask2d = basinArrayBaseline==basinNum
   mask3d = np.repeat(mask2d[np.newaxis,:,:], runoffBaseline.shape[0], axis=0)
   runoff_bl_basin = mask3d * np.flip(runoffBaseline, 1)
   mask2d = basinArray==basinNum
   mask3d = np.repeat(mask2d[np.newaxis,:,:], runoff.shape[0], axis=0)
   runoff_pd_basin = mask3d * np.flip(runoff, 1)
   runoff_bl_basin_timeseries = np.sum(runoff_bl_basin, axis=(1,2))
   runoff_pd_basin_timeseries = np.sum(runoff_pd_basin, axis=(1,2))

   # Yearly sums
   runoff_bl_basin_annual = 12. * moving_average.moving_average_downsample(range(0,len(runoff_bl_basin_timeseries)), runoff_bl_basin_timeseries, 12)
   runoff_pd_basin_annual = 12. * moving_average.moving_average_downsample(range(0,len(runoff_pd_basin_timeseries)), runoff_pd_basin_timeseries, 12)

   # Differences in annual means
   runoff_annual_diffs = runoff_pd_basin_annual - runoff_bl_basin_annual

   runoffBias = np.mean(runoff_annual_diffs)

   return runoffBias
#}}}

def runoff_unit_conversion(inputUnits, outputUnits):
#{{{
   # This assumes that the grid resolution is 1 km x 1 km
   if inputUnits == 'mmWE' and outputUnits == 'm3':
      unitConversion = 1000.
   elif inputUnits == 'mmWE yr-1' and outputUnits == 'kg s-1':
      unitConversion = 1000000./31557600.
   else:
      print('runoff unit conversion from ' + inputUnits + ' to ' + outputUnits + ' not supported!')
      return None
   return unitConversion
#}}}

def submerged_area(submergedAreaArray, basinArray, basinNum):
#{{{
   submergedAreaMasked = (basinArray==basinNum) * submergedAreaArray
   submergedArea = np.unique(submergedAreaMasked)
   submergedArea = submergedArea[submergedArea > 0.]
   if len(submergedArea) != 1:
      return None
   return submergedArea[0]
#}}}

# -------------
# --- Setup ---
# -------------
netcdf_dirs = list()

# Runoff
netcdf_dirs.append(os.environ['RACMO_downscaled_DIR']); netcdf_filename_template = 'runoff.????_???.BN_RACMO2.3p2_ERA5_3h_FGRN055.1km.DD.nc';
#netcdf_dirs.append(os.environ['RACMO_downscaled_DIR']); netcdf_filename_template = 'runoff_WJB_int.????.BN_1958_2013_1km.DD.nc';
tVar = 'time'
xVar = 'x'
yVar = 'y'
RACMOvar = 'runoffcorr'
projection = None # default is EPSG 3413, a.k.a. '+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs'
xOverride = None
yOverride = None
# Conversion from mm W.E. to cubic m (using RACMO grid resolution)
outputUnits = 'm3'
runoffUnitConversion = 1000.

# # SMB
# netcdf_dirs.append(os.environ['RACMO_downscaled_DIR']); netcdf_filename_template = 'smb_rec.????_???.BN_RACMO2.3p2_ERA5_3h_FGRN055.1km.DD.nc';
# #netcdf_dirs.append(os.environ['RACMO_downscaled_DIR']); netcdf_filename_template = 'SMB_rec_WJB_int.????.BN_1958_2013_1km.DD.nc';
# tVar = 'time'
# xVar = 'x'
# yVar = 'y'
# RACMOvar = 'smb_rec'
# #RACMOvar = 'SMB_rec'
# projection = None # default is EPSG 3413, a.k.a. '+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs'
# xOverride = None
# yOverride = None
# # Conversion from mm W.E. to cubic m (using RACMO grid resolution)
# outputUnits = 'mmWE' # 'kg s-1' | 'm3'
# runoffUnitConversion = 1

# Output file
outputDirectory = '.'

# Drainage basins
#{{{
clipfiledir = '~/ISMIP6/GrIS_Ocean_Forcing/Tidewater_Basins/tif'
#if basinsRaster == 'tidewaterbasins_rignotid':
#   clipfile = clipfiledir + '/tidewaterbasins_rignotid.mat_tidewaterbasins.tif'
#   submergedAreaRaster = None
#if basinsRaster == 'basins4highres':
clipfile = clipfiledir + '/basins4highres_xy.mat_basins_select.tif'
submergedAreaRaster = clipfiledir + '/basins4highres_xy.mat_submergedarea.tif'
outputUnits = outputUnits + ' m-2'

# HEL
basinfilter = 'basin=26';
basinValue = 26

# All basins
#clipfile = clipfiledir + '/tidewaterbasins_rignotid.mat_tidewaterbasins.tif'; basinfilter = 'basin=*';
basinfilter = 'basin=*'

# Parallel processing
nprocesses = 8
#}}}

# ------------------
# --- Processing ---
# ------------------
print("Setup:")
for netcdf_dir in netcdf_dirs:
      print(" input netcdf dir(s): " + netcdf_dir)
print(" basin file         : " + clipfile)
print(" output directory   : " + outputDirectory)
print(" basin filter       : " + basinfilter)
print(" ")

print("Processing:")
print(" reading netcdf(s)")
dts = np.array([])
time_bounds = list()
runoffs = np.array([])

x = None
y = None
lat = None
lon = None

# Obtain full list of NetCDF files
ncfilenames = []
for netcdf_dir in sorted(netcdf_dirs):
   ncfilenames.extend(sorted(glob.glob(netcdf_dir + '/' + netcdf_filename_template)))

# Read the first file to get spatial coordinates
ncfile = Dataset(ncfilenames[0], 'r')
if xOverride is None: x = ncfile.variables[xVar][:]
else:                 x = xOverride
if yOverride is None: y = ncfile.variables[yVar][:]
else:                 y = yOverride
ncfile.close()

# Geotransform
xStep = x[1]-x[0]
yStep = y[1]-y[0]
geoTransform = [x[0]-xStep/2, xStep, 0., y[-1]+yStep/2, 0., -yStep]

# Create basin mask(s) from clip file
# Prelim checks {{{
if clipfile.split('.')[-1] == 'shp' and not checkForField(clipfile, basinfilter.split('=')[0]):
   print('ERROR in ' + __file__ + ': Attribute "' + basinfilter.split('=')[0] + '" not found in ' + clipfile)
   sys.exit()
#}}}
if clipfile.split('.')[-1] == 'shp': #{{{
   print(" ")
   print("\033[93mWARNING! This script uses the exterior of each polygon! \033[0m") 
   print(" ")

   if '*' in basinfilter:
      basinValues = getUniqueFieldValues(clipfile, basinfilter.split('=')[0])
   else:
      basinValues = [basinfilter.split('=')[1]]

   # Setup mask array(s)
   if projection:
      print("\033[93mERROR! This script is not set up to handle clipfiles that are shp but have a different projection than climate data! Exiting. \033[0m") 
      sys.exit()
   maskArray = np.zeros( (len(basinValues), runoff_nrows, runoff_ncols) )

   # Create basin mask
   print(" creating basin masks")
   driver = ogr.GetDriverByName("ESRI Shapefile")
   dataSource = driver.Open(clipfile, 0)
   layer = dataSource.GetLayer()

   for iValue, basinValue in enumerate(basinValues):
      AF = basinfilter.split('=')[0] + '=' + str(basinValue)
      layer.SetAttributeFilter(AF)
      for feature in layer:
         geometry = feature.GetGeometryRef()
         ring = geometry.GetGeometryRef(0)
         xClip = list()
         yClip = list()
         for ixy, xy in enumerate(ring.GetPoints()):
            xClip.append(xy[0])
            yClip.append(xy[1])
            
         maskArray_feature = clipImage(np.ones( (runoff_nrows, runoff_ncols) ), xClip, yClip, geoTransform)
         maskArray_feature[np.isnan(maskArray_feature)] = 0.
         maskArray[iValue,:,:] = maskArray[iValue,:,:] + maskArray_feature
#}}}
elif clipfile.split('.')[-1] == 'tif': #{{{
   # Resample mask to reslution of climate data
   print(" ")
   print("\033[93mWARNING! This script resamples the basin mask to the resolution of the climate data using nearest-neighbor! \033[0m") 
   print(" ")

   # Setup mask array(s)
   files = list()
   files.append(clipfile)
   if submergedAreaRaster: files.append(submergedAreaRaster)
   for f in files:
      f_resampled = os.path.basename(f).replace('.tif','_resampled.tif')
      if os.path.isfile(f_resampled):
         print("\033[93mWARNING! Overwriting " + f_resampled + "\033[0m") 
         os.remove(f_resampled)
      if projection:
         print("\033[93mWARNING! Reprojecting to specified projection of climate data:\n\t" + projection + "\033[0m") 
         cmdStr = 'gdalwarp {:s} {:s} -r near -t_srs "{:s}" -ot UInt32 -of GTiff -te {:f} {:f} {:f} {:f} -tr {:f} {:f}'.format( \
               f, f_resampled, \
               projection, \
               np.min(x)-xStep/2, np.min(y)-yStep/2, np.max(x)+xStep/2, np.max(y)+yStep/2, \
               geoTransform[1], -geoTransform[5])
      else: 
         cmdStr = 'gdal_translate {:s} {:s} -r near -ot UInt32 -of GTiff -projwin {:f} {:f} {:f} {:f} -tr {:f} {:f}'.format( \
               f, f_resampled, \
               np.min(x)-xStep/2, np.max(y)+yStep/2, np.max(x)+xStep/2, np.min(y)-yStep/2, \
               geoTransform[1], -geoTransform[5])
      os.system(cmdStr)

   basinArray         = raster.readRasterBandAsArray(os.path.basename(clipfile).replace('.tif','_resampled.tif'), 1)
   if submergedAreaRaster: submergedAreaArray = raster.readRasterBandAsArray(os.path.basename(submergedAreaRaster).replace('.tif','_resampled.tif'), 1)
   else: submergedAreaArray = np.ones(basinArray.shape)
   
   #print(" ")
   #print("DEBUG")
   #basinArrayOrig = raster.readRasterBandAsArray(clipfile, 1)
   #submergedAreaArrayOrig = raster.readRasterBandAsArray(submergedAreaRaster, 1)
   #basinArrayOrig[basinArrayOrig==0.] = np.nan
   #basinValuesOrig = sorted(set(basinArrayOrig[~np.isnan(basinArrayOrig)]))
   #bar = Bar('Processing', max=len(basinValuesOrig))
   #for b in basinValuesOrig:
   #   check = np.unique(submergedAreaArrayOrig[basinArrayOrig==b])
   #   if len(check[~np.isnan(check)]) != 1:
   #      bar.finish()
   #      print('PROBLEM with basin {:.0f}'.format(b))
   #      import pdb; pdb.set_trace()
   #   bar.next()
   #bar.finish()
   #import pdb; pdb.set_trace()

   if '*' in basinfilter:
      basinValues = ['{:.0f}'.format(b) for b in sorted(set(basinArray[~np.isnan(basinArray)]))]
      basinValues.remove('0')
   else:
      basinValues = [basinfilter.split('=')[1]]

   ## DEBUG -- missing basins are the ones that are too small and discarded in the nearest-neighbor
   #origBasins = ['{:.0f}'.format(b) for b in sorted(set(basinArray[~np.isnan(basinArray)]))]
   #origBasins.remove('0')
   #origBasinArray = raster.readRasterBandAsArray(clipfile, 1)
   #for b in set(origBasins).symmetric_difference(set(basinValues)):
   #   print('{:s} {:f}'.format(b, 150.*150.*len(np.where(origBasinArray==float(b))[0])/1000./1000.))
   #import pdb; pdb.set_trace()

   ## Create basin mask
   #print(" creating basin masks")
   #maskArray = np.zeros( (len(basinValues), runoff_nrows, runoff_ncols) )
   #for iValue, basinValue in enumerate(basinValues):
   #   basinValue = float(basinValue)
   #   maskArray[iValue,:,:] = np.where(basinArray==basinValue, 1, 0)
#}}}
else: #{{{
   print(" clipfile file type not supported")
   sys.exit()
# }}}

# Read all netCDFs
dts = list()
runoffSums = list()
#{{{
for ncfilename in ncfilenames:
   print(ncfilename)
   # Read file
   ncfile = Dataset(ncfilename, 'r')
   time = ncfile.variables[tVar][:]
   runoff = ncfile.variables[RACMOvar][:,:,:]

   if hasattr(ncfile.variables[tVar], 'units'):
      time_epoch = datetime.strptime(ncfile.variables[tVar].units.lower().replace('days since ',''), '%Y-%m-%d %H:%M:%S')
   else:
      time_epoch = datetime.strptime('2015-01-01 00:00:00', '%Y-%m-%d %H:%M:%S')

   # Spatial processing
   basinNum = float(basinValue)
   runoffSum = np.empty(runoff.shape[0])
   for itime in range(runoff.shape[0]):
      runoffMasked = (basinArray==basinNum) * np.flipud(runoff[itime,:,:])

      # Sum runoff over the basin
      runoffSums.append(np.nansum(runoffMasked))
      dts.append(time_epoch + timedelta(days=int(time[itime])))

   ncfile.close()
#}}}

# Sort by time
dts = np.array(dts)
runoffSums = np.array(runoffSums)
sort_idx = np.argsort(dts)
dts_sorted = dts[sort_idx]
runoffSums_sorted = runoffSums[sort_idx]

# Write output text file
outputFile = outputDirectory + '/' + netcdf_filename_template.replace('?','').replace('.nc','.csv')
f = open(outputFile, 'w')
f.write('date, runoff ({:s})'.format(outputUnits))
for i in range(len(dts)):
   f.write('{:s}, {: 16.3f}\n'.format(dts_sorted[i].strftime('%Y-%m-%d %H:%M:%S'), runoffSums_sorted[i]))
f.close()


#!/Users/dfelikso/Research/Software/miniconda/miniconda3/envs/freshwater/bin/python

import sys
sys.path.append('/Users/dfelikso/Research/Data/GlacierDrainageBasins/Mankoff/freshwater')
from discharge import discharge

from shapely.geometry import mapping, Polygon
import fiona

from matplotlib import pyplot as plt

# Outlets / basins
df = discharge(base='/Users/dfelikso/Research/Data/GlacierDrainageBasins/Mankoff/basins1.0/freshwater', roi='292789,-2561824', quiet=True).outlets()

# Write shapefile
schema = {
    'geometry': 'Polygon',
    'properties': {'id': 'int'},
}
poly = df.geometry[1]
with fiona.open('data/helheim_ice_catchment_mankoff.shp', 'w', 'ESRI Shapefile', schema) as c:
   ## If there are multiple geometries, put the "for" loop here
   c.write({
       'geometry': mapping(poly),
       'properties': {'id': -9999},
   })

#x, y = df.geometry[0].exterior.xy
#plt.plot(x, y)
#plt.show()

# Discharge
#ds = discharge(base='/Users/dfelikso/Research/Data/GlacierDrainageBasins/Mankoff/basins1.0/freshwater', roi='292789,-2561824', quiet=True).discharge()


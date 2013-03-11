"""
Access data from the NECOFS (New England Coastal Ocean Forecast System) via OPeNDAP

Demonstration using the NetCDF4-Python library to access velocity data
from a triangular grid ocean model (FVCOM) via OPeNDAP, specifying the
desired URL, time, layer and lat/lon region of interest.  The
resulting plot of forecast velocity vectors over color-shaded
bathymetry is useful for a variety of recreational and scientific
purposes.

NECOFS (Northeastern Coastal Ocean Forecast System) is run by groups
at the University of Massachusetts Dartmouth and the Woods Hole
Oceanographic Institution, led by Drs. C. Chen, R. C. Beardsley,
G. Cowles and B. Rothschild. Funding is provided to run the model by
the NOAA-led Integrated Ocean Observing System and the State of
Massachusetts.

NECOFS is a coupled numerical model that uses nested weather models, a
coastal ocean circulation model, and a wave model. The ocean model is
a volume-mesh model with horizontal resolution that is finer in
complicated regions. It is layered (not depth-averaged) and includes
the effects of tides, winds, and varying water densities caused by
temperature and salinity changes.

* Model description: http://fvcom.smast.umassd.edu/research_projects/NECOFS/model_system.html
* THREDDS server with other forecast and archive products: http://www.smast.umassd.edu:8080/thredds/catalog.html
"""

from pylab import *
import matplotlib.tri as Tri
import netCDF4
import datetime as dt

url = 'http://www.smast.umassd.edu:8080/thredds/dodsC/FVCOM/NECOFS/Forecasts/NECOFS_FVCOM_OCEAN_MASSBAY_FORECAST.nc'
# Open DAP
nc = netCDF4.Dataset(url)
nc.variables.keys()

# take a look at the "metadata" for the variable "u"
print nc.variables['u']

# Desired time for snapshot
# ....right now (or some number of hours from now) ...
start = dt.datetime.utcnow() + dt.timedelta(hours=0)
# ... or specific time (UTC)
start = dt.datetime(2013,3,2,15,0,0)

# Get desired time step  
time_var = nc.variables['time']
itime = netCDF4.date2index(start,time_var,select='nearest')

# Get lon,lat coordinates for nodes (depth)
lat = nc.variables['lat'][:]
lon = nc.variables['lon'][:]
# Get lon,lat coordinates for cell centers (depth)
latc = nc.variables['latc'][:]
lonc = nc.variables['lonc'][:]
# Get Connectivity array
nv = nc.variables['nv'][:].T - 1 
# Get depth
h = nc.variables['h'][:]  # depth 

dtime = netCDF4.num2date(time_var[itime],time_var.units)
daystr = dtime.strftime('%Y-%b-%d %H:%M')
print daystr

tri = Tri.Triangulation(lon,lat, triangles=nv)

# get current at layer [0 = surface, -1 = bottom]
ilayer = 0
u = nc.variables['u'][itime, ilayer, :]
v = nc.variables['v'][itime, ilayer, :]

levels=arange(-32,2,1)   # depth contours to plot
ax= [-70.97, -70.82, 42.25, 42.35] # region to plot

# find velocity points in bounding box
ind = argwhere((lonc >= ax[0]) & (lonc <= ax[1]) & (latc >= ax[2]) & (latc <= ax[3]))

subsample=3
np.random.shuffle(ind)
Nvec = int(len(ind) / subsample)
idv = ind[:Nvec]

# tricontourf plot of water depth with vectors on top
fig1 = figure(figsize=(18,10))
ax1 = fig1.add_subplot(111,aspect=(1.0/cos(mean(lat)*pi/180.0)))
plt.tricontourf(tri, -h,levels=levels,shading='faceted',cmap=plt.cm.gist_earth)
plt.axis(ax)
ax1.patch.set_facecolor('0.5')
cbar=colorbar()
cbar.set_label('Water Depth (m)', rotation=-90)
Q = ax1.quiver(lonc[idv],latc[idv],u[idv],v[idv],scale=20)
qk = quiverkey(Q,0.92,0.08,0.50,'0.5 m/s',labelpos='W')
title('NECOFS Velocity, Layer %d, %s' % (ilayer, daystr))

plt.show()

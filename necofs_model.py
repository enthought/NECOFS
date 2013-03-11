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

# standard library imports
import datetime as dt

# Major library imports
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as Tri
import netCDF4

# Enthought library imports
from traits.api import HasTraits, List, Str, Int, Instance


class OceanModel(HasTraits):

    url = Str('http://www.smast.umassd.edu:8080/thredds/dodsC/FVCOM/NECOFS/Forecasts/NECOFS_FVCOM_OCEAN_MASSBAY_FORECAST.nc')
    keys = List
    ilayer = Int(0)
    start = Instance(dt.datetime)
    daystr = Str

    def open(self):
        # Open DAP
        nc = netCDF4.Dataset(self.url)
        self.keys = nc.variables.keys()
        start = dt.datetime.utcnow()

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
        self.daystr = dtime.strftime('%Y-%b-%d %H:%M')

        tri = Tri.Triangulation(lon,lat, triangles=nv)

        # get current at layer [0 = surface, -1 = bottom]
        u = nc.variables['u'][itime, self.ilayer, :]
        v = nc.variables['v'][itime, self.ilayer, :]

        levels=np.arange(-32,2,1)   # depth contours to plot
        ax= [-70.97, -70.82, 42.25, 42.35] # region to plot

        # find velocity points in bounding box
        ind = np.argwhere((lonc >= ax[0]) & (lonc <= ax[1]) & (latc >= ax[2]) & (latc <= ax[3]))

        subsample=3
        np.random.shuffle(ind)
        Nvec = int(len(ind) / subsample)
        idv = ind[:Nvec]
        self.lat, self.lon = lat, lon
        self.latc, self.lonc = latc, lonc
        self.tri = tri
        self.h = h
        self.idv = idv
        self.u, self.v = u, v
        self.levels = levels
        self.ax = ax

    def plot(self):
        # tricontourf plot of water depth with vectors on top
        fig1 = plt.figure(figsize=(18,10))
        ax1 = fig1.add_subplot(111,aspect=(1.0/np.cos(np.mean(self.lat)*np.pi/180.0)))
        plt.tricontourf(self.tri, -self.h,
                        levels=self.levels,
                        shading='faceted',
                        cmap=plt.cm.gist_earth)
        plt.axis(self.ax)
        ax1.patch.set_facecolor('0.5')
        cbar=plt.colorbar()
        cbar.set_label('Water Depth (m)', rotation=-90)
        Q = ax1.quiver(self.lonc[self.idv],
                       self.latc[self.idv],
                       self.u[self.idv],
                       self.v[self.idv],
                       scale=20)
        qk = plt.quiverkey(Q,0.92,0.08,0.50,'0.5 m/s',labelpos='W')
        plt.title('NECOFS Velocity, Layer %d, %s' % (self.ilayer, self.daystr))

        plt.show()

if __name__ == '__main__':
    model = OceanModel()
    model.open()
    model.plot()

    

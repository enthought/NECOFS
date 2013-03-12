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
import matplotlib
# We want matplotlib to use a QT backend
matplotlib.use('Qt4Agg')
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import matplotlib.tri as Tri
import netCDF4

# Enthought library imports
from traits.api import (HasTraits, List, Str, Int, Instance, Array, Property,
                        Range, Any, Float, cached_property, on_trait_change)
from enthought.traits.ui.api import View, Item
import enaml
from enaml.qt.qt_application import QtApplication
from enaml.stdlib.sessions import SimpleSession
from PySide import QtGui, QtCore
#from pyface.qt import QtGui, QtCore

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from enthought.traits.ui.qt4.editor import Editor
from enthought.traits.ui.qt4.basic_editor_factory import BasicEditorFactory


class _MPLFigureEditor(Editor):

    scrollable  = True

    def init(self, parent):
        self.control = self._create_canvas(parent)
        self.set_tooltip()

    def update_editor(self):
        pass

    def _create_canvas(self, parent):
        """ Create the MPL canvas. """
        # matplotlib commands to create a canvas
        mpl_canvas = FigureCanvas(self.value)
        return mpl_canvas


class MPLFigureEditor(BasicEditorFactory):

    klass = _MPLFigureEditor


class OceanModel(HasTraits):

    #url = Str('http://www.smast.umassd.edu:8080/thredds/dodsC/FVCOM/NECOFS/Forecasts/NECOFS_FVCOM_OCEAN_MASSBAY_FORECAST.nc')
    url = Str('10step_nc4.nc')
    nc = Any()
    keys = List
    ilayer = Int(0)
    start = Instance(dt.datetime)
    daystr = Property(Str, depends_on='itime')
    itime = Range(0,10)
    time_var = Any()
    lat = Array
    lon = Array
    latc = Array
    lonc = Array
    tri = Any()
    nv = Array
    h = Array
    u = Property(Array)
    v = Property(Array)
    levels = Array
    west = Float(-70.97)
    east = Float(-70.82)
    south = Float(42.25)
    north = Float(42.35)
    ax = Property(List, depends_on='west, east, south, north')
    ind = Property(Array, depends_on='west, east, south, north')
    idv = Property(Array, depends_on='ind')
    figure = Instance(Figure)
    #figure = Property(Instance(Figure), depends_on='ax, idv')
    quiver = Any()
    # TraitsUI view
    view = View(Item('figure', editor=MPLFigureEditor(),
                     show_label=False),
                Item('itime'),
                Item('north'),
                Item('south'),
                Item('west'),
                Item('east'),
                width=800,
                height=800,
                resizable=True)

    def _nc_default(self):
        """ Open DAP
        """
        return netCDF4.Dataset(self.url)

    def _keys_default(self):
        self.keys = self.nc.variables.keys()

    def _set_start(self, start):
        """ Desired time slice
        """
        self.itime = netCDF4.date2index(self.start, self.time_var, select='nearest')
        self._start = start

    def _get_start(self):
        try:
            return self._start
        except AttributeError:
            return dt.datetime.utcnow()

    def _time_var_default(self):
        return self.nc.variables['time']

    def _lat_default(self):
        """ Latitude of grid nodes
        """
        return self.nc.variables['lat'][:]

    def _lon_default(self):
        """ Longitude of grid nodes
        """
        return self.nc.variables['lon'][:]

    def _latc_default(self):
        """ Latitude of cell centers (depth)
        """
        return self.nc.variables['latc'][:]

    def _lonc_default(self):
        """ Longitude of cell centers (depth)
        """
        return self.nc.variables['lonc'][:]

    def _nv_default(self):
        """ Get Connectivity array
        """
        return self.nc.variables['nv'][:].T - 1

    def _h_default(self):
        """ Get depth
        """
        return self.nc.variables['h'][:]

    def _get_daystr(self):
        """ String representation of current time stamp
        """
        dtime = netCDF4.num2date(self.time_var[self.itime],
                                 self.time_var.units)
        return dtime.strftime('%Y-%b-%d %H:%M')

    def _tri_default(self):
        return Tri.Triangulation(self.lon, self.lat, triangles=self.nv)

    def _get_u(self):
        # get current at layer [0 = surface, -1 = bottom]
        return self.nc.variables['u'][self.itime, self.ilayer, :]

    def _get_v(self):
        return self.nc.variables['v'][self.itime, self.ilayer, :]

    @on_trait_change('itime')
    def update_quiver(self):
        """ set the quiver plot data without redrawing
        """
        self.quiver.set_UVC(self.u[self.idv], self.v[self.idv])
        self.figure.canvas.draw()

    def _levels_default(self):
        """ depth contours to plot
        """
        return np.arange(-32, 2, 1)

    def _get_ax(self):
        """ region to plot
        """
        return [self.west, self.east, self.south, self.north]

    @cached_property
    def _get_ind(self):
        """ find velocity points in bounding box
        """
        print 'get ind'
        return np.argwhere((self.lonc >= self.west) &
                           (self.lonc <= self.east) &
                           (self.latc >= self.south) &
                           (self.latc <= self.north))

    @cached_property
    def _get_idv(self):
        print 'get idv'
        subsample = 3
        np.random.shuffle(self.ind)
        Nvec = int(len(self.ind) / subsample)
        return self.ind[:Nvec]

    @on_trait_change('west, east, south, north')
    def update_plot(self, axis=None):
        """ 
        """
        if axis is None and self.figure is not None:
            axis = self.figure.axes[0]
        #print self.figure
        axis.axis(self.ax)
        axis.patch.set_facecolor('0.5')
        #cbar = plt.colorbar()
        #cbar.set_label('Water Depth (m)', rotation=-90)
        self.quiver = axis.quiver(self.lonc[self.idv],
                       self.latc[self.idv],
                       self.u[self.idv],
                       self.v[self.idv],
                       scale=20)
        axis.quiverkey(self.quiver, 0.92, 0.08, 0.50, '0.5 m/s', labelpos='W')
        #plt.title('NECOFS Velocity, Layer %d, %s' % (self.ilayer, self.daystr))

    def _figure_default(self):
        # tricontourf plot of water depth with vectors on top
        print 'get figure'
        fig1 = Figure(figsize=(600,600), dpi=72, facecolor=(1,1,1), edgecolor=(0,0,0))
        ax1 = fig1.add_subplot(111, aspect=(1.0/np.cos(np.mean(self.lat)*np.pi/180.0)))
        ax1.tricontourf(self.tri, -self.h,
                        levels=self.levels,
                        shading='faceted',
                        cmap=plt.cm.gist_earth)
        self.update_plot(ax1)
        return fig1

if __name__ == '__main__':

    model = OceanModel()
    model.configure_traits()

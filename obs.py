#!/usr/bin/env python
#encoding:utf-8

import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
try:
    from shapely import geometry
    has_shapely = True
except:
    has_shapely = False

class Obs(object):

    def __init__(self):
        self.data = dict()

    def d(self, name):
        return self.data[name]

    def set_data(self, name, data):
        self.data[name] = data

    def set_invalid(self, invalid):
        for name in self.data:
            data = self.data[name]
            self.data[name] = ma.masked_where(data==invalid, data)

class SatelliteObs(Obs):

    def set_coords(self, lon, lat):
        self.lon, self.lat = lon, lat

    def coords(self):
        return self.lon, self.lat

    def set_invalid(self, invalid):
        self.invalid = invalid
        Obs.set_invalid(self, invalid)
        # self.lon = ma.masked_where(self.lon==invalid, self.lon)
        # self.lat = ma.masked_where(self.lat==invalid, self.lat)

    def points_in_domain(self, domain):
        '''
        Returns a boolean field that identifies the points inside the given domain.
        Domain: a polygon from shapely.geometry.Polygon
        '''

        idx = (self.lon > -9000) & (self.lat > -9000)
        traj = geometry.asMultiPoint(list(zip(self.lon[idx], self.lat[idx])))
        hits = domain.intersection(traj)
        if hits.is_empty:
            return None

        hits_array = np.array(hits)
        if np.shape(hits_array)==() or hits_array.ndim==1:
            return None

        idx = np.in1d(self.lat, hits_array[:,1])
        return idx


    def plot_lonlat(self, proj='cyl', boundinglat=None, lon_0=None, delay_show=False, idx=None, m=None, color='b'):
        if m is None:
            m = Basemap(projection=proj, boundinglat=boundinglat, lon_0=lon_0)
        if idx is not None:
            x, y = m(self.lon[idx], self.lat[idx])
        else:
            x, y = m(self.lon, self.lat)
        m.plot(x, y, '.', color=color)
        m.drawcoastlines()
        m.drawmapboundary()
        if not delay_show:
            plt.show()
        return m


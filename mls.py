#!/usr/bin/env python
# encoding: utf-8

'''

MLS functions
for MLS H2O files

V. Noel - Thu Jan 27 17:21:26 CET 2011
'''

import tables

class MLS(object):
    
    def __init__(self, filename, type):
        '''
        filename : name of the MLS file
        type : for now, 'HNO3' or 'H2O'
        '''
        self.h5 = tables.openFile(filename, mode='r')
        self.type = type
        self.swathnode = '/HDFEOS/SWATHS/'+self.type
        
    def close(self):
        self.h5.close()
        
    def coords(self):
        geo = self.h5.getNode(self.swathnode, 'Geolocation Fields')
        lat = geo.Latitude.read()
        lon = geo.Longitude.read()
        return lon, lat
        
    def levels(self):
        '''
        returns pressure levels
        '''
        geo = self.h5.getNode(self.swathnode, 'Geolocation Fields')
        levels = geo.Pressure.read()
        return levels
        
    def _L2gpValue(self):
        '''
        returns data array in MLS level 2 file
        shape (time, levels)
        '''
        data = self.h5.getNode(self.swathnode, 'Data Fields')
        h2o = data.L2gpValue.read()
        return h2o
        
    def precision(self):
        '''
        returns h2o precision
        shape h2o_precision(time, levels)
        '''
        data = self.h5.getNode(self.swathnode, 'Data Fields')
        precision = data.L2gpPrecision.read()
        return precision
        
    def quality(self):
        '''
        returns h2o quality
        shape quality(time)
        '''
        data = self.h5.getNode(self.swathnode, 'Data Fields')
        precision = data.Quality.read()
        return precision
        
    def time(self):
        geo = self.h5.getNode(self.swathnode, 'Geolocation Fields')
        time = geo.Time.read()
        return time
        
        
class MLSHNO3(MLS):

    def __init__(self, filename):
        MLS.__init__(self, filename, 'HNO3')

    def HNO3(self):
        data =self._L2gpValue()
        return data

class MLSH2O(MLS):

    def __init__(self, filename):
        MLS.__init__(self, filename, 'H2O')

    def H2O(self):
        data = self._L2gpValue()
        return data


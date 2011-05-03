#!/usr/bin/env python
# encoding: utf-8

'''

Module for reading MLS Level 2 file
Supported products : H2O, HNO3

V. Noel - Created on Thu Jan 27 17:21:26 CET 2011

'''

import tables

class MLS(object):
    
    '''
    Generic parent class to read MLS objects.
    Use product-specific classes instead.
    Example usage:
        m = MLS('MLS-Aura_L2GP-H2O_v02-23-c01_2010d238.he5', type='H2O')
        lon, lat = m.coords()
        levels = m.levels()
        m.close()
    '''
    
    def __init__(self, filename, type):
        '''
        filename : name of the MLS file
        type : product to read - supported products: 'HNO3', 'H2O'
        '''
        
        self.h5 = tables.openFile(filename, mode='r')
        self.type = type
        self.swathnode = '/HDFEOS/SWATHS/'+self.type
        
    def close(self):
        '''
        close an opened MLS object.
        '''
        self.h5.close()
        
    def coords(self):
        '''
        returns the coordinates of MLS profiles.
        shape: ntime
        '''
        geo = self.h5.getNode(self.swathnode, 'Geolocation Fields')
        lat = geo.Latitude.read()
        lon = geo.Longitude.read()
        return lon, lat
        
    def levels(self):
        '''
        returns pressure levels of MLS profiles
        shape: nlevels
        '''
        geo = self.h5.getNode(self.swathnode, 'Geolocation Fields')
        levels = geo.Pressure.read()
        return levels
        
    def _L2gpValue(self):
        '''
        returns data product containing MLS profiles
        shape: (ntime, nlevels)
        '''
        data = self.h5.getNode(self.swathnode, 'Data Fields')
        h2o = data.L2gpValue.read()
        return h2o
        
    def precision(self):
        '''
        returns MLS precision field
        shape: (ntime, nlevels)
        '''
        data = self.h5.getNode(self.swathnode, 'Data Fields')
        precision = data.L2gpPrecision.read()
        return precision
        
    def quality(self):
        '''
        returns MLS quality field
        shape: time
        '''
        data = self.h5.getNode(self.swathnode, 'Data Fields')
        precision = data.Quality.read()
        return precision
        
    def time(self):
        '''
        returns time for MLS profiles.
        shape: ntime
        '''
        geo = self.h5.getNode(self.swathnode, 'Geolocation Fields')
        time = geo.Time.read()
        return time
        
class MLSH2O(MLS):
    '''
    Class to read MLS H2O products.
    Example usage:
        m = MLSH2O('MLS-Aura_L2GP-H2O_v02-23-c01_2010d238.he5')
        lon, lat = m.coords()
        levels = m.levels()
        h2o = m.H2O()
        m.close()
    '''

    def __init__(self, filename):
        MLS.__init__(self, filename, 'H2O')

    def H2O(self):
        '''
        return numpy array containing the H2O product from the MLS file.
        shape: (ntime, nlevels)
        '''
        data = self._L2gpValue()
        return data

        
class MLSHNO3(MLS):
    '''
    Class to read MLS HNO3 products.
    '''

    def __init__(self, filename):
        MLS.__init__(self, filename, 'HNO3')

    def HNO3(self):
        data =self._L2gpValue()
        return data


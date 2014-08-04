#!/usr/bin/env python
#encoding:utf-8

# Created by VNoel Fri May 09

from pyhdf.SD import SD, SDC

class GeoProf(object):
    
    def __init__(self, filename):
        
        self.hdf = SD(filename)
        self.filename = filename
        self.z = filename[-2:]
        self.orbit = filename[-15:-4]
        self.id = filename[-25:-4]
        
    def close(self):
        self.hdf.end()
        self.hdf = None
        
    def _read_var(self, varname):
        hdfvar = self.hdf.select(varname)
        data = hdfvar[:]
        hdfvar.endaccess()
        return data
    
    def coords(self):
        lat = self._read_var('Latitude')
        lon = self._read_var('Longitude')
        return lon, lat
    
    def time(self):
        time = self._read_var('Time')
        return time
    
    
    def altitude(self):
        ''' altitude in kilometers '''
        
        alt = self._read_var('Height') / 1e3
        return alt
    
    def cloudmask(self):
        '''
        "0 = No cloud detected\n",
        "1 = likely bad data\n",
        "5 = likely ground clutter\n",
        "5-10 = week detection found using along track integration\n",
        "20 to 40 = Cloud detected .. increasing values represents clouds with lower chance of a being a false detection" ;
        _FillValue = '\200' ;
        
        shape [nprof, nalt]
        '''
        
        cm = self._read_var('CPR_Cloud_mask')
        
        return cm
    

def test_geoprof():
    
    path = '/bdd/CFMIP/OBS_LOCAL/ATRAIN_COLOC/CLOUDSAT_COLOC/CALTRACK-GEOPROF/2007/2007_01_01/'
    geofile = 'CALTRACK-5km_CS-2B-GEOPROF_V1-00_2007-01-01T02-01-45ZN.hdf'
    geo = GeoProf(path + geofile)

    assert geo.id == '2007-01-01T02-01-45ZN'
    assert geo.orbit == 'T02-01-45ZN'

    lon, lat = geo.coords()
    assert lon.shape == (3728,)
    assert lat.shape == (3728,)
    
    alt = geo.altitude()
    assert alt.max() < 40.
    assert alt.min() > -10.
    
    cm = geo.cloudmask()
    assert alt.shape == (3728, 125)
    assert cm.shape == (3728, 125)
        
    geo.close()
        
    
def show_example():

    import matplotlib.pyplot as plt
    
    path = '/bdd/CFMIP/OBS_LOCAL/ATRAIN_COLOC/CLOUDSAT_COLOC/CALTRACK-GEOPROF/2008/2008_01_01/'
    geofile = 'CALTRACK-5km_CS-2B-GEOPROF_V1-00_2008-01-01T01-30-23ZN.hdf'
    geo = GeoProf(path + geofile)

    plt.pcolormesh(geo.cloudmask().T[:,::-1])
    plt.ylim(125,0)
    plt.colorbar()
    plt.show()
    
if __name__ == '__main__':
    show_example()

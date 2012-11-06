#!/usr/bin/env python
#encoding:utf-8

import netCDF4
import numpy as np

# find matrix indexes that match coordinates the closest in a WRF domain

def main(wrffile='/mnt/cfmipfs/lov/homedata/noel/Projects/RedPanda/run.20060627.79571/WPS/geo_em.d01.nc', xlon=-60, xlat=-71):

    print 'Looking in ', wrffile
    print 'Requested coordinates : '
    print '  lon = ', xlon
    print '  lat = ', xlat

    nc = netCDF4.Dataset(wrffile)
    lon = nc.variables['XLONG_M'][:]
    lon = lon.squeeze()
    lat = nc.variables['XLAT_M'][:]
    lat = lat.squeeze()
    
    dist = np.power(lon-xlon, 2) + np.power(lat-xlat, 2)
    mindist = 1e10
    
    ni, nj = lon.shape
    imin, jmin = -1, -1
    
    for i in np.r_[0:ni]:
        for j in np.r_[0:nj]:
            if dist[i,j] < mindist:
                imin = i
                jmin = j
                mindist = dist[i,j]
                
    print 'Found'
    print '  lon = ', lon[imin,jmin]
    print '  lat = ', lat[imin,jmin]
    print 'at indexes ', imin, jmin

if __name__ == '__main__':
    import plac
    plac.call(main)
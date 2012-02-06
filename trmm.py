#!/usr/bin/env python
# encoding: utf-8
"""
trmm_ir

Created by Vincent Noel - LMD/CNRS on 2012-02-01.
"""

import numpy as np


def remap_lon(lon):
    imid = 9896/2
    lon2 = np.zeros_like(lon)
    lon2[imid:] = lon[:imid]
    lon2[:imid] = lon[imid:] - 360.
    return lon2
    

def remap_data(data):
    imid = 9896/2
    data2 = np.zeros_like(data)
    data2[:,imid:] = data[:,:imid]
    data2[:,:imid] = data[:,imid:]
    return data2


def read_trmm_ir(f):
    
    '''
    Read TRMM's ancillary IRBT products
    '''
    
    nlon = 9896
    nlat = 3298
    lonstep = 0.036378335
    latstep = 0.036383683
    lon0 = 0.0182
    lat0 = 59.982

    s = nlon * nlat
    
    if f.endswith('.Z'):

        # .Z files are compressed using the Unix compress command and are a PITA
        # to uncompress through 3rd party libraries
        print 'Opening compressed TRMM file through gunzip'

        import subprocess

        commandline = 'gunzip --stdout %s' % f

        p = subprocess.Popen(commandline, shell=True, stdout=subprocess.PIPE)
        uncompressed_data, stderr = p.communicate()

        data = np.fromstring(uncompressed_data, dtype=np.uint8, count=s)
    else:
        fd = open(f, 'rb')
        data = np.fromfile(file=fd, dtype=np.uint8, count=s)
        fd.close()
    
    data = data.reshape((nlat, nlon))
    data = remap_data(data) + 75.
        
    lon = np.linspace(lon0, lon0 + nlon * lonstep, num=nlon)
    lon = remap_lon(lon)
    
    lat = np.linspace(lat0, lat0 - nlat * latstep, num=nlat)
    
    return lon, lat, data
    

def main():
    pass


if __name__ == '__main__':
    main()


#!/usr/bin/env python
# encoding: utf-8

'''
Map utility functions
VNoel 2012
'''

from mpl_toolkits.basemap import Basemap

def northpole(boundinglat=60):

    bm = Basemap(projection='npstere', boundinglat=boundinglat, lon_0=0, resolution='i')
    return bm
    
def peninsula():
    
    bm = Basemap(projection='lcc', lat_0=-70, lon_0=-60, width=2e6, height=2e6, resolution='i')
    return bm
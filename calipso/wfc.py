#!/usr/bin/env python
#encoding:utf-8

# Created by V. Noel [LMD/CNRS] on 2016-06-09

import numpy as np
from calipso_hdf import _Cal
import seaborn

class WFC(_Cal):
    
    def __init__(self, filename):
        _Cal.__init__(self, filename)
        
    def time(self):
        time = self._read_var('Scan_UTC_Time')
        return time
        
    def coords(self): # 120560 x 40
        lon, lat = self._read_var('Longitude'), self._read_var('Latitude')
        return lon, lat
        
    def radiance(self): # 120560 x 40  watts/m^2/micrometer/steradian
        rad = self._read_var('Radiance')
        return rad
        
    def reflectance(self): #  120560 x 40
        ref = self._read_var('Reflectance')
        return ref
        
    def qc(self):
        qcx = self._read_var('Pixel_QC_Flag')
        return qcx


def main(filename='/bdd/CFMIP/OBS_LOCAL/WFC/2007_03_01/CAL_WFC_L1_125m-ValStage1-V3-01.2007-03-01T00-55-42ZD.hdf'):
    
    c = WFC(filename)
    time = c.time()
    lon, lat = c.coords()
    radiance = c.radiance()
    reflectance = c.reflectance()
    qc = c.qc()
    c.close()
    
    import matplotlib.pyplot as plt
    
    plt.figure(figsize=[25,3])
    plt.pcolormesh(radiance.T, cmap='viridis')
    cb = plt.colorbar()
    cb.set_label('radiance')

    plt.figure(figsize=[25,3])
    plt.pcolormesh(radiance.T, cmap='viridis')
    cb = plt.colorbar()
    cb.set_label('radiance')
    plt.xlim(35e3, 36e3)
    
    plt.figure()
    bins = np.r_[0:500:20]
    plt.hist(radiance.flatten(), bins=bins)
    plt.xlim(0, 500)
    plt.xlabel('Radiance')
    idx = (radiance.flatten() < 40)
    print '%f percents of points have radiance < 40' % (100. * np.sum(idx) / radiance.size)
    

    plt.figure(figsize=[25,3])
    plt.pcolormesh(reflectance.T, cmap='viridis')
    plt.clim(0,1)
    cb = plt.colorbar()
    cb.set_label('reflectance')

    plt.figure(figsize=[25,3])
    plt.pcolormesh(reflectance.T, cmap='viridis')
    plt.clim(0,1)
    cb = plt.colorbar()
    cb.set_label('reflectance')
    plt.xlim(35e3, 36e3)
    
    plt.figure()
    bins = np.r_[0:1.0:0.04]
    plt.hist(reflectance.flatten(), bins=bins)
    plt.xlim(0, 1.0)
    plt.xlabel('Reflectance')

    idx = (reflectance.flatten() < 0.2)
    print '%f percents of points have reflectance < 0.1' % (100. * np.sum(idx) / radiance.size)
    
    plt.show()

if __name__ == '__main__':
    import plac
    plac.call(main)
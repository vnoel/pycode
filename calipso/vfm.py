#!/usr/bin/env python
#encoding:utf-8

# Created by VNoel 2014-05-16 12:40

from calipso_hdf import _Cal

class VFM(_Cal):
    """
    Class to process CALIOP Level 2 Vertical Feature Masks.
    example use:
        
    >>> from calipso import VFM
        
    >>> vfm = VFM('CAL_LID_L2_VFM-ValStage1-V3-01.2010-03-04T01-32-02ZN.hdf')
    >>> time = vfm.time()
    >>> lon, lat = vfm.coords()
    >>> flags = vfm.flags()
        ...
    >>> vfm.close()
    """

    def __init__(self, filename):
        
        _Cal.__init__(self, filename)
        
    def time(self):
        
        time = self._read_var('Profile_Time')[:,0]
        return time
    
    def utc_time(self):
        
        time = self._read_var('Profile_UTC_Time')[:,0]
        return time
    
    def coords(self):
        
        lat = self._read_var('Latitude')[:,0]
        lon = self._read_var('Longitude')[:,0]
        return lon, lat
    
    def day_night_flags(self):
        
        flag = self._read_var('Day_Night_Flag')[:,0]
        return flag
    
    def land_water_mask(self):
        
        mask = self._read_var('Land_Water_Mask')[:,0]
        return mask
        
    def flags(self):
        '''
        reads the Feature_Classification_Flags array.
        This array needs to be decomposed to a 2D grid to be usable
        '''
        
        flags = self._read_var('Feature_Classification_Flags')[:,:]
        return flags
    
    
def test_read():
    
    testfile = '/users/noel/data/Data/VFM/CAL_LID_L2_VFM-ValStage1-V3-01.2010-03-04T01-32-02ZN.hdf'
    
    vfm = VFM(testfile)
    time = vfm.time()
    lon, lat = vfm.coords()
    flags = vfm.flags()
    vfm.close()
        
    assert time.shape[0] == 3744
    assert lon.shape[0] == 3744
    assert flags.shape[0] == 3744
    assert flags.shape[1] == 5515

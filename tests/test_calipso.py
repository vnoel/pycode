#!/usr/bin/env python
#encoding:utf-8

import pytest
import numpy as np
import calipso

test_l1_file='/bdd/CALIPSO/Lidar_L1/CAL_LID_L1.v3.30/2013/2013_04_07/CAL_LID_L1-ValStage1-V3-30.2013-04-07T10-26-26ZN.hdf'

nz = 583
nz_met = 33
nprof = 56190

@pytest.fixture(scope="module")
def cal1file():
    c = calipso.Cal1(test_l1_file)
    return c


def test_date(cal1file):
    
    assert cal1file.year==2013
    assert cal1file.month == 4
    assert cal1file.day == 7
    

def load_vectors(cal1file, navg=1):
    
    lon, lat = cal1file.coords(navg=navg)
    time = cal1file.time(navg=navg)
    utc  = cal1file.utc_time(navg=navg)
    elev = cal1file.surface_elevation(navg=navg)
    return lon, lat, time, utc, elev

    
def test_nprof(cal1file):
    
    vectors = load_vectors(cal1file,navg=1)
    for vector in vectors:
        assert len(vector)==nprof
    vectors = load_vectors(cal1file, navg=30)
    for vector in vectors:
        assert len(vector)==np.floor(1.*nprof/30)
        
    
def test_nz(cal1file):
    
    alt = cal1file.lidar_alt
    assert len(cal1file.lidar_alt)==nz
    assert len(cal1file.met_alt)==nz_met
    
    
def test_arrays(cal1file, navg=30):
    
    atb = cal1file.atb(navg=navg)
    assert atb.shape[0]==np.floor(1.*nprof/navg)
    assert atb.shape[1]==nz
    
    mol = cal1file.mol(navg=navg)
    temp = cal1file.temperature(navg=navg)
    for var in mol,temp:
    	assert var.shape[0]==np.floor(1.*nprof/navg)
    	assert var.shape[1]==nz_met
    

if __name__=='__main__':
    test_read_l1()
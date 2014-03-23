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
def c1():
    c = calipso.Cal1(test_l1_file)
    return c


def test_date(c1):
    
    assert c1.year==2013
    assert c1.month == 4
    assert c1.day == 7
    

def load_vectors(c1, navg=1):
    
    lon, lat = c1.coords(navg=navg)
    vectors = [lon, lat]
    for func in c1.time, c1.utc_time, c1.datetimes, c1.surface_elevation:
        vectors.append(func(navg=navg))
    for func in c1.tropopause_height, c1.tropopause_temperature, c1.parallel_rms_baseline:
        vectors.append(func(navg=navg))
    return vectors

    
def test_nprof(c1):
    vectors = load_vectors(c1,navg=1)
    for vector in vectors:
        assert len(vector)==nprof
        
        
def test_nprof_avg(c1):
    vectors = load_vectors(c1, navg=30)
    for vector in vectors:
        assert len(vector)==np.floor(1.*nprof/30)
        
        
def test_nprof_too_many(c1):
    vectors = load_vectors(c1, navg=1000000)
    for vector in vectors:
        assert vector is None
        
        
def test_nprof_not_enough(c1):
    vectors = load_vectors(c1,navg=0)
    for vector in vectors:
        assert vector == []
        
    
def test_nz(c1):
    
    alt = c1.lidar_alt
    assert len(c1.lidar_alt)==nz
    assert len(c1.met_alt)==nz_met
    
        
def test_arrays_shape(c1, navg=30):
    
    for func in c1.atb, c1.perp, c1.atb_std, c1.atb1064:
        var = func(navg=navg)
        assert var.shape[0]==np.floor(1.*nprof/navg)
        assert var.shape[1]==nz
    
    for func in c1.mol,  c1.temperature, c1.rh, c1.pressure:
        var = func(navg=navg)
        assert var.shape[0]==np.floor(1.*nprof/navg)
        assert var.shape[1]==nz_met
            

if __name__=='__main__':
    test_read_l1()
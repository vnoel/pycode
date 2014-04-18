#!/usr/bin/env python
#encoding:utf-8

import pytest
import calipso_l2 as calipso
import numpy as np
from pyhdf.SD import HDF4Error


l2_file='/bdd/CALIPSO/Lidar_L2/05kmCLay.v3.02/2012/2012_01_01/CAL_LID_L2_05kmCLay-Prov-V3-02.2012-01-01T10-04-03ZN.hdf'
nonexistent_file='/really_wrong_path/one_more_time/this_file.hdf'
not_hdf_file='/users/noel/pycode/calipso.py'
folder_instead_of_file='/users/noel/'

nprof = 3728

@pytest.fixture(scope='module')
def c2():
    c = calipso.Cal2(l2_file)
    return c
    

def test_files():
    with pytest.raises(HDF4Error):
        c = calipso.Cal2(nonexistent_file)
    with pytest.raises(HDF4Error):
        c = calipso.Cal2(not_hdf_file)
    with pytest.raises(HDF4Error):
        c = calipso.Cal2(folder_instead_of_file)


def load_vectors(c2):
    lon, lat = c2.coords()
    vectors = [lon, lat]
    for func in c2.time, c2.utc_time, c2.datetime, c2.datetime2, c2.off_nadir_angle:
        vectors.append(func())
    for func in c2.tropopause_height, c2.tropopause_temperature, c2.dem_surface_elevation, c2.lidar_surface_elevation:
        vectors.append(func())
    return vectors


def test_vectors(c2):
    vectors = load_vectors(c2)
    for vector in vectors:
        assert np.size(vector, 0) == nprof
        
        
def load_layers(c2):
    nlay, base, top = c2.layers()
    layers = [base, top]
    for func in (c2.layers_pressure, c2.midlayer_temperature):
        layers.append(func())
    return nlay, layers
        
        
def test_layers(c2):
    nlay, layers = load_layers(c2)
    assert len(nlay)==nprof
    for layerinfo in layers:
        assert np.size(layerinfo,0)==nprof
        assert np.size(layerinfo,1)==10

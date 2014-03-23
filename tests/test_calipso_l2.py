#!/usr/bin/env python
#encoding:utf-8

import pytest
import calipso_l2 as calipso
import numpy as np


l2_file='/bdd/CALIPSO/Lidar_L2/05kmCLay.v3.02/2012/2012_01_01/CAL_LID_L2_05kmCLay-Prov-V3-02.2012-01-01T10-04-03ZN.hdf'

@pytest.fixture(scope='module')
def c2():
    c = calipso.Cal2(l2_file)
    return c
    

def load_vectors(c2):
    lon, lat = c2.coords()
    print np.size(lon, 0)
    assert np.size(lon, 0) == 12
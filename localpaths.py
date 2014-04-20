#!/usr/bin/env python
#encoding:utf-8

"""
Utility functions to find data files on Climserv or ICARE servers
Created by VNoel on 2014-02-27
"""

import socket

hostname = socket.gethostname()

icare_id = 'icare.univ-lille1.fr'
climserv_id = 'ipsl.polytechnique.fr'

# known data directories
l1dir = None        # CALIOP L1
l2dir = None        # CALIOP L2
eradir = None       # ECMWF

# host-dependent paths pointing to CALIOP level 1 and level 2 data
if hostname.endswith(icare_id):
    # ICARE
    l1dir = ('/DATA/LIENS/CALIOP/CAL_LID_L1.v3.01',)
    l2dir = ('/DATA/LIENS/CALIOP/05kmCLay.v3.02',
             '/DATA/LIENS/CALIOP/05kmCLay.v3.01')
    l2adir = ('/DATA/LIENS/CALIOP/05kmALay.v3.01',)
elif hostname.endswith(climserv_id):
    # CLIMSERV
    l1dir = ('/bdd/CALIPSO/Lidar_L1/CAL_LID_L1.v3.30',
             '/bdd/CALIPSO/Lidar_L1/CAL_LID_L1.v3.02',
             '/bdd/CALIPSO/Lidar_L1/CAL_LID_L1.v3.01')
    l2dir = ('/bdd/CALIPSO/Lidar_L2/05kmCLay.v3.02',
             '/bdd/CALIPSO/Lidar_L2/05kmCLay.v3.01',
             '/bdd/CALIPSO/Lidar_L2/05kmCLay.v2.01',
             '/bdd/CALIPSO/Lidar_L2/05kmCLay.v2.01')

# ECMWF data
if hostname.endswith(climserv_id):
    eradir = '/bdd/OPERA/NETCDF/GLOBAL_1125/4xdaily/'
    

def main():
    print ('CALIPSO L1 data dir : ',l1dir)
    print ('CALIPSO L2 data dir : ',l2dir)
    print ('ECMWF data dir : ', eradir)


if __name__ == '__main__':
    import plac
    plac.call(main)
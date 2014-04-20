#!/usr/bin/env python
#encoding:utf-8

"""
Utility functions to find data files on Climserv or ICARE servers
Created by VNoel on 2014-02-27
"""

import socket

hostname = socket.gethostname()

# host-dependent paths pointing to CALIOP level 1 and level 2 data
if hostname.endswith('icare.univ-lille1.fr'):
    # ICARE
    l1dir = ('/DATA/LIENS/CALIOP/CAL_LID_L1.v3.01',)
    l2dir = ('/DATA/LIENS/CALIOP/05kmCLay.v3.02',
             '/DATA/LIENS/CALIOP/05kmCLay.v3.01')
    l2adir = ('/DATA/LIENS/CALIOP/05kmALay.v3.01',)
elif hostname.endswith('ipsl.polytechnique.fr'):
    # CLIMSERV
    l1dir = ('/bdd/CALIPSO/Lidar_L1/CAL_LID_L1.v3.30',
             '/bdd/CALIPSO/Lidar_L1/CAL_LID_L1.v3.02',
             '/bdd/CALIPSO/Lidar_L1/CAL_LID_L1.v3.01')
    l2dir = ('/bdd/CALIPSO/Lidar_L2/05kmCLay.v3.02',
             '/bdd/CALIPSO/Lidar_L2/05kmCLay.v3.01',
             '/bdd/CALIPSO/Lidar_L2/05kmCLay.v2.01',
             '/bdd/CALIPSO/Lidar_L2/05kmCLay.v2.01')
else:
    l1dir = None
    l2dir = None


def main():
    print ('CALIPSO L1 data dir : ',l1dir)
    print ('CALIPSO L2 data dir : ',l2dir)


if __name__ == '__main__':
    import plac
    plac.call(main)
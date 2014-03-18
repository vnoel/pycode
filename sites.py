#!/usr/bin/env python
#encoding:utf-8

"""
Utility functions to find CALIOP files on Climserv or ICARE servers
Created by VNoel on 2014-02-27
"""

import socket
import glob


# server-dependent paths pointing to CALIOP level 1 and level 2 data

hostname = socket.gethostname()
if hostname == 'access.icare.univ-lille1.fr':
    l1dir = ('/DATA/LIENS/CALIOP/CAL_LID_L1.v3.01', None)
    l2dir = ('/DATA/LIENS/CALIOP/05kmCLay.v3.02',
             '/DATA/LIENS/CALIOP/05kmCLay.v3.01')
    l2adir = ('/DATA/LIENS/CALIOP/05kmALay.v3.01',)
else:
    l1dir = ('/bdd/CALIPSO/Lidar_L1/CAL_LID_L1.v3.30/',
             '/bdd/CALIPSO/Lidar_L1/CAL_LID_L1.v3.02/',
             '/bdd/CALIPSO/Lidar_L1/CAL_LID_L1.v3.01/')
    l15dir = ('/bdd/CFMIP/CAL_LID_L1.5',)
    # for level 2, look first in our own data, then on climserv data
    l2dir = ('/bdd/CFMIP/Lidar_L2',
             '/bdd/CALIPSO/Lidar_L2/05kmCLay.v3.02',
             '/bdd/CALIPSO/Lidar_L2/05kmCLay.v3.01',
             '/bdd/CALIPSO/Lidar_L2/05kmCLay.v2.01',
             '/bdd/CALIPSO/Lidar_L2/05kmCLay.v2.01')


def l1_files(y, m, d, mask='*'):
    """
    returns the list of available CALIOP L1 files matching a date and mask
    """
    goodpath = '/'
    for l1d in l1dir:
        calpath = l1d + '/%04d/%04d_%02d_%02d/' % (y, y, m, d)
        if os.path.isdir(calpath):
            goodpath = calpath
    files = glob.glob(goodpath + mask + '.hdf')
    return files


def l1_night_files(y, m, d):
    """
    returns the list of available CALIOP nighttime L1 files matching a date
    """
    files = l1_files(y, m, d, '*ZN*')
    return files


def l1_file_from_orbit(y, m, d, orbit):
    """returns full path to a CALIOP L1 orbit file"""
    files = l1_files(y, m, d, '*' + orbit + '*')
    if not files:
        return None
    return files[0]


def l1_file_from_l2_file(y, m, d, l2file):
    orbit = l2file[-15:-4]
    l1file = l1_file_from_orbit(y, m, d, orbit)
    return l1file


def l15_files(y, m, d, mask):
    calpath = l15dir[0] + '/%04d/%04d_%02d_%02d/' % (y, y, m, d)
    files = glob.glob(calpath + mask + '.nc')
    return files


def l15_night_files(y, m, d):
    files = l15_files(y, m, d, '*ZN*')
    return files


def l15_file_from_orbit(y, m, d, orbit):
    files = l15_files(y, m, d, '*' + orbit + '*')
    if not files:
        return None
    return files[0]


def l2_files(y, m, d, mask):
    datepath = '/%04d/%04d_%02d_%02d/' % (y, y, m, d)
    for l2d in l2dir:
        calpath = l2d + datepath
        if os.path.isdir(calpath):
            files = glob.glob(calpath + mask + '.hdf')
            return files
    return None


def l2_cfiles(y, m, d, mask):
    return l2_files(y, m, d, mask)


def l2_night_files(y, m, d):
    files = l2_files(y, m, d, '*ZN*')
    return files


def l2_day_files(y, m, d):
    files = l2_files(y, m, d, '*ZD*')
    return files


def l2_file_from_orbit(y, m, d, orbit):
    files = l2_files(y, m, d, '*' + orbit + '*')
    if not files:
        return None
    return files[0]


def l2_afiles(y, m, d, mask):
    datepath = '/%04d/%04d_%02d_%02d/' % (y, y, m, d)
    calpath = l2adir[0] + datepath
    files = glob.glob(calpath + mask + '.hdf')
    if not files:
        calpath = l2adir[1] + datepath
        files = glob.glob(calpath + mask + '.hdf')
    return files


def l2bis_files(y, m, d, mask):
    datepath = '/%04d/%04d_%02d_%02d/' % (y, y, m, d)
    calpath = '/bdd/CFMIP/CAL_LID_L2bis' + datepath
    files = glob.glob(calpath + mask + '.nc')
    return files


def l2bis_night_files(y, m, d):
    return l2bis_files(y, m, d, '*ZN*')


def l2bis_file_from_orbit(y, m, d, orbit):
    files = l2bis_files(y, m, d, '*' + orbit + '*')
    if not files:
        return None
    return files[0]


def read_l2bis(calfile):
    from scipy.io.matlab import mio

    try:
        x = mio.loadmat(calfile)
    except:
        print 'Cannot read file ' + calfile
        return [None] * 5
    t = x['T_lay']
    lat = x['lat'][0, :]
    lon = x['lon'][0, :]
    zbase = x['zbase']
    ztop = x['ztop']
    return lon, lat, zbase, ztop, t

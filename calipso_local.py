#!/usr/bin/env python
#encoding:utf-8

# Created by VNoel on 2014-04-20

'''
Utility functions to list available CALIOP files on local host
'''

import os
import glob
from localpaths import l1dir, l2dir, l2dir_333

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


def l2_files(y, m, d, mask='*', havg=5):
    if havg==0.333:
        l2dir_list = l2dir_333
    else:
        l2dir_list = l2dir
    goodpath = '/'
    datepath = '/%04d/%04d_%02d_%02d/' % (y, y, m, d)
    for l2d in l2dir_list:
        calpath = l2d + datepath
        if os.path.isdir(calpath):
            goodpath = calpath
    fullmask = goodpath + mask + '.hdf'
    print 'Looking in ' + fullmask
    files = glob.glob(fullmask)
    return files
    

def l2_cfiles(y, m, d, mask):
    return l2_files(y, m, d, mask)


def l2_night_files(y, m, d, havg=5):
    files = l2_files(y, m, d, '*ZN*', havg=havg)
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


def main(year=2007, month=5, day=4):
    
    print('Searching CALIOP files for %04d/%02d/%02d/' % (year, month, day))
    l1files = l1_files(year, month, day)
    l2files = l2_files(year, month, day)
    print('Found %d L1 files' % len(l1files))
    print('Found %d L2 files' % len(l2files))
    print('Available CALIOP L1 files : ', l1files)
    print('Available CALIOP l2 files : ', l2files)


if __name__ == '__main__':
    import plac
    plac.call(main)
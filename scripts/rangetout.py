#!/usr/bin/env python

import os
import glob
import shutil

srcdir='./'
destdir='./'

year = 2010
month = 4

print('On range %04d-%02d' % (year, month))

days = list(range(32))

'''
for day in days:
    newdir = destdir+('%04d_%02d_%02d' % (year,month,day))
    print 'Creating '+newdir
    os.mkdir(newdir)
    '''

mask = '*%04d-%02d-*.hdf' % (year, month)
files = glob.glob(srcdir+mask)
for file in files:
    day = int(file[-17:-15])
    daydir = destdir+'%04d_%02d_%02d' % (year, month, day)
    print('Moving '+file+' in '+daydir)
    shutil.move(file, daydir)

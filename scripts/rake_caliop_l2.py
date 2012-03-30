#!/usr/bin/env python
# encoding: utf-8
"""
rake_caliop.py

finds CALIOP files in a required time frame
links them in a folder

Created by Vincent Noel - LMD/CNRS on 2012-01-25.
"""

from datetime import datetime, timedelta
import calipso
import os

def rake_caliop_l2(y='2008', m='8', target='out'):
    
    if not os.path.isdir(target):
        exit('Error : ' + target + ' is not a path')
    
    year = int(y)
    month = int(m)
    start = datetime(year, month, 1)
    end = datetime(year, month+1, 1)
    
    current = start
    
    n = 0
    while current <= end:
        files = calipso.l2_files(current.year, current.month, current.day, '*')
        if files:
            n+=len(files)
            for f in files:
                try:
                    linkedfile = target + '/' + os.path.basename(f)
                    os.symlink(f, linkedfile)
                except OSError:
                    # catch already-linked-error silently
                    pass
        current += timedelta(days=1)
    
    print n, 'files linked in ', target
    return


if __name__ == '__main__':
    import plac
    plac.call(rake_caliop_l2)


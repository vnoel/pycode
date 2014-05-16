#!/usr/bin/env python
#encoding:utf-8

"""
draft module to read NASA AMES data files.
This is NOT FINISHED.
It won't even load.
"""

import numpy as np

def ames_read(f):

    fid = open(f, 'r')
    lines = fid.readlines(7)
    fid.close()
    
    header = lines[0:7]
    nlhead = int(header[0].split()[0])
    print 'N lines : ', nlhead
    data = lines[nlhead:-1]
    
    nmes = 0
    imes = 0
    file_end = False
    while not file_end:
        d = datetime(data[imes])
        npts = data[imes+1].split()[0]
        for i in np.r_[5:5+npts]:
            
        
    # this is bullshit
        
    
    
def main(f):
    
    ames_read(f)
    
if __name__ == '__main__':
    import plac
    plac.call(main)
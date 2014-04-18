# Fonctions utiles pour l'analyse de donnees MODIS
# VN 2008
# LMDX

from pyhdf.SD import SD
import glob
import numpy as np

modisbasedir = '/bdd/CFMIP/ATRAIN_COLOC/MODIS_COLOC/'

def planck(wl, radiance):
    '''Renvoie la temperature de brillance pour des radiances, en connaissant la longueur
    d'onde wl en microns'''
    c1 = 14390.
    c2 = 1.191e8
    temp = c1 / (wl * np.log(1. + c2/(radiance * (wl**5))))
    return temp 
    

def modispath (y, m, d):
    return modisbasedir + '/%04d/%04d_%02d_%02d/' % (y, y, m, d)


def orbit(mfile):
    return mfile[-15:-4]

def find_file_for_orbit (y, m, d, orbit):
    p = modispath (y, m, d)
    files = glob.glob(p + '*' + orbit + '*.hdf')
    if not files:
        print("No MODIS data found")
        return None

    if len(files)>1:
        print("More than 1 MODIS file has been found : ", files)

    return files[0]


def nightfiles (y, m, d):
    p = modispath (y, m, d)
    files = glob.glob(p + '*ZN.hdf')
    return files


## Non-layer info

def t11(mfile):
    hdffile = SD(mfile)
    x = hdffile.select('MYD021KM_EV_1KM_Emissive_Band31')
    scale_factor = x.attributes()['scale_factor']
    add_offset = x.attributes()['add_offset']
    x = x[:]
    hdffile.end()
    
    idx = x > 65534
    x = (x - add_offset) * scale_factor
    t = planck (11.03, x)
    t[idx] = np.nan
    return t


def t12(mfile):
    hdffile = SD(mfile)
    x = hdffile.select('MYD021KM_EV_1KM_Emissive_Band32')
    scale_factor = x.attributes()['scale_factor']
    add_offset = x.attributes()['add_offset']
    x = x[:]
    hdffile.end()    
    
    idx = x > 65534
    x = (x - add_offset) * scale_factor
    t = planck (12.02, x)
    t[idx] = np.nan
    return t


def latitude(mfile):
    hdffile = SD(mfile)
    x = hdffile.select('Latitude')[:]
    hdffile.end()
    return x
    
def longitude(mfile):
    hdffile = SD(mfile)
    x = hdffile.select('Longitude')[:]
    hdffile.end()
    return x

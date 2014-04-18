#!/usr/bin/env python
# encoding: utf-8

import numpy as np
import netCDF4

# for backward compatibility
# pre-class

def ecmwf_file (year, month, prefix):
    p = '/bdd/OPERA/NETCDF/GLOBAL_1125/4xdaily/AN_PL/' + str(year) + '/'
    name = p + prefix + ('.%04d%02d.aph.GLOBAL_1125.nc' % (year, month))
    return name

def read_named_var(year, month, varname,level=-1, lon360=None):
    f = ecmwf_file (year, month, varname)
    print(('Reading '+f))
    nc = netCDF4.Dataset(f)
    
    lon = nc.variables['lon'][:]
    lat = nc.variables['lat'][:]
    levels = nc.variables['level'][:]
    
    if level != -1:
        levels = levels[level]
    
    # passer de 0->360 a -180->180 (comme ds fichiers calipso)
    if not lon360:
        lon[lon >= 180] -= 360.

    add_offset = 0.
    scale_factor = 1.
    
    # getattr ne marche pas avec les structures netCDF, pfff
    if hasattr(nc.variables[varname], 'add_offset'):
        add_offset = nc.variables[varname].add_offset
        
    if hasattr(nc.variables[varname], 'scale_factor'):
        scale_factor = nc.variables[varname].scale_factor
        
    if level==-1:
        data = (nc.variables[varname][:] * scale_factor) + add_offset
    else:
        data = (nc.variables[varname][:,level,:,:] * scale_factor) + add_offset
    
    nc.close()
    
    return lon, lat, levels, data
    
def ecmwf_surface_file(year, month, prefix):
    p = '/bdd/OPERA/NETCDF/GLOBAL_1125/4xdaily/AN_SF/' + str(year) + '/'
    name = p + prefix + ('.%04d%02d.ash.GLOBAL_1125.nc' % (year, month))
    return name
    
def read_surface_var(year, month, varname):
    f = ecmwf_surface_file(year, month, varname)
    print(('Reading '+f))
    nc = netCDF4.Dataset(f)
    
    lon = nc.variables['lon'][:]
    lat = nc.variables['lat'][:]
    
    # passer aux coord calipso
    lon[lon >= 180] -= 360.
    
    x = nc.variables[varname]
    
    add_offset = 0.
    scale_factor = 1.
    
    if hasattr(x, 'add_offset'):
        add_offset = x.add_offset
        
    if hasattr(x, 'scale_factor'):
        scale_factor = x.scale_factor
        
    data = (x[:] * scale_factor) + add_offset
    nc.close()
    
    return lon, lat, data

#bla bla bla

# la classe ecmwf ne correspond pas a un fichier
# mais a un des sous-ensembles ecmwf dispos sur climserv

class Ecmwf:

    # pour initialiser la classe, il faut passer
    # une chaine de caractere contenant la resolution
    # voulue pour les fichiers ECMWF
    # e.g. '075' ou '1125'
    def __init__(self, deg_res):
        if deg_res=='1125':
            self.path = '/bdd/OPERA/NETCDF/GLOBAL_1125/4xdaily/'
            self.identifier = ''
            self.desc = 'OPERA Global dataset, 1.125° res'
            res = 1.125
        elif deg_res=='075':
            self.path = '/bdd/ERAI/NETCDF/GLOBAL_075/4xdaily/'
            self.identifier = 'ei'
            self.desc = 'ERA-Interim Global dataset, 0.75° res'
            res = 0.75
        else:
            raise NameError('Unknown resolution: ' + deg_res)
        self.res = deg_res
        self.lat = np.r_[90.:-90.-res:-res]
        self.lon = np.r_[0:360:res]
        
            
    def pl_file(self, year, month, varname):
        # aph for era40 pressure level files
        # aphei for era-interim pressure level files
        identifier = 'aph' + self.identifier
        filename = varname + '.%04d%02d.%s.GLOBAL_%s.nc' % (year, month, identifier, self.res)
        path = self.path + 'AN_PL/' + '%04d' % year + '/'
        name = path + filename
        return name

    def pl_var(self, year, month, varname,level=-1, lon180=None):
        f = self.pl_file (year, month, varname)
        print(('Reading '+f))
        try:
            nc = netCDF4.Dataset(f)
        except:
            raise NameError('No such ECMWF file: ' + f)

        lon = nc.variables['lon'][:]
        lat = nc.variables['lat'][:]
        if not(all(self.lon == lon) and all(lat==self.lat)):
            raise NameError('Longitude or latitude vectors not of expected size')
        
        levels = nc.variables['level'][:]

        if level > -1:
            levels = levels[level]

        # passer de 0->360 a -180->180 (comme ds fichiers calipso)
        if lon180:
            lon[lon >= 180] -= 360.

        add_offset = 0.
        scale_factor = 1.

        # getattr ne marche pas avec les structures netCDF, pfff
        if hasattr(nc.variables[varname], 'add_offset'):
            add_offset = nc.variables[varname].add_offset

        if hasattr(nc.variables[varname], 'scale_factor'):
            scale_factor = nc.variables[varname].scale_factor

        if level==-1:
            data = (nc.variables[varname][:] * scale_factor) + add_offset
        else:
            data = (nc.variables[varname][:,level,:,:] * scale_factor) + add_offset

        nc.close()
        return lon, lat, levels, data
    
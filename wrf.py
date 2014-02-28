#!/usr/bin/env python
#encoding: utf-8

'''
module wrf 
V. Noel 2010
'''

import netCDF4
import numpy as np
import matplotlib.pyplot as plt

# private functions

def _tk (p, tpot):
    f = 8.314472/29.19
    tk = (tpot + 300)/((1000/(p/100))**f)
    return tk

def _remap (var, lon, lat, lonorb, latorb):
    '''Extrait des profils verticaux WRF sur une liste de coordonnees lon/lat.'''

    n = lonorb.size                   # nombre de coordonnees
    var2 = np.zeros([var.shape[0], n])

    for i in np.arange(n):
        # pour chaque profil, on calcule une grille de distance avec les coords WRF
        d = np.sqrt((lon-lonorb[i])**2 + (lat-latorb[i])**2)
        #im = argmin(d)
        #imin, jmin = unravel_index(im, lon.shape)
        # the code below is much faster
        imin, jmin = divmod(d.argmin(), d.shape[1])
        var2[:,i] = var[:,imin,jmin]

    return var2.T

class wrf(object):
    
    def __init__(self, wrfout):
        self.nc = netCDF4.Dataset(wrfout)
    
    def close(self):
        self.nc.close()
        self.nc = None        
    
    def coords(self, it=0):
        if it >= self.nc.variables['XLONG'].shape[0]:
            print '*ERROR* : Requested timestamp not in WRF file'
            print 'Requested timestamp index : ', it
            print 'Number of time indexes in WRF file', self.nc.variables['XLONG'].shape[0]
            return None, None
        lon = self.nc.variables['XLONG'][it,...]
        lat = self.nc.variables['XLAT'][it,...]
        lon = np.squeeze(lon)
        lat = np.squeeze(lat)
        return lon, lat
        
    def ntime(self):
        return self.nc.variables['Times'].shape[0]

    def times(self):
        t = self.nc.variables['Times']
        ts = []
        ts = [ttxt.tostring() for ttxt in t[:,:]]
        return t

    def time(self, it=0):
        t = self.nc.variables['Times'][it,:]
        t = t.tostring()
        return t
        
    def p_top(self, it=0):
        p = self.nc.variables['P_TOP'][it]
        return p
    
    def pressure(self, it=0, pascals=False, on_orbit=None):
        '''
        Lit le champ de pression du fichier WRF. 
        Parametres:  
        it: indice temporel du champ a extraire [default: 0]
        pascals: la pression sera renvoyee en pascals si True, en hPa si False.
        on_orbit: Liste facultative de lon et lat sur lesquels extraire les profils (lon_orbit, lat_orbit)
        Renvoie:
        champ de pression 3D [X, Y, Z] ([PROFIL, Z] si on_orbit)
        '''
        p = self.nc.variables['P'][it,...] + self.nc.variables['PB'][it,...]
        if not pascals:
            p /= 100.
        if on_orbit:
            wrflon, wrflat = self.coords(it=it)
            p = _remap(p, wrflon, wrflat, on_orbit[0], on_orbit[1])
        p = np.squeeze(p)
        return p

    def temperature(self, it=0, kelvins=False, on_orbit=None):
        '''
        Lit le champ de temperature du fichier WRF. 
        Parametres:  
        it: indice temporel du champ a extraire [default: 0]
        kelvins: la temperature sera renvoyee en K si True, en Celsius si False.
        on_orbit: Liste facultative de lon et lat sur lesquels extraire les profils (lon_orbit, lat_orbit)
        Renvoie:
        champ de temperature 3D [X, Y, Z] ([PROFIL, Z] si on_orbit)
        '''
        
        p = self.nc.variables['P'][it,...] + self.nc.variables['PB'][it,...]
        tpot = self.nc.variables['T'][it,...]
        t = _tk (p, tpot)
        if not kelvins:
            t -= 273.
        if on_orbit:
            wrflon, wrflat = self.coords(it=it)
            t = _remap(t, wrflon, wrflat, on_orbit[0], on_orbit[1])
            
        t = np.squeeze(t)
        return t


    def u(self, it=0):
        u = self.nc.variables['U'][it,...]
        u = np.squeeze(u)
        return u
        
        
    def v(self, it=0):
        v = self.nc.variables['V'][it,...]
        v = np.squeeze(v)
        return v
    

    def w(self, it=0, on_orbit=None):
        '''
        Lit le champ de vitesse de vent vertical du fichier WRF. 
        Parametres:  
        it: indice temporel du champ a extraire [default: 0]
        on_orbit: Liste facultative de lon et lat sur lesquels extraire les profils (lon_orbit, lat_orbit)
        Renvoie:
        champ de vitesse de vent verticale 3D [X, Y, Z] ([PROFIL, Z] si on_orbit)
        '''
        w = self.nc.variables['W'][it,...]
        if on_orbit:
            wrflon, wrflat = self.coords(it=it)
            w = _remap(w, wrflon, wrflat, on_orbit[0], on_orbit[1])
        return w

def map_peninsula (lon, lat, xvar, centerlon=-60, w=2000,h=2000, cb=False, cl=None, cbtitle=None, ec='None', ax=None, fp=None):
    from mpl_toolkits import basemap
    from matplotlib.pyplot import colorbar
    import niceplots
    
    map = basemap.Basemap(projection='lcc',lat_0=-70,lon_0=centerlon,width=w*1000,height=h*1000, resolution='i', ax=ax)
    x, y = map(lon, lat)
    if cl is None:
        pc = map.pcolormesh(x, y, xvar, edgecolors=ec)
    else:
        pc = map.pcolormesh(x, y, xvar, vmin=cl[0], vmax=cl[1], edgecolors=ec)
    map.drawmapboundary(color='grey')
    if fp:
        map.drawmeridians(np.arange(-105,-15,15), labels=[0,0,0,1], fontproperties=fp) # left, right, top, bottom
        map.drawparallels(np.arange(-90,90,10), labels=[1,0,0,0], fontproperties=fp) # left, right, top, bottom
    else:
        map.drawmeridians(np.arange(-105,-15,15), labels=[0,0,0,1]) # left, right, top, bottom
        map.drawparallels(np.arange(-90,90,10), labels=[1,0,0,0]) # left, right, top, bottom
    map.drawcoastlines(color='grey')
    
    if cb:
        cbar=colorbar(pc)
        if fp:
            niceplots.beautify_colorbar(cbar, title=cbtitle)
        else:
            if cbtitle is not None:
                cbar.set_label(cbtitle)
    
    return map, pc

def antarctica_map(maxlat=-50, lon_0=180):
    from mpl_toolkits.basemap import Basemap
    
    m = Basemap(projection='splaea', boundinglat=maxlat, lon_0=lon_0, resolution='i')
    
    return m

# def map_antarctica_temp(wrfout, maxlat=-50, it=0):
#     from mpl_toolkits.basemap import Basemap
#     
#     m = Basemap(projection='splaea', boundinglat=maxlat, lon_0=180)
#     
#     lon, lat = coords(wrfout)
#     temp = temperature(wrfout)
#     plt.figure()
#     x, y = m(lon, lat)
#     m.pcolormesh(x, y, temp[it,:,:])
#     m.drawcoastlines()
#     plt.colorbar()

# def geo_show(f='geo_em.d01.nc'):
#     print 'Reading file ' + f
#     try:
#         nc = netCDF4.Dataset(f)
#     except:
#         raise NameError(f + ' is not a netCDF file')
#     
#     try:
#         lat = nc.variables['XLAT_M'][0,:,:]
#         lon = nc.variables['XLONG_M'][0,:,:]
#         hgt = nc.variables['HGT_M'][0,:,:]
#     except:
#         raise NameError('Cannot find variables XLAT_M, XLONG_M, HGT_M in ' + f)
#     
#     #m = Basemap(projection='lcc', lat_0=np.min(lat.ravel()), lon_0=np.min(lon.ravel()), width=10000*1000, height=6000*1000)
#     m = Basemap(projection='robin', lon_0=0)
#     x, y = m(lon, lat)
#     
#     print 'Plotting the damn thing'
#     plt.figure()
#     
#     m.pcolor(x, y, hgt)
#         
#     m.drawmapboundary()
#     m.drawmeridians(np.r_[-180:180:30])
#     m.drawparallels(np.r_[-90:90:15])
#     
#     plt.show()

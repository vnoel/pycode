# Fonctions utiles pour l'analyse de donnees CALIOP
# VN 2008-2010
# LMDX/CNRS

from pyhdf.SD import SD
from scipy.stats import nanmean
from scipy.integrate import trapz
import numpy as np
import numpy.ma as ma
import glob
import datetime
import os
import socket
import warnings

# cannot use pytables since calipso is hdf4...

datapath = os.path.dirname(__file__) + '/staticdata/'
lidar_alt = np.loadtxt(datapath+'lidaralt.asc')
met_alt = np.loadtxt(datapath+'metalt.asc')


# the following ranges have been calibrated more carefully for nighttime than daytime

# maximum atb to consider for calibration in the 26-28 km range
atb_max = {'ZN':1e-4, 'ZD':1}

# maximum integrated atb for calibrationi in the 26-28 km range
iatb_bounds = {'ZN':[1e-5, 8e-5], 'ZD':[-8e-3, 8e-3]}


# server-dependent paths
hostname = socket.gethostname()
if hostname=='access.icare.univ-lille1.fr':
    l1dir = ('/DATA/LIENS/CALIOP/CAL_LID_L1.v3.01',)
    l2dir = ('/DATA/LIENS/CALIOP/05kmCLay.v3.01', '/DATA/LIENS/CALIOP/05kmCLay.v3.01')
    l2adir = ('/DATA/LIENS/CALIOP/05kmALay.v3.01',)
else:
    # on climserv, look first in the subset on v3.01 data available, then use the 2.01
    l1dir = ('/bdd/CALIPSO/Lidar_L1/CAL_LID_L1.v3.01/','/bdd/CALIPSO/Lidar_L1/CAL_LID_L1.v2.01')
    l15dir = ('/bdd/CFMIP/CAL_LID_L1.5',)
    # for level 2, look first in our own data, then on climserv data
    l2dir = ('/bdd/CFMIP/Lidar_L2', '/bdd/CALIPSO/Lidar_L2/05kmCLay.v2.01')

def l1_files (y, m, d, mask):
    goodpath = '/'
    for l1d in l1dir:
        calpath = l1d + '/%04d/%04d_%02d_%02d/' % (y, y, m, d)
        if os.path.isdir(calpath):
            goodpath = calpath
    files = glob.glob(goodpath + mask + '.hdf')
    return files

def l1_night_files (y, m, d):
    files = l1_files(y, m, d, '*ZN*')
    return files

def l1_file_from_orbit (y, m, d, orbit):
    files = l1_files(y, m, d, '*'+orbit+'*')
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

def l15_night_files (y, m, d):
    files = l15_files(y, m, d, '*ZN*')
    return files

def l15_file_from_orbit (y, m, d, orbit):
    files = l15_files(y, m, d, '*'+orbit+'*')
    if not files:
        return None
    return files[0]

def l2_files(y, m, d, mask):
    goodpath = '/'
    datepath = '/%04d/%04d_%02d_%02d/' % (y, y, m, d)
    for l2d in l2dir:
        calpath = l2d + datepath
        if os.path.isdir(calpath):
            files = glob.glob(calpath + mask + '.hdf')
            return files
    return None

def l2_cfiles(y, m, d, mask):
    return l2_files(y, m, d, mask)

def l2_night_files (y, m, d):
    files = l2_files(y, m, d, '*ZN*')
    return files

def l2_file_from_orbit(y, m, d, orbit):
    files = l2_files(y, m, d, '*'+orbit+'*')
    if not files: return None
    return files[0]

def l2_afiles(y, m, d, mask):
    datepath = '/%04d/%04d_%02d_%02d/' % (y, y, m, d)
    calpath = l2adir[0] + datepath
    files = glob.glob(calpath + mask + '.hdf')
    if not files:
        calpath = l2dirbis + datepath
        files = glob.glob(calpath + mask + '.hdf')
    return files

def l2bis_files(y, m, d, mask):
    datepath = '/%04d/%04d_%02d_%02d/' % (y, y, m, d)
    calpath = '/bdd/CFMIP/CAL_LID_L2bis' + datepath
    files = glob.glob(calpath + mask + '.nc')
    return files

def l2bis_night_files (y, m, d):
    return l2bis_files(y, m, d, '*ZN*')
    
def l2bis_file_from_orbit(y, m, d, orbit):
    files = l2bis_files(y, m, d, '*'+orbit+'*')
    if not files: return None
    return files[0]

def read_l2bis (calfile):
    from scipy.io.matlab import mio
    try:
        x = mio.loadmat(calfile)
    except:
        print 'Cannot read file '+calfile
        return [None] * 5
    t = x['T_lay']
    lat = x['lat'][0,:]
    lon = x['lon'][0,:]
    zbase = x['zbase']
    ztop = x['ztop']
    return lon, lat, zbase, ztop, t


# Fonctions maths utiles

def _integrate_signal(data, alt):
    integrale = trapz(data[::-1], x=alt[::-1])
    return integrale

def _vector_average (v0, navg):
    ''' v = _vector_average (v0, navg)
        moyenne le vector v0 tous les navg points.'''

    v0 = v0.squeeze()

    assert v0.ndim == 1, 'in _vector_average, v0 should be a vector'
    if navg == 1: return v0
        
    n = np.size(v0,0)/navg
    v = np.zeros(n)
    for i in np.arange(n):
        n0 = i * navg
        v[i] = v0[n0:n0+navg-1].mean()
    return v

def _array_average (a0, navg, weighted = False):
    ''' a = _array_average (a0, navg, weighted=False)
        moyenne le tableau a0 le long des x tous les navg profils.'''

    #a0 = a0.squeeze()

    assert a0.ndim == 2, 'in _array_average, a0 should be a 2d array'
    if navg == 1: return a0

    if weighted and (navg % 2 != 1):
        weighted = False
        print "_array_average: navg is even, turning weights off"

    # create triangle-shaped weights
    if weighted:
        w = np.zeros(navg)
        w[:navg/2.+1] = np.r_[1:navg/2.+1]
        w[navg/2.+1:] = np.r_[int(navg/2.):0:-1]

    # create averaged array a with number of averaged profiles n
    n = np.size(a0,0)/navg
    a = np.zeros([n, np.size(a0,1)])
    
    for i in np.arange(n):
        n0 = i * navg
        if weighted:
            a[i,:] = np.average(a0[n0:n0+navg,:], axis=0, weights=w)
        else:
            a[i,:] = a0[n0:n0+navg-1,:].mean(axis=0)
    return a

def _remap_y (z0, y0, y):
    ''' z = remap (z0, y0, y)
            interpole les donnees du tableau z0 sur un nouveau y.'''
    z = np.zeros([np.size(z0,0), np.size(y)])
    for i in np.arange(np.size(z,0)):
        z[i,:] = np.interp(y, np.flipud(y0), np.flipud(z0[i,:]))
    return z

def ncep_remap (var, lat, lon, latorb, lonorb):
    n = np.size(latorb,0)
    var2 = np.zeros([n, np.size(var,0)])

    for i in np.arange(n):
        d = np.sqrt((lat-latorb[i])**2)
        xlat = d.argmin()
        d = np.sqrt((lon-lonorb[i])**2)
        xlon = d.argmin()
        var2[i,:] = var[:,xlat,xlon]

    return var2


# Generic calipso file class
class _Cal:
        
    # il faut faire un try de l'instantiation de Cal1 ou Cal2
    # car hdf.SD balance une exception si le fichier n'existe pas.
    # la bonne methode en python n'est pas d'empecher l'instantiation
    # mais d'attrapper l'exception levee durant l'instantiation
    def __init__(self, filename):
        warnings.simplefilter('ignore', DeprecationWarning)
        self.hdf = SD(filename)
        self.filename = filename
        self.orbit = filename[-15:-4]
        self.z = self.orbit[-2:]            # zn or zd
        self.year = int(filename[-25:-21])
        self.month = int(filename[-20:-18])
        self.day = int(filename[-17:-15])
        hour = int(filename[-14:-12])
        minutes = int(filename[-11:-9])
        seconds = int(filename[-8:-6])
        self.date = datetime.datetime(self.year, self.month, self.day, hour, minutes, seconds)
        
    def __repr__(self):
        return self.filename
        
    def close(self):
        self.hdf.end()
        self.hdf = None


# classe pour level 1
class Cal1(_Cal):

    def _read_var(self, var, navg, idx=(0,-1)):
        '''
        Read a variable in an hdf file, averaging the data if navg is > 1
        '''
        try:
            var = self.hdf.select(var)
        except:
            print var+' not present in file '+self.filename+' (version too old, maybe ?)'
            return None
            
        data = var[:].squeeze()
        if data.ndim == 1:
            data = var[idx[0]:idx[1]]
            if navg > 1: 
                data = _vector_average(data, navg)
        else:
            data = var[idx[0]:idx[1],:]
            if navg  > 1: 
                data = _array_average(data, navg)

        var.endaccess()

        return data
        
    def time(self, navg=30):
        var = self.hdf.select('Profile_Time')
        time  = var[0:-1].squeeze()
        
        if navg > 1:
            n = np.size(time,0)/navg
            time2 = np.zeros(n)
            for i in np.arange(n):
                time2[i] = time[i * navg]
            time = time2
            
        return time

    # Coordoonees
    def coords (self, navg=30, idx=(0,-1)):
        ''' lon, lat = read_coords (filename, navg)
                lit les coordonnees profil calipso, moyennes sur navg'''
        lat = self._read_var('Latitude', navg, idx=idx)
        lon = self._read_var('Longitude', navg, idx=idx)

        return lon, lat

    def find_profiles_in_lonlat_box (self, lonmin, lonmax, latmin, latmax):
        ''' idx = find_profiles_in_lonlat_box (self, lonmin, lonmax, latmin, latmax)
            Trouve les index min et max correspondant a une boite lon, lat dans un fichier CALIPSO.
            Renvoie None si aucun profil n'est trouve.
            Devrait marcher independemment avec les fichiers niveau 1 ou niveau 2.
            '''
        lon = hdf.select('Longitude').get()
        mask1 = (lon >= lonmin) & (lon < lonmax)
        if mask1.sum()==0:
            return None

        lat = hdf.select('Latitude').get()
        mask2 = (lat >= latmin) & (lat < latmax)
        idx = np.where(mask1 & mask2)[0]

        if np.size(idx)==0:
            return None
        else:
            return idx[0], idx[-1]

    def find_profiles_south_of (self, nlimit):
        lat = hdf.select('Latitude').get()
        ind = np.where(lat < nlimit)[0]
        if np.size(ind)==0:
            return None
        else:
            return ind[0], ind[-1]

    # Data
    
    def surface_elevation (self, navg=30, prof=None, idx=(0,-1)):
        '''Lit elevation'''
        elev = self._read_var('Surface_Elevation', navg, idx=idx)
        if prof is not None:
            elev = elev[prof,:]
        return elev

    def perp (self, navg=30, prof=None, idx=(0,-1)):
        '''perp = read_perp (filename, navg)'''
        perp = self._read_var('Perpendicular_Attenuated_Backscatter_532', navg, idx=idx)
        if prof:
            perp = perp[prof,:]

        return perp

    def atb (self, navg=30, prof=None, idx=(0,-1)):
        ''' atb = cal_read_atb (filename, navg)
                lit les donnees calipso latitude, longitude et atb'''
        atb = self._read_var('Total_Attenuated_Backscatter_532', navg, idx=idx)
        if prof:
            atb = atb[prof,:]
        return atb

    def atb1064 (self, navg=30, prof=None, idx=(0,-1)):
        ''' atb = cal_read_atb (filename, navg)
                lit les donnees calipso latitude, longitude et atb'''
        atb = self._read_var('Attenuated_Backscatter_1064', navg, idx=idx)
        if prof:
            atb = atb[prof,:]
        return atb

    def pressure (self, navg=30, idx=(0,-1)):
        ''' p = read_pressure (filename, navg)
            lit les donnees calipso pressure'''
        p = self._read_var("Pressure", navg, idx=idx)
        return p

    def mol (self, navg=30, prof=None,idx=(0,-1)):
        ''' mol = read_mol (filename, navg)
            lit les profils moleculaires dans un fichier calipso'''
        mol = self._read_var('Molecular_Number_Density', navg, idx=idx)
        if prof:
            mol = mol[prof,:]
        return mol

    def temperature (self, navg=30,idx=(0,-1)):
        ''' t = read_temperature (filename, navg)
                lit les temperatures dans un fichier calipso'''
        t = self._read_var('Temperature', navg, idx=idx)
        return t

    def rh (self, navg=30, idx=(0,-1)):
        rh = self._read_var('Relative_Humidity', navg, idx=idx)
        return rh
        
    def pressure_on_lidar_alt (self, navg=30, alt=lidar_alt, metalt=met_alt, idx=(0,-1)):
        ''' lit la pression calipso en la moyennant sur navg profils,
            et en l'interpolant sur les niveaux d'altitude lidar.'''
        p0 = self.pressure (navg=navg,idx=idx)
        p = _remap_y (p0, metalt, alt)

        return p

    def mol_on_lidar_alt (self, navg=30,prof=None, alt=lidar_alt, metalt=met_alt,idx=(0,-1)):
        """Lit le moleculaire et l'interpole sur les altitudes lidar"""
        mol0 = self.mol (navg=navg, prof=prof,idx=idx)
        mol = _remap_y (mol0, metalt, alt)
        
        return mol

    def mol_calibration_coef(self, mol=None, atb=None, navg=30, prof=None, alt=lidar_alt, metalt=met_alt, idx=(0,-1), navgh=50, zmin=30, zmax=34):
        '''
        mol_calibration_coef(mol=mol, atb=atb, navg=navg, prof=prof, alt=lidar_alt, metalt=met_alt, idx=idx, navgh=100)
            Returns the molecular calibration coefficient, computed from atb 532 nm and molecular density profiles,
            averaged between zmin and zmax km vertically and over profiles [i-navgh:i+navgh] horizontally (moving average horizontally)
            
            Same strategy as CALIPSO NASA algorithm
        '''
        if mol is None and atb is None:            
            mol = self.mol_on_lidar_alt(navg=navg, prof=prof, alt=alt, metalt=metalt, idx=idx)
            atb = self.atb(navg=navg, prof=prof, idx=idx)

        # remove atb and molecular unfit for calibration purposes
        # this level of backscattering is most probably due to noise in the lower stratosphere
        # (and if it's not noise we don't want it anyway)
        atb[np.abs(atb) > atb_max[self.z]] = np.nan
        mol[mol < 0] = np.nan
            
        idx = (alt >= zmin) & (alt <= zmax)
        atb_calib_profile = nanmean(atb[:,idx], axis=1)
        mol_calib_profile = nanmean(mol[:,idx], axis=1)
         
        # now do a moving average, weeding out bad profiles
        atbbounds = iatb_bounds[self.z]

        nprof = mol.shape[0]        
        coef = np.zeros([nprof])
        
        for i in np.arange(nprof):
            
            idxh = np.r_[np.max([0, i-navgh]):np.min([nprof-1, i+navgh])]
            atbslice = atb_calib_profile[idxh]
            molslice = mol_calib_profile[idxh]
            idx = (atbslice > atbbounds[0]) & (atbslice < atbbounds[1])
            coef[i] = nanmean(atbslice[idx])/nanmean(molslice)
            
        return coef
        
    def calibration_data(self, mol=None, atb=None, navg=20, prof=None, alt=lidar_alt, metalt=met_alt, idx=(0,-1), navgh=5, zmin=30, zmax=34):
        '''
        for debugging purposes
        '''
        if atb is None:            
            mol = self.mol_on_lidar_alt(navg=navg, prof=prof, alt=alt, metalt=metalt, idx=idx)
            atb = self.atb(navg=navg, prof=prof, idx=idx)

        # remove atb unfit for calibration purposes
        # this level of backscattering is most probably due to noise in the lower stratosphere
        # (and if it's not noise we don't want it anyway)
        atb[np.abs(atb) > atb_max[self.z]] = np.nan

        nprof = atb.shape[0]        

        mol[mol < 0] = np.nan
        idx = (alt >= zmin) & (alt <= zmax)
        
        atb_for_normalization = atb[:,idx]
        mol_for_normalization = mol[:,idx]
        atb_for_normalization_by_profile = nanmean(atb[:,idx], axis=1)
        mol_for_normalization_by_profile = nanmean(mol[:,idx], axis=1)

        # now do a moving average, weeding out bad profiles
        atb_for_normalization_averaged = np.zeros(nprof)
        mol_for_normalization_averaged = np.zeros(nprof)
        
        bounds = iatb_bounds[self.z]

        for i in np.arange(nprof):
            
            idxh = np.r_[np.max([0, i-navgh]):np.min([nprof, i+navgh])]
            atbslice = atb_for_normalization_by_profile[idxh]
            molslice = mol_for_normalization_by_profile[idxh]
            idx = (atbslice > bounds[0]) & (atbslice < bounds[1])
            atb_for_normalization_averaged[i] = nanmean(atbslice[idx])
            mol_for_normalization_averaged[i] = nanmean(molslice)

        return (atb_for_normalization, mol_for_normalization, atb_for_normalization_averaged,
                atb_for_normalization_by_profile, mol_for_normalization_averaged, mol_for_normalization_by_profile)

    def mol_attenuate(self, mol, alt):
        '''
        Attenuate the molecular backscatter signal
        Make it a "attenuated molecular backscatter"
        Might be really slow...
        '''
        smol = 8.*np.pi/3.
        
        nprof, nalt = mol.shape

        for i in np.arange(nprof):
            molprof = mol[i,:]
            for ialt in np.arange(nalt):
                # compute integrated molecular backscatter down to ialt
                imol = _integrate_signal(molprof[0:ialt], alt[0:ialt])
                # integration extinction
                ialphamol = smol * imol
                # attenuate molecular
                mol[i,ialt] *= np.exp(-2.* ialphamol)
        
        return mol

    def mol_on_lidar_alt_calibrated(self, navg=30, prof=None, alt=lidar_alt, metalt=met_alt, idx=(0,-1), zcal=(30,34), attenuated=False):
        '''
        mol_on_lidar_alt_calibrated()
            returns profiles of visible molecular backscatter calibrated on the 532 nm channel
        '''
        mol = self.mol_on_lidar_alt(navg=navg, prof=prof, alt=alt, metalt=metalt, idx=idx)
        atb = self.atb(navg=navg, prof=prof, idx=idx)
        # taille de mol [nprof, nalt]
        
        coef = self.mol_calibration_coef(mol=mol, atb=atb, zmin=zcal[0], zmax=zcal[1])
        for i in np.arange(mol.shape[0]):
            mol[i,:] *= coef[i]
            
        if attenuated:
            mol = self.mol_attenuate(mol, alt)
        
        return mol
    
    def scattering_ratio(self, navg=30, prof=None, alt=lidar_alt, metalt=met_alt, idx=(0,-1)):
        '''
        scattering_ratio()
            returns the scattering ratio (atb532/mol) computed from atb532 and calibrated molecular density profiles
        '''
        mol_calib = self.mol_on_lidar_alt_calibrated(navg=navg, prof=prof, alt=lidar_alt, metalt=met_alt, idx=idx)
        atb = self.atb(navg=navg, prof=prof, idx=idx)
        r = atb/mol_calib
        return r
    
    def temperature_on_lidar_alt (self, navg=30, prof=None, alt=lidar_alt, metalt=met_alt, idx=(0,-1)):
        ''' t = read_temperature_remap (filename, navg)
                lit les temperatures dans un fichier calipso,
                les moyenne tous les navg profils,
                et les interpole sur les niveaux d'altitude lidar'''
        t0 = self.temperature (navg=navg, idx=idx)
        t = _remap_y (t0, metalt, alt)
        if prof:
            t = t[prof,:]
        return t

    def tropopause_height(self, navg=30, prof=None, idx=(0,-1)):
        tropoz = self._read_var('Tropopause_Height', navg, idx=idx)
        return tropoz
        
    def tropopause_temperature(self, navg=30, prof=None, idx=(0,-1)):
        tropot = self._read_var('Tropopause_Temperature', navg, idx=idx)
        return tropot
            

    def rh_on_lidar_alt (self, navg=30, prof=None, alt=lidar_alt, metalt=met_alt, idx=(0,-1)):
        rh0 = self.rh(navg=navg, idx=idx)
        rh = _remap_y(rh0, metalt, alt)
    
        if prof:
            rh = rh[prof,:]
        return rh

    def peek(self, lat, alt, data, latrange, dataname, vmin, vmax):
        import matplotlib.pyplot as plt

        print 'Showing ' + dataname
        print 'Lat range : ', latrange

        plt.ioff()
        plt.figure(figsize=(14,8))
        plt.pcolormesh(lat, alt, data.T, vmin=vmin, vmax=vmax)
        plt.xlim(latrange[0], latrange[1])
        plt.ylim(0, 22)
        plt.colorbar().set_label(dataname)
        plt.xlabel('Latitude')
        plt.ylabel('Altitude [km]')
        plt.title(dataname + ' - CALIOP '+ str(self.date) + ' - '+ self.orbit)
        
    def peek_atb(self, latrange=[-84,84], navg=120):
        
        lon,lat = self.coords(navg=navg)
        atb = self.atb(navg=navg)
        idx = (lat > latrange[0]) & (lat < latrange[1])
        lat = lat[idx]
        atb = atb[idx,:]
        print 'Averaging every %d profiles = %d shown profiles' % (navg, len(lat))
        self.peek(lat, lidar_alt, np.log10(atb), latrange, 'Attenuated Total Backscatter', -4, -2)
                
    def peek_depol(self, latrange=[-84,84], navg=120):

        lon,lat = self.coords(navg=navg)
        total = self.atb(navg=navg)
        perp = self.perp(navg=navg)
        para = total - perp
        depol = perp/para
        idx = (lat > latrange[0]) & (lat < latrange[1])
        lat = lat[idx]
        depol = depol[idx,:]
        print 'Averaging every %d profiles = %d shown profiles' % (navg, len(lat))

        self.peek(lat, lidar_alt, depol, latrange, 'Volumic Depolarization Ratio', 0, 0.8)

    def peek_colorratio(self, latrange=[-84,84], navg=120):

        lon,lat = self.coords(navg=navg)
        total = self.atb(navg=navg)
        total1064 = self.atb1064(navg=navg)
        cr = total1064 / total
        idx = (lat > latrange[0]) & (lat < latrange[1])
        lat = lat[idx]
        cr = cr[idx,:]
        print 'Averaging every %d profiles = %d shown profiles' % (navg, len(lat))

        self.peek(lat, lidar_alt, cr, latrange, 'Volumic Color Ratio', 0, 1.2)

# classe pour level 2

class Cal2(_Cal):

    def _read_var(self, var, idx=(0, -1)):
        var = self.hdf.select(var)
        if len(var.dimensions()) == 1:
            data = var[idx[0]:idx[1]]
        else:
            data = var[idx[0]:idx[1],:]
        var.endaccess()
        return data

    def lat(self, idx=(0, -1)):
        return self._read_var('Latitude', idx=idx)[:,1]

    def lon(self, idx=(0, -1)):
        return self._read_var('Longitude', idx=idx)[:,1]

    def coords(self, idx=(0,-1)):
        lat =self._read_var('Latitude', idx=idx)[:,1]
        lon =self._read_var('Longitude', idx=idx)[:,1]
        return lon, lat

    def off_nadir_angle (self, idx=(0, -1)):
        return self._read_var('Off_Nadir_Angle', idx=idx)

    def tropopause_height (self, idx=(0,-1)):
        return self._read_var('Tropopause_Height', idx=idx)[:,0]

    def tropopause_temperature(self, idx=(0,-1)):
        return self._read_var('Tropopause_Temperature', idx=idx)

    def dem_surface_elevation(self, idx=(0, -1)):
        return self._read_var('DEM_Surface_Elevation', idx=idx)

    ## Layer info

    def layers(self, idx=(0,-1)):
        nl = self._read_var('Number_Layers_Found', idx=idx)[:,0]
        top = self._read_var('Layer_Top_Altitude', idx=idx)
        base = self._read_var('Layer_Base_Altitude', idx=idx)
        return nl, base, top

    def layers_number(self, idx=(0,-1)):
        return self._read_var('Number_Layers_Found', idx=idx)[:,0]

    # in retrospect, this 'cloud_' prefix was a bad idea, I keep this for compatibility
    def cloud_layers (self, idx=(0,-1)):
        nl, base, top = self.layers(idx=idx)        

    def cloud_layers_number(self, idx=(0,-1)):
        return self.layers_number(idx=idx)
        
    def midlayer_temperature (self, idx=(0,-1)):
        return self._read_var('Midlayer_Temperature', idx=idx)

    def flag (self, idx=(0,-1)):
        return self._read_var('Feature_Classification_Flags', idx=idx)

    def layer_type (self, idx=(0,-1)):
        f = self.flag(idx=idx) 
        # type flag : bits 1 to 3
        typeflag = (f & 7)
        return typeflag
        
    def layer_subtype(self, idx=(0, -1)):
        f = self.flag(idx=idx)
        # subtype flag : bits 10 to 12
        subtype = (f & 3584) >> 9
        return subtype

    def layer_type_qa(self, idx=(0,-1)):
        f = self.flag(idx=idx)
        typeflag_qa = (f & 24) >> 3
        return typeflag_qa

    def phase (self, idx=(0,-1)):
        f = self.flag(idx=idx)
        # 96 = 0b1100000, bits 6 to 7
        cloudflag = (f & 96) >> 5
        return cloudflag

    def phase_qa (self, idx=(0,-1)):
        f = self.flag(hdf, idx=idx)
        cloudflag_qa = (f & 384) >> 7
        return cloudflag_qa

    def cloud_type (self, idx=(0,-1)):
        return self.phase(idx=idx)

    def opacity_flag(self, idx=(0,-1)):
        return self._read_var('Opacity_Flag', idx=idx)

    def horizontal_averaging(self, idx=(0, -1)):
        return self._read_var('Horizontal_Averaging', idx=idx)

    def iatb532 (self, idx=(0,-1)):
        return self._read_var('Integrated_Attenuated_Backscatter_532', idx=idx)

    def ivdp (self, idx=(0,-1)):
        return self._read_var('Integrated_Volume_Depolarization_Ratio', idx=idx)

    def ipdp (self, idx=(0,-1)):
        return self._read_var('Integrated_Particulate_Depolarization_Ratio', idx=idx)

    def icr (self, idx=(0,-1)):
        return self._read_var('Integrated_Attenuated_Total_Color_Ratio', idx=idx)
        
    def ipcr (self, idx=(0,-1)):
        return self._read_var('Integrated_Particulate_Color_Ratio', idx=idx)

    def od(self, idx=(0, -1)):
        return self._read_var('Feature_Optical_Depth_532', idx=idx)


import unittest
import random

def randomdate():
    year = random.randrange(2006, 2009)
    month = random.randrange(13)
    day = random.randrange(32)
    return year, month, day

class TestCal1(unittest.TestCase):    
    
    def setUp(self):
        # choisir un fichier CALIPSO au hasard...
        files = []
        while len(files)==0:
            files = l1_night_files(*randomdate())
        self.filename = files[random.randrange(len(files))]

        print 'Using file '+self.filename+', because why not'

        self.cal = Cal1(self.filename)
        
    def testFilename(self):
        self.assertEqual(self.filename, self.cal.__repr__())
        
    def testAverageTooMuchProfiles(self):
        lon, lat = self.cal.coords(100000000)
        
    def testAverageTooLittleProfiles(self):
        lon, lat = self.cal.coords(0)
        
    def testAveraging(self):
        lon, lat = self.cal.coords(1)
        nproftotal = lon.shape[0]
        lon, lat = self.cal.coords(15)
        nprof = lon.shape[0]
        self.assertEqual(nproftotal/15, nprof)
        
    def testNumberOfProfilesIsTheSameInVectors(self):
        navg = (1, 15, 19, 1000)
        lon, lat = self.cal.coords()
        elev = self.cal.surface_elevation()
        self.assertEqual(lon.shape[0], lat.shape[0])
        self.assertEqual(lon.shape[0], elev.shape[0])
        
    def testNumberOfProfilesIsTheSameInArrays(self):
        navg = (1, 15, 19, 1000)
        for n in navg:
            lon, lat = self.cal.coords(n)
            nprof = lon.shape[0]
            atb = self.cal.atb(n)
            perp = self.cal.perp(n)
            atb1064 = self.cal.atb1064(n)
            pressure = self.cal.pressure(n)
            mol = self.cal.mol(n)
            temperature = self.cal.temperature(n)
            rh = self.cal.rh(n)
            self.assertEqual(atb.shape[0], nprof)
            self.assertEqual(perp.shape[0], nprof)
            self.assertEqual(atb1064.shape[0], nprof)
            self.assertEqual(pressure.shape[0], nprof)
            self.assertEqual(mol.shape[0], nprof)
            self.assertEqual(temperature.shape[0], nprof)
            self.assertEqual(rh.shape[0], nprof)
        
if __name__ == '__main__':
    unittest.main()

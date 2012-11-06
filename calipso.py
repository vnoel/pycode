#!/usr/bin/env python
# encoding: utf-8

"""
Read CALIPSO lidar data files
level 1 and level 2.

V. Noel 2008-2011
LMD/CNRS
"""

# cannot use pytables since calipso files are hdf4.
from pyhdf.SD import SD, SDC
from scipy.integrate import trapz
import numpy as np
import numpy.ma as ma
import glob
import datetime
import os
import socket
import warnings


# these should be in the same directory as the module
datapath = os.path.dirname(__file__) + '/staticdata/'
lidar_alt = np.loadtxt(datapath + 'lidaralt.asc')
met_alt = np.loadtxt(datapath + 'metalt.asc')

# the following ranges have been calibrated more carefully for nighttime than
# daytime maximum atb to consider for calibration in the 26-28 km range
atb_max = {'ZN': 1e-4, 'ZD': 1}

# maximum integrated atb for calibrationi in the 26-28 km range
iatb_bounds = {'ZN': [1e-5, 8e-5], 'ZD': [-8e-3, 8e-3]}

# server-dependent paths
# these need to point to the path of CALIOP level 1 and level 2 data
# they might be better in a config file...
hostname = socket.gethostname()
if hostname == 'access.icare.univ-lille1.fr':
    l1dir = ('/DATA/LIENS/CALIOP/CAL_LID_L1.v3.01', None)
    l2dir = ('/DATA/LIENS/CALIOP/05kmCLay.v3.02',
        '/DATA/LIENS/CALIOP/05kmCLay.v3.01')
    l2adir = ('/DATA/LIENS/CALIOP/05kmALay.v3.01',)
else:
    l1dir = ('/bdd/CALIPSO/Lidar_L1/CAL_LID_L1.v3.02/', 
            '/bdd/CALIPSO/Lidar_L1/CAL_LID_L1.v3.01/')
    l15dir = ('/bdd/CFMIP/CAL_LID_L1.5',)
    # for level 2, look first in our own data, then on climserv data
    l2dir = ('/bdd/CFMIP/Lidar_L2',
        '/bdd/CALIPSO/Lidar_L2/05kmCLay.v3.01',
        '/bdd/CALIPSO/Lidar_L2/05kmCLay.v2.02',
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

def l2_night_files (y, m, d):
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

def l2bis_night_files (y, m, d):
    return l2bis_files(y, m, d, '*ZN*')

def l2bis_file_from_orbit(y, m, d, orbit):
    files = l2bis_files(y, m, d, '*' + orbit + '*')
    if not files:
        return None
    return files[0]

def read_l2bis (calfile):
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


# Fonctions maths utiles

def _integrate_signal(data, alt):
    integrale = trapz(data[::-1], x=alt[::-1])
    return integrale


def _vector_average(v0, navg):
    """ v = _vector_average (v0, navg)
        moyenne le vector v0 tous les navg points."""

    v0 = v0.squeeze()

    assert v0.ndim == 1, 'in _vector_average, v0 should be a vector'
    if navg == 1:
        return v0
        
    n = np.size(v0, 0) / navg
    v = np.zeros(n)
    for i in np.arange(n):
        n0 = i * navg
        v[i] = v0[n0:n0 + navg - 1].mean()
    return v


def _array_average(a0, navg, weighted=False):
    """ a = _array_average (a0, navg, weighted=False)
        moyenne le tableau a0 le long des x tous les navg profils."""

    a0 = a0.squeeze()

    assert a0.ndim == 2, 'in _array_average, a0 should be a 2d array'
    if navg == 1:
        return a0

    if weighted and (navg % 2 != 1):
        weighted = False
        print "_array_average: navg is even, turning weights off"

    # create triangle-shaped weights
    if weighted:
        w = np.zeros(navg)
        w[:navg / 2. + 1] = np.r_[1:navg / 2. + 1]
        w[navg / 2. + 1:] = np.r_[int(navg / 2.):0:-1]

    # create averaged array a with number of averaged profiles n
    n = np.size(a0, 0) / navg
    a = np.zeros([n, np.size(a0, 1)])
    
    for i in np.arange(n):
        n0 = i * navg
        if weighted:
            a[i, :] = np.average(a0[n0:n0 + navg, :], axis=0, weights=w)
        else:
            a[i, :] = np.mean(a0[n0:n0 + navg - 1, :],axis=0)
    return a


def _remap_y(z0, y0, y):
    """ z = remap (z0, y0, y)
            interpole les donnees du tableau z0 sur un nouveau y."""

    z = np.zeros([np.size(z0, 0), np.size(y)])
    for i in np.arange(np.size(z, 0)):
        z[i, :] = np.interp(y, np.flipud(y0), np.flipud(z0[i, :]))
    return z


def ncep_remap(var, lat, lon, latorb, lonorb):

    n = np.size(latorb, 0)
    var2 = np.zeros([n, np.size(var, 0)])

    for i in np.arange(n):
        d = np.sqrt((lat - latorb[i]) ** 2)
        xlat = d.argmin()
        d = np.sqrt((lon - lonorb[i]) ** 2)
        xlon = d.argmin()
        var2[i, :] = var[:, xlat, xlon]

    return var2


# Generic calipso file class
# Use Cal1 and Cal2 classes instead


class _Cal:
    
    """
    Trying to open a non-existing CALIOP file gives an exception
    """

    def __init__(self, filename):
        warnings.simplefilter('ignore', DeprecationWarning)
        self.hdf = SD(filename, SDC.READ)
        self.filename = filename
        self.orbit = filename[-15:-4]
        self.z = self.orbit[-2:]            # zn or zd
        self.year = int(filename[-25:-21])
        self.month = int(filename[-20:-18])
        self.day = int(filename[-17:-15])
        self.hour = int(filename[-14:-12])
        self.minutes = int(filename[-11:-9])
        self.seconds = int(filename[-8:-6])
        self.date = datetime.datetime(self.year, self.month, self.day,
            self.hour, self.minutes, self.seconds)
        
    def __repr__(self):
        return self.filename
        
    def close(self):
        self.hdf.end()
        self.hdf = None


class Cal1(_Cal):
    
    """
    Class to process CALIOP Level 1 file.
    Most of the reading function accept an argument navg 
    which requests horizontal
    averaging along the provided number of profiles.
    example use:
        
        from calipso import Cal1
        
        c = Cal1('CAL_LID_L1-ValStage1-V3-01.2010-02-05T01-47-40ZN.hdf')
        time = c.time(navg=15)
        lon, lat = c.coords(navg=15)
        atb = c.atb(navg=15)
        ...
        c.close()
        
    """

    def _read_var(self, var, navg, idx=(0, -1)):
        """
        Read a variable in an hdf file, averaging the data if navg is > 1
        """
        try:
            var = self.hdf.select(var)
        except:
            print var+' not present in file '+self.filename+' (version too old, maybe ?)'
            return None
            
        this_var = var[:].squeeze()
        if this_var.ndim == 1:
            
            if idx[0] is 0 and idx[1] is -1:
                data = this_var[:]
            else:
                data = this_var[idx[0]:idx[1]]
                
            if navg > 1:
                data = _vector_average(data, navg)
        else:
            
            if idx[0] is 0 and idx[1] is -1:
                data = this_var[...]
            else:
                data = this_var[idx[0]:idx[1], :]
            if navg > 1:
                data = _array_average(data, navg)

        var.endaccess()

        return data

    def utc_time(self, navg=30, idx=(0, -1)):
        var = self.hdf.select('Profile_UTC_Time')
        time = var[:].squeeze()
        time = time[idx[0]:idx[1]]
        if navg > 1:
            n = np.size(time, 0) / navg
            time2 = np.zeros(n)
            for i in np.arange(n):
                time2[i] = time[i * navg]
            time = time2

        var.endaccess()
        return time

    def datetimes(self, navg=30):

        utc = self.utc_time(navg=navg)
        datetimes = []
        for u in utc:
            y = np.floor(u / 10000.)
            m = np.floor((u - y * 10000) / 100.)
            d = np.floor(u - y * 10000 - m * 100)
            y = int(y) + 2000
            m = int(m)
            d = int(d)
            seconds_into_day = np.int((u - np.floor(u)) * 24. * 3600.)
            profile_datetime = datetime.datetime(y, m, d, 0, 0, 0) + \
                datetime.timedelta(seconds=seconds_into_day)
            datetimes.append(profile_datetime)
        return np.array(datetimes)

    def time(self, navg=30, idx=(0, -1)):
        """
        Reads time for CALIOP profiles, averaged on navg
        shape [nprof]
        Example:
            time = c.time(navg=15)
        """
        var = self.hdf.select('Profile_Time')
        time = var[:].squeeze()
        time = time[idx[0]:idx[1]]

        if navg > 1:
            n = np.size(time, 0) / navg
            time2 = np.zeros(n)
            for i in np.arange(n):
                time2[i] = time[i * navg]
            time = time2

        var.endaccess()

        return time

    def coords(self, navg=30, idx=(0, -1)):
        """
        Reads coordinates for CALIOP profiles, averaged on navg
        shape [nprof]
        Example:
            lon, lat = c.coords(navg=15)
        """
        lat = self._read_var('Latitude', navg, idx=idx)
        lon = self._read_var('Longitude', navg, idx=idx)

        return lon, lat

    def surface_elevation(self, navg=30, prof=None, idx=(0, -1)):
        """
        Reads the surface elevation from CALIOP file
        shape [nprof]
        """
        elev = self._read_var('Surface_Elevation', navg, idx=idx)
        if prof is not None:
            elev = elev[prof, :]
        return elev
        
    def perp(self, navg=30, prof=None, idx=(0, -1)):
        """
        Reads the perpendicular signal from CALIOP file
        shape [nprof, nz]
        """
        perp = self._read_var('Perpendicular_Attenuated_Backscatter_532',
            navg, idx=idx)
        if prof:
            perp = perp[prof, :]
        return perp

    def atb(self, navg=30, prof=None, idx=(0, -1)):
        """
        Reads the Attenuated Total Backscatter 532nm from CALIOP file
        shape [nprof, nz]
        """
        atb = self._read_var('Total_Attenuated_Backscatter_532', navg, idx=idx)
        if prof:
            atb = atb[prof, :]
        return atb
        
    def ccu_532(self, navg=30, prof=None, idx=(0, -1)):
        au = self._read_var('Calibration_Constant_Uncertainty_532', navg, idx=idx)
        if prof:
            au = au[prof,:]
        return au

    def atb1064(self, navg=30, prof=None, idx=(0, -1)):
        """
        Reads the Attenuated Total Backscatter 1064nm from CALIOP file
        """
        atb = self._read_var('Attenuated_Backscatter_1064', navg, idx=idx)
        if prof:
            atb = atb[prof, :]
        return atb

    def pressure(self, navg=30, idx=(0, -1)):
        """
        Reads the ancillary pressure field from CALIOP file
        shape [nprof, nlevels]
        """
        p = self._read_var("Pressure", navg, idx=idx)
        return p

    def mol(self, navg=30, prof=None, idx=(0, -1)):
        """
        Reads the ancillary molecular number density from CALIOP file
        shape [nprof, nlevels]
        """
        mol = self._read_var('Molecular_Number_Density', navg, idx=idx)
        if prof:
            mol = mol[prof, :]
        return mol

    def temperature(self, navg=30, idx=(0, -1)):
        """
        Reads the ancillary temperature field in degC from CALIOP file
        shape [nprof, nlevels]
        """
        t = self._read_var('Temperature', navg, idx=idx)
        return t

    def rh(self, navg=30, idx=(0, -1)):
        """
        Reads the ancillary relative humidity field from CALIOP file
        shape [nprof, nlevels]
        """
        rh = self._read_var('Relative_Humidity', navg, idx=idx)
        return rh
        
    def pressure_on_lidar_alt(self, navg=30, alt=lidar_alt,
        metalt=met_alt, idx=(0, -1)):
        """
        Reads the ancillary pressure field from CALIOP file, interpolated on
        CALIOP altitude levels.
        shape [nprof, nz]
        """
        p0 = self.pressure(navg=navg, idx=idx)
        p = _remap_y(p0, metalt, alt)

        return p

    def mol_on_lidar_alt(self, navg=30, prof=None, alt=lidar_alt,
        metalt=met_alt, idx=(0, -1)):
        """
        Reads the ancillary molecular number density from CALIOP file,
        interpolated on
        CALIOP altitude levels.
        shape [nprof, nz]
        """
        mol0 = self.mol(navg=navg, prof=prof, idx=idx)
        mol = _remap_y(mol0, metalt, alt)

        return mol

    def mol_calibration_coef(self, mol=None, atb=None, navg=30, prof=None,
            alt=lidar_alt, metalt=met_alt, idx=(0, -1), navgh=50, zmin=30,
            zmax=34):
        """
        Returns the molecular calibration coefficient, computed from
        atb 532 nm and molecular density profiles,
        averaged between zmin and zmax km vertically, and [i-navgh:i+navgh]
        profiles horizontally using a moving average.
        shape [nprof]
        """
        if mol is None and atb is None:
            mol = self.mol_on_lidar_alt(navg=navg, prof=prof, alt=alt,
                metalt=metalt, idx=idx)
            atb = self.atb(navg=navg, prof=prof, idx=idx)

        # remove atb and molecular unfit for calibration purposes
        # this level of backscattering is most probably due to noise
        # in the lower stratosphere
        # (and if it's not noise we don't want it anyway)
        
        atb[np.abs(atb) > atb_max[self.z]] = np.nan
        mol[mol < 0] = np.nan
        
        atb = ma.masked_invalid(atb)
        mol = ma.masked_invalid(mol)

        idx = (alt >= zmin) & (alt <= zmax)
        
        atb_calib_profile = ma.mean(atb[:, idx], axis=1)
        mol_calib_profile = ma.mean(mol[:, idx], axis=1)
         
        # now do a moving average, weeding out bad profiles
        atbbounds = iatb_bounds[self.z]

        nprof = mol.shape[0]
        coef = np.zeros([nprof])
        
        for i in np.arange(nprof):
            
            idxh = np.r_[np.max([0, i - navgh]):np.min([nprof - 1, i + navgh])]
            atbslice = atb_calib_profile[idxh]
            molslice = mol_calib_profile[idxh]
            idx = (atbslice > atbbounds[0]) & (atbslice < atbbounds[1])
            coef[i] = ma.mean(atbslice[idx]) / ma.mean(molslice)
            
        return coef
        
    def mol_on_lidar_alt_calibrated(self, navg=30, prof=None, alt=lidar_alt,
         navgh=50, metalt=met_alt, idx=(0, -1), zcal=(30, 34), atb=None):
        """
        Returns an estimate of the molecular backscatter at 532 nm, computed
        from the molecular number density calibrated on
        the attenuated total backscatter at 532 nm (both from the CALIOP file),
        using the calibration coefficient returned
        from mol_calibration_coef().
        Shape: [nprof, nz]
        
        atb can be passed as an argument if it's been read and averaged already.
        """
        mol = self.mol_on_lidar_alt(navg=navg, prof=prof, alt=alt,
            metalt=metalt, idx=idx)
        if atb is None:
            atb = self.atb(navg=navg, prof=prof, idx=idx)
        
        coef = self.mol_calibration_coef(mol=mol, atb=atb, zmin=zcal[0],
            zmax=zcal[1], navgh=navgh)
            
        # x * y is equivalent to x[:,i] * y[i]
        # mol is [i,:], so we need to rotate it twice
        mol = mol.T
        mol *= coef
        mol = mol.T
            
        return mol
    
    def temperature_on_lidar_alt(self, navg=30, prof=None, alt=lidar_alt,
        metalt=met_alt, idx=(0, -1)):
        """
        Returns the ancillary temperature field in degC, interpolated on
        CALIOP altitude levels.
        """
        t0 = self.temperature(navg=navg, idx=idx)
        t = _remap_y(t0, metalt, alt)
        if prof:
            t = t[prof, :]
        return t

    def rh_on_lidar_alt(self, navg=30, prof=None, alt=lidar_alt,
        metalt=met_alt, idx=(0, -1)):
        """
        Returns the ancillary relative humidity field from the CALIOP file,
        interpolated on CALIOP altitude levels.
        """
        rh0 = self.rh(navg=navg, idx=idx)
        rh = _remap_y(rh0, metalt, alt)

        if prof:
            rh = rh[prof, :]
        return rh

    def scattering_ratio(self, navg=30, prof=None, alt=lidar_alt,
        metalt=met_alt, idx=(0, -1)):
        """
        Returns the scattering ratio, i.e. the ratio between the attenuated
        total backscatter and the molecular
        backscatter, both at 532 nm, both read from the CALIOP file.
        """
        mol_calib = self.mol_on_lidar_alt_calibrated(navg=navg,
            prof=prof, alt=lidar_alt, metalt=met_alt, idx=idx)
        atb = self.atb(navg=navg, prof=prof, idx=idx)
        r = atb / mol_calib
        return r

    def tropopause_height(self, navg=30, prof=None, idx=(0, -1)):
        """
        Reads the ancillary tropopause height, in km, from the CALIOP file.
        shape [nprof]
        """
        tropoz = self._read_var('Tropopause_Height', navg, idx=idx)
        return tropoz
        
    def tropopause_temperature(self, navg=30, idx=(0, -1)):
        """
        Reads the ancillary tropopause temperature, in degC,
        from the CALIOP file.
        shape [nprof]
        """
        tropot = self._read_var('Tropopause_Temperature', navg, idx=idx)
        return tropot

    # some display utility functions
    
    def _peek(self, lat, alt, data, latrange, dataname, vmin, vmax,
        datetime=False, ymin=0, ymax=25, cmap=None):

        import niceplots
        
        import matplotlib.pyplot as plt
        from pylab import get_cmap

        print 'Showing ' + dataname
        print 'Lat range : ', latrange
        
        if cmap is None:
            cmap = get_cmap('gist_stern_r')
        
        plt.ioff()
        fig = plt.figure(figsize=(20, 6))
        ax = plt.gca()
        plt.pcolormesh(lat, alt, data.T, vmin=vmin, vmax=vmax, cmap=cmap)
        plt.xlim(latrange[0], latrange[1])
        plt.ylim(ymin, ymax)
        if datetime:
            niceplots.axis_set_date_format(ax, format='%H:%M')
            fig.autofmt_xdate()
        else:
            plt.xlabel('Latitude')
            
        plt.colorbar().set_label(dataname)
        plt.ylabel('Altitude [km]')
        plt.title(dataname + ' - CALIOP ' + str(self.date))
            
    def peek_atb_time(self, navg=30, mintime=None, maxtime=None):

        from pylab import get_cmap

        lon, lat = self.coords(navg=navg)
        atb = self.atb(navg=navg)
        time = self.datetimes(navg=navg)
        import matplotlib.dates as mdates
        numtime = mdates.date2num(time)
        
        if mintime is None:
            mintime = np.min(numtime)
        else:
            mintime = mdates.date2num(mintime)
        
        if maxtime is None:
            maxtime = np.max(numtime)
        else:
            maxtime = mdates.date2num(maxtime)
        
        self._peek(numtime, lidar_alt, np.log10(atb), [mintime, maxtime],
            'Coefficient de retrodiffusion [log10]', -3.5, -2, datetime=True,
            cmap=get_cmap('gist_stern_r'))

    def peek_cr_time(self, navg=30, mintime=None, maxtime=None):

        from pylab import get_cmap

        lon, lat = self.coords(navg=navg)
        atb532 = self.atb(navg=navg)
        atb1064 = self.atb1064(navg=navg)
        cr = atb1064 / atb532
        time = self.datetimes(navg=navg)
        import matplotlib.dates as mdates
        numtime = mdates.date2num(time)

        if mintime is None:
            mintime = np.min(numtime)
        else:
            mintime = mdates.date2num(mintime)
        if maxtime is None:
            maxtime = np.max(numtime)
        else:
            maxtime = mdates.date2num(maxtime)

        self._peek(numtime, lidar_alt, cr, [mintime, maxtime],
            'Volumic Attenuated Color Ratio', 0, 1.4, datetime=True,
            cmap=get_cmap('jet'))
        
    def peek_psc_time(self, navg=30, mintime=None, maxtime=None):
        atb = self.atb(navg=navg)
        time = self.datetimes(navg=navg)
        import matplotlib.dates as mdates
        
        print np.min(time), np.max(time)
        print mintime, maxtime

        numtime = mdates.date2num(time)

        if mintime is None:
            mintime = np.min(numtime)
        else:
            mintime = mdates.date2num(mintime)
            if mintime < np.min(numtime):
                print 'Error : mintime < min(time)'
            
        if maxtime is None:
            maxtime = np.max(numtime)
            if maxtime > np.max(numtime):
                print 'Error : maxtime > max(time)'
        else:
            maxtime = mdates.date2num(maxtime)

        tropo = self.tropopause_height(navg=navg)
        tropo[tropo > 13] = 11
        for i, z in enumerate(tropo):
            idx = (lidar_alt < (z + 3))
            atb[i, idx] = atb[i, idx] * 0.5
        atb[atb < 0] = 1e-6

        self._peek(numtime, lidar_alt, np.log10(atb), [mintime, maxtime],
        'Coefficient de retrodiffusion (log10)', -3.5, -2.25, datetime=True,
        ymin=5, ymax=28)
        import matplotlib.pyplot as plt
        plt.plot(numtime, tropo + 1, color='green')

    def peek_atb(self, latrange=(-84, 84), navg=120):
        """
        Display a quick look at the attenuated total backscatter as a function of latitude
        and altitude.
        """

        lon, lat = self.coords(navg=navg)
        atb = self.atb(navg=navg)
        idx = (lat > latrange[0]) & (lat < latrange[1])
        lat = lat[idx]
        atb = atb[idx, :]
        print 'Averaging every %d profiles = %d shown profiles' % (navg, len(lat))
        self._peek(lat, lidar_alt, np.log10(atb), latrange,
            'Attenuated Total Backscatter', -4, -2)

    def peek_depol(self, latrange=(-84, 84), navg=120):
        """
        Display a quick look at the depolarization ratio as a function of latitude
        and altitude.
        """

        lon, lat = self.coords(navg=navg)
        total = self.atb(navg=navg)
        perp = self.perp(navg=navg)
        para = total - perp
        depol = perp / para
        idx = (lat > latrange[0]) & (lat < latrange[1])
        lat = lat[idx]
        depol = depol[idx, :]
        print 'Averaging every %d profiles = %d shown profiles' % (navg, len(lat))

        self._peek(lat, lidar_alt, depol, latrange,
            'Volumic Depolarization Ratio', 0, 0.8)

    def peek_colorratio(self, latrange=(-84, 84), navg=120):
        """
        Display a quick look at the color ratio as a function of latitude and altitude.
        """

        lon, lat = self.coords(navg=navg)
        total = self.atb(navg=navg)
        total1064 = self.atb1064(navg=navg)
        cr = total1064 / total
        idx = (lat > latrange[0]) & (lat < latrange[1])
        lat = lat[idx]
        cr = cr[idx, :]
        print 'Averaging every %d profiles = %d shown profiles' % (navg, len(lat))

        self._peek(lat, lidar_alt, cr, latrange,
            'Volumic Color Ratio', 0.2, 1.)  # , cmap=get_cmap('jet'))


class Cal2(_Cal):

    """
    Class to process CALIOP Level 2 files.
    No averaging is possible here given the nature of some variables.
    example use:
        
        from calipso import Cal2
        
        c = Cal2('CAL_LID_L2_05kmCLay-Prov-V3-01.2010-12-31T01-37-30ZN.hdf')
        lon, lat = c.coords()
        nl, base, top = c.layers()
        ...
        c.close()
        
    """

    def _read_var(self, var, idx=(0, -1)):
        hdfvar = self.hdf.select(var)
        if idx[0] is 0 and idx[1] is -1:
            data = hdfvar[:]
        else:
            if len(hdfvar.dimensions()) == 1:
                data = hdfvar[idx[0]:idx[1]]
            else:
                data = hdfvar[idx[0]:idx[1], :]
        hdfvar.endaccess()
        return data

    def lat(self, idx=(0, -1)):
        """
        Returns latitude for profiles.
        shape [nprof]
        """
        return self._read_var('Latitude', idx=idx)[:, 1]

    def lon(self, idx=(0, -1)):
        """
        Returns longitude for profiles.
        shape [nprof]
        """
        return self._read_var('Longitude', idx=idx)[:, 1]

    def coords(self, idx=(0, -1)):
        """
        Returns longitude and latitude for profiles.
        shape [nprof]
        """
        lat = self._read_var('Latitude', idx=idx)[:, 1]
        lon = self._read_var('Longitude', idx=idx)[:, 1]
        return lon, lat

    def time(self, idx=(0, -1)):
        '''
        returns profile time (TAI)
        shape [nprof]
        '''
        return self._read_var('Profile_Time', idx=idx)[:]

    def utc_time(self, idx=(0, -1)):
        '''
        Returns utc time value (decimal time)
        '''
        time = self._read_var('Profile_UTC_Time', idx=idx)[:]
        return time

    def datetime(self, idx=(0, -1)):
        '''
        Returns an array of datetime objects based on utc_time values
        '''
        utc = self.utc_time(idx=idx)
        datetimes = []
        for u in utc:
            y = np.floor(u / 10000.)
            m = np.floor((u - y * 10000) / 100.)
            d = np.floor(u - y * 10000 - m * 100)
            y = int(y) + 2000
            m = int(m)
            d = int(d)
            seconds_into_day = np.int((u - np.floor(u)) * 24. * 3600.)
            profile_datetime = datetime.datetime(y, m, d, 0, 0, 0) + \
                datetime.timedelta(seconds=seconds_into_day)
            datetimes.append(profile_datetime)
        return np.array(datetimes)

    def off_nadir_angle(self, idx=(0, -1)):
        """
        Returns the off-nadir-angle, in deg, for profiles.
        shape [nprof]
        """
        return self._read_var('Off_Nadir_Angle', idx=idx)

    def tropopause_height(self, idx=(0, -1)):
        """
        Returns the ancillary tropopause height, in km, for profiles.
        shape [nprof]
        """
        return self._read_var('Tropopause_Height', idx=idx)[:, 0]

    def tropopause_temperature(self, idx=(0, -1)):
        """
        Returns the ancillary tropopause temperature, in degC, for profiles
        shape [nprof]
        """
        return self._read_var('Tropopause_Temperature', idx=idx)

    def dem_surface_elevation(self, idx=(0, -1)):
        """
        Returns the ancillary surface elevation (from digital elevation model)
        in km, for profiles
        shape [nprof]
        """
        return self._read_var('DEM_Surface_Elevation', idx=idx)

    def lidar_surface_elevation(self, idx=(0, -1)):
        """
        returns 8 values per profile
        min, max, mean, stddev for upper boundary of the surface echo.
        min, max, mean, stddev for lower boundary of the surface echo.
        en kilometres 
        """
        
        return self._read_var('Lidar_Surface_Elevation', idx=idx)

    # Layer info

    def layers(self, idx=(0, -1)):
        """
        Returns layer information by profile:
        nl = number of layers, shape [nprof]
        top = layer top, shape [nprof, nlaymax]
        base = layer base, shape [nprof, nlaymax]
        """

        nl = self._read_var('Number_Layers_Found', idx=idx)[:, 0]
        top = self._read_var('Layer_Top_Altitude', idx=idx)
        base = self._read_var('Layer_Base_Altitude', idx=idx)
        return nl, base, top

    def layers_number(self, idx=(0, -1)):
        """
        Returns the number of layer found by profile
        shape [nprof]
        """
        return self._read_var('Number_Layers_Found', idx=idx)[:, 0]

    def midlayer_temperature(self, idx=(0, -1)):
        """
        Returns the midlayer temperature by layer, in km
        shape [nprof, nlaymax]
        """
        return self._read_var('Midlayer_Temperature', idx=idx)

    def flag(self, idx=(0, -1)):
        """
        Returns the feature classification flag by layer
        shape [nprof, nlaymax]
        """
        return self._read_var('Feature_Classification_Flags', idx=idx)

    def layer_type(self, idx=(0, -1)):
        """
        Returns the layer type from the feature classification flag
        shape [nprof, nlaymax]
        """
        f = self.flag(idx=idx)
        # type flag : bits 1 to 3
        typeflag = (f & 7)
        return typeflag

    def layer_subtype(self, idx=(0, -1)):
        """
        Returs the layer subtype, as identified from the feature
        classification flag
        shape [nprof, nlaymax]
        """
        f = self.flag(idx=idx)
        # subtype flag : bits 10 to 12
        subtype = (f & 3584) >> 9
        return subtype

    def layer_type_qa(self, idx=(0, -1)):
        """
        Returs the quality flag for the layer type, as identified from the
        feature classification flag
        shape [nprof, nlaymax]
        """
        f = self.flag(idx=idx)
        typeflag_qa = (f & 24) >> 3
        return typeflag_qa

    def phase(self, idx=(0, -1)):
        """
        Returs the layer thermodynamical phase, as identified from the
        feature classification flag
        shape [nprof, nlaymax]
        """
        f = self.flag(idx=idx)
        # 96 = 0b1100000, bits 6 to 7
        cloudflag = (f & 96) >> 5
        return cloudflag

    def phase_qa(self, idx=(0, -1)):
        """
        Returns the quality flag for the layer thermodynamical phase,
        as identified from the feature classification flag
        shape [nprof, nlaymax]
        """

        # f = self.flag(hdf, idx=idx)
        # cloudflag_qa = (f & 384) >> 7
        # return cloudflag_qa

    def opacity_flag(self, idx=(0, -1)):
        """
        Returns the opacity flag by layer.
        shape [nprof, nlaymax]
        """

        return self._read_var('Opacity_Flag', idx=idx)

    def horizontal_averaging(self, idx=(0, -1)):
        """
        Returns the horizontal averaging needed to detect a given layer.
        shape [nprof, nlaymax]
        """
        return self._read_var('Horizontal_Averaging', idx=idx)

    def iatb532(self, idx=(0, -1)):
        """
        Returns the integrated attenuated total backscatter at 532 nm
        along the layer thickness.
        shape [nprof, nlaymax]
        """
        return self._read_var('Integrated_Attenuated_Backscatter_532', idx=idx)

    def ivdp(self, idx=(0, -1)):
        """
        Returns the volumic depolarization ratio for the entire layer
        thickness, obtained by the
        ratio of integrated perpendicular backscatter on the integrated
        parallel backscatter at 532 nm
        shape [nprof, nlaymax]
        """
        return self._read_var(
            'Integrated_Volume_Depolarization_Ratio', idx=idx)

    def ipdp(self, idx=(0, -1)):
        """
        Returns the particulate depolarization ratio for the entire
        layer thickness, i.e. the volumic
        depolarization ratio corrected to remove its molecular component.
        shape [nprof, nlaymax]
        """
        return self._read_var(
            'Integrated_Particulate_Depolarization_Ratio', idx=idx)

    def icr(self, idx=(0, -1)):
        """
        Returns the integrated attenuated color ratio for the entire
        layer thickness.
        shape [nprof, nlaymax]
        """
        return self._read_var(
            'Integrated_Attenuated_Total_Color_Ratio', idx=idx)

    def ipcr(self, idx=(0, -1)):
        """
        Returns the integrated color ratio for the entire layer thickness
        shape [nprof, nlaymax]
        """
        return self._read_var('Integrated_Particulate_Color_Ratio', idx=idx)

    def od(self, idx=(0, -1)):
        """
        Returns the optical depth found by layer.
        shape [nprof, nlaymax]
        """
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
        while len(files) == 0:
            files = l1_night_files(*randomdate())
        self.filename = files[random.randrange(len(files))]

        print 'Using file ' + self.filename + ', because why not'

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
        self.assertEqual(nproftotal / 15, nprof)

    def testNumberOfProfilesIsTheSameInVectors(self):
        # navg = (1, 15, 19, 1000)
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

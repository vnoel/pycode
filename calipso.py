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
import datetime
import os
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


# Useful maths

def _integrate_signal(data, alt):
    integrale = trapz(data[::-1], x=alt[::-1])
    return integrale


def _vector_average(v0, navg, missing=None, valid=None):
    """ v = _vector_average (v0, navg)
        moyenne le vector v0 tous les navg points."""

    v0 = v0.squeeze()

    assert v0.ndim == 1, 'in _vector_average, v0 should be a vector'
    if navg == 1:
        return v0

    n = np.size(v0, 0) / navg
    v = np.zeros(n)
    if valid is None:
        valid = np.ones_like(v0)

    for i in np.arange(n):
        n0 = i * navg
        vslice = v0[n0:n0 + navg - 1]
        validslice = valid[n0:n0 + navg - 1]
        if missing is None:
            idx = (validslice != 0) & (vslice > -999.)
            if idx.sum() == 0:
                v[i] = -9999.
            else:
                v[i] = vslice[idx].mean()
        else:
            idx = (vslice != missing) & (validslice != 0) & (vslice > -999.)
            if idx.sum() == 0:
                v[i] = None
            else:
                v[i] = vslice[idx].mean()
    return v


def _array_std(a0, navg, valid=None):
    a0 = a0.squeeze()
    assert a0.ndim == 2, 'in _array_std, a0 should be a 2d array'
    if navg == 1:
        return np.zeros_like(a0)
    n = np.size(a0, 0) / navg
    a = np.zeros([n, np.size(a0, 1)])

    if valid is None:
        valid = np.ones(np.size(a0, 0))

    for i in np.arange(n):
        n0 = i * navg
        aslice = a0[n0:n0 + navg - 1, :]
        validslice = valid[n0:n0 + navg - 1]
        idx = (validslice > 0)
        if idx.sum() == 0:
            a[i, :] = -9999.
        else:
            a[i, :] = np.std(aslice[idx, :], axis=0)
    return a


def _array_average(a0, navg, weighted=False, valid=None, missing=None):
    """
    a = _array_average (a0, navg, weighted=False)
    moyenne le tableau a0 le long des x tous les navg profils.
    missing = valeur a ignorer (genre -9999), ou None  
    precising a missing value might slow things down...      
    """

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
    else:
        w = None

    # create averaged array a with number of averaged profiles n
    n = np.size(a0, 0) / navg
    a = np.zeros([n, np.size(a0, 1)])
    if valid is None:
        # si on n'a pas d'info sur les profils valides, on suppose qu'ils le sont tous
        valid = np.ones(np.size(a0, 0))

    for i in np.arange(n):
        n0 = i * navg
        aslice = a0[n0:n0 + navg - 1, :]
        validslice = valid[n0:n0 + navg - 1]

        if missing is None:

            idx = (validslice != 0)  # & (np.all(aslice > -999., axis=1))

            if idx.sum() == 0:
                # no valid profiles in the slice according to valid profiles data
                a[i, :] = -9999.
            else:
                aslice = aslice[idx, :]
                # there might be invalid points somewhere in those profiles
                # find number of valid points along the vertical
                npts = np.sum(aslice > -999., axis=0)
                # sum only the valid points along the vertical
                aslice[aslice < -999.] = 0
                asliceavg = np.sum(aslice, axis=0) / npts
                # mark averaged points with no valid raw points as invalid
                asliceavg[npts == 0] = -9999.
                a[i, :] = asliceavg
        else:
            aslice = np.ma.masked_where(missing == aslice, aslice)
            a[i, :] = np.ma.mean(aslice, axis=0, weights=w)
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
        self.z = self.orbit[-2:]  # zn or zd
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

    def __init__(self, filename, max_rms=None):

        _Cal.__init__(self, filename)
        self.valid_rms_profiles = None
        if max_rms is not None:
            rms = self.parallel_rms_baseline(navg=1)
            self.valid_rms_profiles = (rms < max_rms)

    def _read_std(self, var, navg, idx=(0, -1)):
        """
        Reads a variable in an hdf file, and computes the standard deviation
        of the variable over navg profiles
        """

        try:
            var = self.hdf.select(var)
        except:
            print 'Cannot read ' + var + ' in ' + self.filename
            return None

        this_var = var[:].squeeze()
        if this_var.ndim == 1:
            print 'sorry, ndim=1 not implemented in _read_std'
            return None
        data = this_var[...]
        data = _array_std(data, navg, valid=self.valid_rms_profiles)

        var.endaccess()

        return data

    def _read_var(self, var, navg, idx=(0, -1), missing=None):
        """
        Read a variable in an hdf file, averaging the data if navg is > 1
        considers only profiles with valid RMS if required at file opening
        """
        try:
            var = self.hdf.select(var)
        except:
            print var + ' not present in file ' + self.filename + ' (version too old, maybe ?)'
            return None

        this_var = var[:].squeeze()
        if this_var.ndim == 1:

            if idx[0] is 0 and idx[1] is -1:
                data = this_var[:]
            else:
                data = this_var[idx[0]:idx[1]]

            if navg > 1:
                data = _vector_average(data, navg, missing=missing, valid=self.valid_rms_profiles)
        else:

            if idx[0] is 0 and idx[1] is -1:
                data = this_var[...]
            else:
                data = this_var[idx[0]:idx[1], :]

            if navg > 1:
                data = _array_average(data, navg, missing=missing, valid=self.valid_rms_profiles)

        var.endaccess()

        return data

    def valid_profiles(self, navg=30):
        """
        returns the percentage of valid profiles used to average
        over navg.
        Invalid profiles depend on the required RMS for pre-averaging filtering
        """
        if navg < 2:
            return self.valid_rms_profiles
        else:
            n = np.size(self.valid_rms_profiles, 0) / navg
            nprof = np.zeros(n)
            for i in np.arange(n):
                validslice = self.valid_rms_profiles[i * navg:i * navg + navg - 1]
                nprof[i] = np.sum(validslice)
            nprof = 100. * nprof / (navg - 1)

        return nprof

    def utc_time(self, navg=30, idx=(0, -1)):
        var = self.hdf.select('Profile_UTC_Time')
        time = var[:].squeeze()
        if idx[0] > 0 and idx[1] > -1:
            time = time[idx[0]:idx[1]]
        # if self.valid_rms_profiles is not None:
        #     time = time[self.valid_rms_profiles]
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

    def atb_std(self, navg=30, prof=None, idx=(0, -1)):
        atbstd = self._read_std('Total_Attenuated_Backscatter_532', navg, idx=idx)
        if prof:
            atbstd = self.atb[prof, :]
        return atbstd

    def atb(self, navg=30, prof=None, idx=(0, -1)):
        """
        Reads the Attenuated Total Backscatter 532nm from CALIOP file
        shape [nprof, nz]
        """
        atb = self._read_var('Total_Attenuated_Backscatter_532', navg, idx=idx)
        if prof:
            atb = atb[prof, :]
        return atb

    def parallel_rms_baseline(self, navg=30, prof=None, idx=(0, -1)):
        """
        Reads the Parallel RMS Baseline at 532nm from CALIOP file
        shape [nprof]
        JP Vernier utilise un seuil a 150 photons pour decider si un profil est bon ou pas
        units = counts
        """
        rms = self._read_var('Parallel_RMS_Baseline_532', navg, idx=idx, missing=-9999.)
        if prof:
            rms = rms[prof, :]
        return rms

    def ccu_532(self, navg=30, prof=None, idx=(0, -1)):
        au = self._read_var('Calibration_Constant_Uncertainty_532', navg, idx=idx)
        if prof:
            au = au[prof, :]
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

    def pressure_on_lidar_alt(self, navg=30, alt=lidar_alt, metalt=met_alt, idx=(0, -1)):
        """
        Reads the ancillary pressure field from CALIOP file, interpolated on
        CALIOP altitude levels.
        shape [nprof, nz]
        """
        p0 = self.pressure(navg=navg, idx=idx)
        p = _remap_y(p0, metalt, alt)

        return p

    def mol_on_lidar_alt(self, navg=30, prof=None, alt=lidar_alt, metalt=met_alt, idx=(0, -1)):
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

    def temperature_on_lidar_alt(self, navg=30, prof=None, alt=lidar_alt, metalt=met_alt, idx=(0, -1)):
        """
        Returns the ancillary temperature field in degC, interpolated on
        CALIOP altitude levels.
        """
        t0 = self.temperature(navg=navg, idx=idx)
        t = _remap_y(t0, metalt, alt)
        if prof:
            t = t[prof, :]
        return t

    def rh_on_lidar_alt(self, navg=30, prof=None, alt=lidar_alt, metalt=met_alt, idx=(0, -1)):
        """
        Returns the ancillary relative humidity field from the CALIOP file,
        interpolated on CALIOP altitude levels.
        """
        rh0 = self.rh(navg=navg, idx=idx)
        rh = _remap_y(rh0, metalt, alt)

        if prof:
            rh = rh[prof, :]
        return rh

    def scattering_ratio(self, navg=30, prof=None, alt=lidar_alt, metalt=met_alt, idx=(0, -1)):
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

    def _peek(self, lat, alt, data, latrange, dataname, vmin, vmax, datetime=False, ymin=0, ymax=25, cmap=None):

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
        """
        internal helping func to read a variable (1D or 2D) in HDF file
        """

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

        # the [:,1] is because level 2 vectors contain 3 values, that describe
        # the parameter for the first profile of the averaged section, the last profile,
        # and an average. The middle value is the average, that's what we use.

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
        """
        returns profile time (TAI)
        shape [nprof]
        """
        return self._read_var('Profile_Time', idx=idx)[:]

    def utc_time(self, idx=(0, -1)):
        """
        Returns utc time value (decimal time)
        """
        time = self._read_var('Profile_UTC_Time', idx=idx)[:, 1]
        return time

    def datetime(self, idx=(0, -1)):
        """
        Returns an array of datetime objects based on utc_time values
        """
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


    def datetime2(self, idx=(0, -1)):
        """
        Returns an array of datetime objects based on utc_time values
        this version is 5 times faster than the datetime function above. 
        Is it worth it ? not sure.
        """

        def _decdate_to_ymd(decdate):

            y = np.floor(decdate / 10000.)
            remainder = decdate - y * 10000
            m = np.floor(remainder / 100.)
            d = np.floor(remainder - m * 100)

            return y + 2000, m, d

        utc = self.utc_time(idx=idx)
        seconds_into_day = ((utc - np.floor(utc)) * 24. * 3600.)

        # first let's check if the orbit spans more than a single date
        y0, m0, d0 = _decdate_to_ymd(utc[0])
        y1, m1, d1 = _decdate_to_ymd(utc[-1])
        if d0 == d1:
            # orbit spans a single date
            # thus we can be a little faster
            datetimes = np.array(datetime.datetime(int(y0), int(m0), int(d0), 0, 0, 0)) + np.array(
                [datetime.timedelta(seconds=int(ss)) for ss in seconds_into_day])
        else:
            # orbits spans more than a day, we have to compute everything
            print 'multi date', y0, m0, d0, y1, m1, d1
            y, m, d = _decdate_to_ymd(utc)
            datetimes = [datetime.datetime(int(yy), int(mm), int(dd), 0, 0, 0) + datetime.timedelta(seconds=int(ss)) for
                         yy, mm, dd, ss in zip(y, m, d, seconds_into_day)]

        return np.array(datetimes)

    def statistics_532(self):

        var = self._read_var('Attenuated_Backscatter_Statistics_532')
        stats = dict()
        stats['min'] = var[:, 0::6]
        stats['max'] = var[:, 1::6]
        stats['mean'] = var[:, 2::6]
        stats['std'] = var[:, 3::6]
        stats['centroid'] = var[:, 4::6]
        stats['skewness'] = var[:, 5::6]

        return stats

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

    def layers_pressure(self, idx=(0, -1)):
        """
        returns layer pressure by profile
        units : hPa
        """

        ptop = self._read_var('Layer_Top_Pressure', idx=idx)
        pbase = self._read_var('Layer_Base_Pressure', idx=idx)
        return pbase, ptop

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

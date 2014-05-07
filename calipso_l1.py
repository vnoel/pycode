#!/usr/bin/env python
#encoding:utf-8

"""
CALIPSO Level 1 data file class

V. Noel 2008-2014
LMD/CNRS
"""

import os
import numpy as np
import numpy.ma as ma
import datetime
from calipso_base import _Cal, _vector_average, _array_std, _array_average, _remap_y
import netCDF4


static_path = os.path.dirname(__file__) + '/staticdata/'
# this is a hack while reading vdata segfaults pyhdf
# altitude vector before November 2007
lidar_alt_pretilt = np.loadtxt(static_path + 'lidaralt_pre2007.asc')
# altitude vector after November 2007
lidar_alt_posttilt = np.loadtxt(static_path + 'lidaralt_post2007.asc')
# this has a higher chance of being the correct altitude
lidar_alt = lidar_alt_posttilt

# FIXME : need to verify if the same thing is needed for met_alt
met_alt = np.loadtxt(static_path + 'metalt.asc')

# maximum molecular atb for normalization
atb_max = {'ZN': 1e-4, 'ZD': 1}

# maximum integrated atb for calibration in the 26-28 km range
iatb_bounds = {'ZN': [1e-5, 8e-5], 'ZD': [-8e-3, 8e-3]}


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
        if self.date < datetime.datetime(2007,12,1):
            self.lidar_alt = lidar_alt_pretilt
        else:
            self.lidar_alt = lidar_alt_posttilt
        self.met_alt = met_alt
        if max_rms is not None:
            rms = self.parallel_rms_baseline(navg=1)
            self.valid_rms_profiles = (rms < max_rms)

    def _read_sds(self, varname):
        try:
            hdfvar = self.hdf.select(varname)
        except:
            print('Cannot read ' + varname)
            return None
            
        var = hdfvar[:].squeeze()
        hdfvar.endaccess()
        return var

    def _read_std(self, varname, navg):
        """
        Reads a variable in an hdf file, and computes the standard deviation
        of the variable over navg profiles
        """

        var = self._read_sds(varname)
        if var.ndim == 1:
            print('sorry, ndim=1 not implemented in _read_std')
            return None
        data = var[...]
        data = _array_std(data, navg, valid=self.valid_rms_profiles)
        
        return data

    def _read_var(self, varname, navg, idx=None, missing=None):
        """
        Read a variable in an hdf file, averaging the data if navg is > 1
        considers only profiles with valid RMS if required at file opening
        """
        
        if navg==0:
            return []
        
        var = self._read_sds(varname)
        if navg > np.size(var, 0):
            return None
        
        if idx is None:
            data = var[...]
        else:
            i0, i1 = idx
            if var.ndim == 1:
                data = var[i0:i1]
            elif var.ndim == 2:
                data = var[i0:i1,:]
        
        if navg > 1:
            avg_function = {1:_vector_average, 2:_array_average}[var.ndim]
            data = avg_function(data, navg, missing=missing, valid=self.valid_rms_profiles)

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

    def utc_time(self, navg=30, idx=None):
        if navg==0:
            return []
        time = self._read_var('Profile_UTC_Time', navg=1)
        if navg > np.size(time, 0):
            return None
        if idx is not None:
            time = time[idx[0]:idx[1]]
        if time is not None and navg > 1:
            n = np.size(time, 0) / navg
            time2 = np.zeros(n)
            for i in np.arange(n):
                time2[i] = time[i * navg]
            time = time2

        return time

    def datetimes(self, navg=30):
        if navg==0:
            return []
        time = self.time(navg=navg)
        datetimes = netCDF4.num2date(time, units='seconds since 1993-01-01')
        return np.array(datetimes)

    def time(self, navg=30, idx=None):
        """
        Reads time for CALIOP profiles, averaged on navg
        shape [nprof]
        Example:
            time = c.time(navg=15)
        Units: seconds from January 1, 1993
        """
        if navg==0:
            return []
        time = self._read_var('Profile_Time', navg=1, idx=idx)
        if navg > len(time):
            return None
        if time is not None and navg > 1:
            n = np.size(time, 0) / navg
            time2 = np.zeros(n)
            for i in np.arange(n):
                time2[i] = time[i * navg]
            time = time2

        return time

    def altitude(self):
        '''
        Reads altitude levels from CALIOP metadata.
        '''
        from pyhdf import VS
        from pyhdf import HDF
        
        hdffile = HDF.HDF(self.filename, HDF.HC.READ)
        vs = hdffile.vstart()
        meta_id = vs.find('metadata')
        vd = vs.attach('metadata', write=0) # same as vs.attach(meta_id, write=0)
        alt = vd.field('Lidar_Data_Altitudes')  # segfaults
        
        # workaround to get the altitude metadata field : use the hdf dump tool
        # hdp dumpvd -d -n metadata -f Lidar_Data_Altitudes HDF_FILENAME
        
        return alt

    def coords(self, navg=30, idx=None):
        """
        Reads coordinates for CALIOP profiles, averaged on navg
        shape [nprof]
        Example:
            lon, lat = c.coords(navg=15)
        """
        lat = self._read_var('Latitude', navg, idx=idx)
        lon = self._read_var('Longitude', navg, idx=idx)

        return lon, lat

    def surface_elevation(self, navg=30, idx=None):
        """
        Reads the surface elevation from CALIOP file
        shape [nprof]
        """
        elev = self._read_var('Surface_Elevation', navg, idx=idx)
        return elev

    def perp(self, navg=30, idx=None):
        """
        Reads the perpendicular signal from CALIOP file
        shape [nprof, nz]
        """
        perp = self._read_var('Perpendicular_Attenuated_Backscatter_532',
                              navg, idx=idx)
        return perp

    def atb_std(self, navg=30):
        """
        Standard deviation from the Attenuated Total Backscatter 532nm from CALIOP file
        shape [nprof/navg, nz]
        """
        atbstd = self._read_std('Total_Attenuated_Backscatter_532', navg)
        return atbstd

    def atb(self, navg=30, idx=None):
        """
        Reads the Attenuated Total Backscatter 532nm from CALIOP file
        shape [nprof, nz]
        """
        atb = self._read_var('Total_Attenuated_Backscatter_532', navg, idx=idx)
        return atb

    def parallel_rms_baseline(self, navg=30, idx=None):
        """
        Reads the Parallel RMS Baseline at 532nm from CALIOP file
        shape [nprof]
        JP Vernier utilise un seuil a 150 photons pour decider si un profil est bon ou pas
        units = counts
        """
        rms = self._read_var('Parallel_RMS_Baseline_532', navg, idx=idx, missing=-9999.)
        return rms

    def ccu_532(self, navg=30, idx=None):
        au = self._read_var('Calibration_Constant_Uncertainty_532', navg, idx=idx)
        return au

    def atb1064(self, navg=30, idx=None):
        """
        Reads the Attenuated Total Backscatter 1064nm from CALIOP file
        """
        atb = self._read_var('Attenuated_Backscatter_1064', navg, idx=idx)
        return atb

    def pressure(self, navg=30, idx=None):
        """
        Reads the ancillary pressure field from CALIOP file
        shape [nprof, nlevels]
        """
        p = self._read_var("Pressure", navg, idx=idx)
        return p

    def mol(self, navg=30, idx=None):
        """
        Reads the ancillary molecular number density from CALIOP file
        shape [nprof, nlevels]
        """
        mol = self._read_var('Molecular_Number_Density', navg, idx=idx)
        return mol

    def temperature(self, navg=30, idx=None):
        """
        Reads the ancillary temperature field in degC from CALIOP file
        shape [nprof, nlevels]
        """
        t = self._read_var('Temperature', navg, idx=idx)
        return t

    def rh(self, navg=30, idx=None):
        """
        Reads the ancillary relative humidity field from CALIOP file
        shape [nprof, nlevels]
        """
        rh = self._read_var('Relative_Humidity', navg, idx=idx)
        return rh

    def pressure_on_lidar_alt(self, navg=30, alt=lidar_alt, metalt=met_alt, idx=None):
        """
        Reads the ancillary pressure field from CALIOP file, interpolated on
        CALIOP altitude levels.
        shape [nprof, nz]
        """
        p0 = self.pressure(navg=navg, idx=idx)
        p = _remap_y(p0, metalt, alt)

        return p

    def mol_on_lidar_alt(self, navg=30, alt=lidar_alt, metalt=met_alt, idx=None):
        """
        Reads the ancillary molecular number density from CALIOP file,
        interpolated on
        CALIOP altitude levels.
        shape [nprof, nz]
        """
        mol0 = self.mol(navg=navg, idx=idx)
        mol = _remap_y(mol0, metalt, alt)

        return mol

    def mol_calibration_coef(self, mol=None, atb=None, navg=30, 
                             alt=lidar_alt, metalt=met_alt, idx=None, navgh=50, zmin=30,
                             zmax=34):
        """
        Returns the molecular calibration coefficient, computed from
        atb 532 nm and molecular density profiles,
        averaged between zmin and zmax km vertically, and [i-navgh:i+navgh]
        profiles horizontally using a moving average.
        shape [nprof]
        """
        if mol is None and atb is None:
            mol = self.mol_on_lidar_alt(navg=navg, alt=alt,
                                        metalt=metalt, idx=idx)
            atb = self.atb(navg=navg, idx=idx)

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

    def mol_on_lidar_alt_calibrated(self, navg=30, alt=lidar_alt,
                                    navgh=50, metalt=met_alt, idx=None, zcal=(30, 34), atb=None):
        """
        Returns an estimate of the molecular backscatter at 532 nm, computed
        from the molecular number density calibrated on
        the attenuated total backscatter at 532 nm (both from the CALIOP file),
        using the calibration coefficient returned
        from mol_calibration_coef().
        Shape: [nprof, nz]
        
        atb can be passed as an argument if it's been read and averaged already.
        """
        mol = self.mol_on_lidar_alt(navg=navg, alt=alt,
                                    metalt=metalt, idx=idx)
        if atb is None:
            atb = self.atb(navg=navg, idx=idx)

        coef = self.mol_calibration_coef(mol=mol, atb=atb, zmin=zcal[0],
                                         zmax=zcal[1], navgh=navgh)

        # x * y is equivalent to x[:,i] * y[i]
        # mol is [i,:], so we need to rotate it twice
        mol = mol.T
        mol *= coef
        mol = mol.T

        return mol

    def temperature_on_lidar_alt(self, navg=30, alt=lidar_alt, metalt=met_alt, idx=None):
        """
        Returns the ancillary temperature field in degC, interpolated on
        CALIOP altitude levels.
        """
        t0 = self.temperature(navg=navg, idx=idx)
        t = _remap_y(t0, metalt, alt)
        return t

    def rh_on_lidar_alt(self, navg=30, alt=lidar_alt, metalt=met_alt, idx=None):
        """
        Returns the ancillary relative humidity field from the CALIOP file,
        interpolated on CALIOP altitude levels.
        """
        rh0 = self.rh(navg=navg, idx=idx)
        rh = _remap_y(rh0, metalt, alt)

        return rh

    def scattering_ratio(self, navg=30, alt=None, metalt=None, idx=None):
        """
        Returns the scattering ratio, i.e. the ratio between the attenuated
        total backscatter and the molecular
        backscatter, both at 532 nm, both read from the CALIOP file.
        """
        mol_calib = self.mol_on_lidar_alt_calibrated(navg=navg,
                                                     alt=alt, metalt=metalt, idx=idx)
        atb = self.atb(navg=navg, idx=idx)
        sr = atb / mol_calib
        return sr

    def tropopause_height(self, navg=30, idx=None):
        """
        Reads the ancillary tropopause height, in km, from the CALIOP file.
        shape [nprof]
        """
        tropoz = self._read_var('Tropopause_Height', navg, idx=idx)
        return tropoz

    def tropopause_temperature(self, navg=30, idx=None):
        """
        Reads the ancillary tropopause temperature, in degC,
        from the CALIOP file.
        shape [nprof]
        """
        tropot = self._read_var('Tropopause_Temperature', navg, idx=idx)
        return tropot

    # some display utility functions

    def _peek(self, lat, alt, data, latrange, dataname, vmin, vmax, axis_datetime=False, ymin=0, ymax=25, cmap=None):

        import niceplots

        import matplotlib.pyplot as plt
        from pylab import get_cmap

        print('Showing ' + dataname)
        print('Lat range : ', latrange)

        if cmap is None:
            cmap = get_cmap('gist_stern_r')

        plt.ioff()
        fig = plt.figure(figsize=(20, 6))
        ax = plt.gca()
        plt.pcolormesh(lat, alt, data.T, vmin=vmin, vmax=vmax, cmap=cmap)
        plt.xlim(latrange[0], latrange[1])
        plt.ylim(ymin, ymax)
        if axis_datetime:
            niceplots.axis_set_date_format(ax, format='%H:%M')
            fig.autofmt_xdate()
        else:
            plt.xlabel('Latitude')

        plt.colorbar().set_label(dataname)
        plt.ylabel('Altitude [km]')
        plt.title(dataname + ' - CALIOP ' + str(self.date))

    def peek_atb_time(self, navg=30, mintime=None, maxtime=None):

        from pylab import get_cmap

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
                   'Coefficient de retrodiffusion [log10]', -3.5, -2, axis_datetime=True,
                   cmap=get_cmap('gist_stern_r'))

    def peek_cr_time(self, navg=30, mintime=None, maxtime=None):

        from pylab import get_cmap

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
                   'Volumic Attenuated Color Ratio', 0, 1.4, axis_datetime=True,
                   cmap=get_cmap('jet'))

    def peek_psc_time(self, navg=30, mintime=None, maxtime=None):
        atb = self.atb(navg=navg)
        time = self.datetimes(navg=navg)
        import matplotlib.dates as mdates

        print(np.min(time), np.max(time))
        print(mintime, maxtime)

        numtime = mdates.date2num(time)

        if mintime is None:
            mintime = np.min(numtime)
        else:
            mintime = mdates.date2num(mintime)
            if mintime < np.min(numtime):
                print('Error : mintime < min(time)')

        if maxtime is None:
            maxtime = np.max(numtime)
            if maxtime > np.max(numtime):
                print('Error : maxtime > max(time)')
        else:
            maxtime = mdates.date2num(maxtime)

        tropo = self.tropopause_height(navg=navg)
        tropo[tropo > 13] = 11
        for i, z in enumerate(tropo):
            # I don't remember why this exists
            idx = (lidar_alt < (z + 3))
            atb[i, idx] *= 0.5
        atb[atb < 0] = 1e-6

        self._peek(numtime, lidar_alt, np.log10(atb), [mintime, maxtime],
                   'Coefficient de retrodiffusion (log10)', -3.5, -2.25, axis_datetime=True,
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
        print('Averaging every %d profiles = %d shown profiles' % (navg, len(lat)))
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
        print('Averaging every %d profiles = %d shown profiles' % (navg, len(lat)))

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
        print('Averaging every %d profiles = %d shown profiles' % (navg, len(lat)))

        self._peek(lat, lidar_alt, cr, latrange,
                   'Volumic Color Ratio', 0.2, 1.)  # , cmap=get_cmap('jet'))


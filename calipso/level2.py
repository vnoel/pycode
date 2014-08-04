#!/usr/bin/env python
#encoding:utf-8

"""
CALIPSO Level 2 data file class

V. Noel 2008-2014
LMD/CNRS
"""

import numpy as np
import datetime
from calipso_hdf import _Cal


class Cal2(_Cal):
    """
    Class to process CALIOP Level 2 files.
    No averaging is possible here given the qualitative nature of variables. (can't average altitudes etc.)
    example use:
        
        from calipso import Cal2
        
        c = Cal2('CAL_LID_L2_05kmCLay-Prov-V3-01.2010-12-31T01-37-30ZN.hdf')
        lon, lat = c.coords()
        nl, base, top = c.layers()
        ...
        c.close()
        
    """

    def __init__(self, filename):
        _Cal.__init__(self, filename)
        lat = self._read_var('Latitude')
        # identify 333m or more level 2 files
        if lat.shape[1] == 1:
            self.havg = 0.333
            self.iavg = 0
        else:
            self.havg = 5.
            self.iavg = 1

    def coords(self, idx=None):
        """
        Returns longitude and latitude for profiles.
        shape [nprof]
        # the [:,1] is because level 2 vectors contain 3 values, that describe
        # the parameter for the first profile of the averaged section, the last profile,
        # and an average. The middle value is the average, that's what we use.
        """
        lat = self._read_var('Latitude', idx=idx)[:, self.iavg]
        lon = self._read_var('Longitude', idx=idx)[:, self.iavg]
        return lon, lat
    
    def coords_bounds(self):
        """
        returns the lat and lon boundaries for each 5-km profile.
        """
        if self.havg < 1.:
            raise BaseException('333m file == no boundaries')
        lat = self._read_var('Latitude')
        lat0, lat1 = lat[:, 0], lat[:,-1]
        lon = self._read_var('Longitude')
        lon0, lon1 = lon[:, 0], lon[:,-1]
        return (lon0, lat0), (lon1, lat1)

    def time(self, idx=None):
        """
        returns profile time (TAI)
        shape [nprof]
        """
        return self._read_var('Profile_Time', idx=idx)[:][:, self.iavg]

    def time_bounds(self):
        '''
        returns the time boundaries for each 5-km profile.
        '''
        if self.havg < 1:
            raise BaseException('333m file == no boundaries')
        time = self._read_var('Profile_Time')
        time0, time1 = time[:,0], time[:,-1]
        return time0, time1
    

    def utc_time(self, idx=None):
        """
        Returns utc time value (decimal time)
        """
        time = self._read_var('Profile_UTC_Time', idx=idx)[:, self.iavg]
        return time

    def datetime(self, idx=None):
        """
        Returns an array of datetime objects based on utc_time values
        """
        time = self.time(navg=navg)
        datetimes = netCDF4.num2date(time, units='seconds since 1993-01-01')
        return np.array(datetimes)

    def datetime2(self, idx=None):
        """
        Returns an array of datetime objects based on utc_time values
        this version is 5 times faster than the datetime function above. 
        Is it worth it ? not sure.
        """

        def _decdate_to_ymd(decdate):

            year = np.floor(decdate / 10000.)
            remainder = decdate - year * 10000
            month = np.floor(remainder / 100.)
            day = np.floor(remainder - month * 100)

            return year + 2000, month, day

        utc = self.utc_time(idx=idx)
        seconds_into_day = ((utc - np.floor(utc)) * 24. * 3600.)

        # first let's check if the orbit spans more than a single date
        y0, m0, d0 = _decdate_to_ymd(utc[0])
        y1, m1, d1 = _decdate_to_ymd(utc[-1])
        if d0 == d1:
            # orbit spans a single date
            # we can be a little faster
            datetimes = np.array(datetime.datetime(int(y0), int(m0), int(d0), 0, 0, 0)) + np.array(
                [datetime.timedelta(seconds=int(ss)) for ss in seconds_into_day])
        else:
            # orbits spans more than a day, we have to compute everything
            print('multi date', y0, m0, d0, y1, m1, d1)
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

    def off_nadir_angle(self, idx=None):
        """
        Returns the off-nadir-angle, in deg, for profiles.
        shape [nprof]
        """
        return self._read_var('Off_Nadir_Angle', idx=idx)

    def tropopause_height(self, idx=None):
        """
        Returns the ancillary tropopause height, in km, for profiles.
        shape [nprof]
        """
        return self._read_var('Tropopause_Height', idx=idx)[:, 0]

    def tropopause_temperature(self, idx=None):
        """
        Returns the ancillary tropopause temperature, in degC, for profiles
        shape [nprof]
        """
        return self._read_var('Tropopause_Temperature', idx=idx)

    def dem_surface_elevation(self, idx=None):
        """
        Returns the ancillary surface elevation (from digital elevation model)
        in km, for profiles
        shape [nprof]
        """
        return self._read_var('DEM_Surface_Elevation', idx=idx)

    def lidar_surface_elevation(self, idx=None):
        """
        returns 8 values per profile
        min, max, mean, stddev for upper boundary of the surface echo.
        min, max, mean, stddev for lower boundary of the surface echo.
        en kilometres 
        """

        return self._read_var('Lidar_Surface_Elevation', idx=idx)

    def IGBP_Surface_Type(self, idx=None):
        """
        IGBP_Surface_Type:format = "Int_8" ;
        IGBP_Surface_Type:valid_range = "1....18" ;
        IGBP_Surface_Type:fillvalue = '\367' ;
        IGBP_Surface_Type:range_value = "evergreen needleleaf forest, evergreen broadleaf forest, deciduous needleleaf forest, deciduous broadleaf forest, mixed forest, closed shrublands, open shrublands,woody savannas, savannas, grasslands, permanent wetlands, croplands, urban and built-up,cropland/natural vegetation mosaic, snow and ice, barren or sparsely vegetated, water bodies, tundra" ;
        """
        return np.squeeze(self._read_var('IGBP_Surface_Type', idx=idx))
        

    # Layer info

    def layers(self, idx=None):
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

    def layers_pressure(self, idx=None):
        """
        returns layer pressure by profile
        units : hPa
        """

        ptop = self._read_var('Layer_Top_Pressure', idx=idx)
        pbase = self._read_var('Layer_Base_Pressure', idx=idx)
        return pbase, ptop

    def layers_number(self, idx=None):
        """
        Returns the number of layer found by profile
        shape [nprof]
        """
        return self._read_var('Number_Layers_Found', idx=idx)[:, 0]

    def midlayer_temperature(self, idx=None):
        """
        Returns the midlayer temperature by layer, in km
        shape [nprof, nlaymax]
        """
        return self._read_var('Midlayer_Temperature', idx=idx)

    def flag(self, idx=None):
        """
        Returns the feature classification flag by layer
        shape [nprof, nlaymax]
        cf https://eosweb.larc.nasa.gov/sites/default/files/project/calipso/quality_summaries/CALIOP_L2VFMProducts_3.01.pdf
        """
        return self._read_var('Feature_Classification_Flags', idx=idx)

    def layer_type(self, idx=None):
        """
        Returns the layer type from the feature classification flag
        shape [nprof, nlaymax]
        
        0 = invalid (bad or missing data) 
        1 = "clear air"
        2 = cloud
        3 = aerosol
        4 = stratospheric feature
        5 = surface
        6 = subsurface
        7 = no signal (totally attenuated)
        
        """
        f = self.flag(idx=idx)
        # type flag : bits 1 to 3
        typeflag = (f & 7)
        return typeflag

    def layer_subtype(self, idx=None):
        """
        Returs the layer subtype, as identified from the feature
        classification flag
        shape [nprof, nlaymax]
        
        for clouds (feature type == layer_type == 2)
        0 = low overcast, transparent
        1 = low overcast, opaque
        2 = transition stratocumulus
        3 = low, broken cumulus
        4 = altocumulus (transparent) 5 = altostratus (opaque)
        6 = cirrus (transparent)
        7 = deep convective (opaque)
        """
        f = self.flag(idx=idx)
        # subtype flag : bits 10 to 12
        subtype = (f & 3584) >> 9
        return subtype

    def layer_type_qa(self, idx=None):
        """
        Returs the quality flag for the layer type, as identified from the
        feature classification flag
        shape [nprof, nlaymax]
        """
        f = self.flag(idx=idx)
        typeflag_qa = (f & 24) >> 3
        return typeflag_qa

    def phase(self, idx=None):
        """
        Returs the layer thermodynamical phase, as identified from the
        feature classification flag
        shape [nprof, nlaymax]
        
        0 = unknown / not determined 1 = randomly oriented ice
        2 = water
        3 = horizontally oriented ice
        """
        f = self.flag(idx=idx)
        # 96 = 0b1100000, bits 6 to 7
        cloudflag = (f & 96) >> 5
        return cloudflag

    def phase_qa(self, idx=None):
        """
        Returns the quality flag for the layer thermodynamical phase,
        as identified from the feature classification flag
        shape [nprof, nlaymax]
        
        0 = none
        1 = low
        2 = medium 3 = high
        """

        f = self.flag(idx=idx)
        cloudflag_qa = (f & 384) >> 7
        return cloudflag_qa

    def opacity_flag(self, idx=None):
        """
        Returns the opacity flag by layer.
        shape [nprof, nlaymax]
        """

        return self._read_var('Opacity_Flag', idx=idx)

    def horizontal_averaging(self, idx=None):
        """
        Returns the horizontal averaging needed to detect a given layer.
        shape [nprof, nlaymax]
        """
        return self._read_var('Horizontal_Averaging', idx=idx)

    def iatb532(self, idx=None):
        """
        Returns the integrated attenuated total backscatter at 532 nm
        along the layer thickness.
        shape [nprof, nlaymax]
        """
        return self._read_var('Integrated_Attenuated_Backscatter_532', idx=idx)

    def ivdp(self, idx=None):
        """
        Returns the volumic depolarization ratio for the entire layer
        thickness, obtained by the
        ratio of integrated perpendicular backscatter on the integrated
        parallel backscatter at 532 nm
        shape [nprof, nlaymax]
        """
        return self._read_var('Integrated_Volume_Depolarization_Ratio', idx=idx)

    def ipdp(self, idx=None):
        """
        Returns the particulate depolarization ratio for the entire
        layer thickness, i.e. the volumic
        depolarization ratio corrected to remove its molecular component.
        shape [nprof, nlaymax]
        """
        return self._read_var('Integrated_Particulate_Depolarization_Ratio', idx=idx)

    def icr(self, idx=None):
        """
        Returns the integrated attenuated color ratio for the entire
        layer thickness.
        shape [nprof, nlaymax]
        """
        return self._read_var('Integrated_Attenuated_Total_Color_Ratio', idx=idx)

    def ipcr(self, idx=None):
        """
        Returns the integrated color ratio for the entire layer thickness
        shape [nprof, nlaymax]
        """
        return self._read_var('Integrated_Particulate_Color_Ratio', idx=idx)

    def od(self, idx=None):
        """
        Returns the optical depth found by layer.
        shape [nprof, nlaymax]
        """
        return self._read_var('Feature_Optical_Depth_532', idx=idx)

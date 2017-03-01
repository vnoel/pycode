#!/usr/bin/env python
#encoding:utf-8

# Created by V. Noel [LMD/CNRS] on 2016-01-11

import numpy as np
import h5py
from datetime import datetime, timedelta

def datetime_from_utc(utc_date, utc_time):
 
    dt = []
    if len(utc_time.shape) > 1:
        utc_time = utc_time[:,1]
    for date, time in zip(utc_date.values, utc_time.values):
        datestr = str(date)
        dt.append(datetime.strptime(datestr, '%Y%m%d.0') + timedelta(seconds=(time*24.*3600.)))

    dt = np.array(dt)

    return dt   


class cats_base(object):

    def __init__(self, file):
        
        self.h5 = h5py.File(file, 'r')
        self.geolocation = self.h5['geolocation']


    def close(self):
        
        self.h5.close()


    def coords(self):
        
        lat = self.geolocation['CATS_Fore_FOV_Latitude'][:]
        lon = self.geolocation['CATS_Fore_FOV_Longitude'][:]
        return lat, lon


class cats_l2(cats_base):
    
    def __init__(self, file):
        cats_base.__init__(self, file)
        self.layer_descriptor = self.h5['layer_descriptor']
        
    def time_utc(self):

        # day
        utc_date = self.layer_descriptor['Profile_UTC_Date'][:]
        # fraction of day
        utc_time = self.layer_descriptor['Profile_UTC_Time'][:]
        dt = datetime_from_utc(utc_date, utc_time)

        return dt        
        
    def layer_top(self):
        # -999.99 when empty
        # shape = [nprof, nlayer=10]
        top = self.layer_descriptor['Layer_Top_Altitude_Fore_FOV'][:]
        return top
    
    def layer_base(self):
        # -999.99 when empty
        base = self.layer_descriptor['Layer_Base_Altitude_Fore_FOV'][:]
        return base
        
    def layer_type(self):
        # same shape as layer_top
        # 0 = empty
        # 1 = cloud
        # 2 = undetermined
        # 3 = aerosol        
        type = self.layer_descriptor['Feature_Type_Fore_FOV'][:]
        return type
        
    def layer_opacity(self):
        # same shape as layer_top
        opacity = self.layer_descriptor['Opacity_Fore_FOV'][:]
        return opacity
        
    def nl(self):
        nlayers = self.layer_descriptor['Number_Layers_Fore_FOV'][:]
        return nlayers
        
    def layers(self):
        nl, type, top, base = self.nl(), self.layer_type(), self.layer_top(), self.layer_base()
        return nl, type, base, top
        

class cats_l1(cats_base):
    
    def __init__(self, file):
        cats_base.__init__(self, file)
        self.ancillary_data = self.h5['ancillary_data']
        self.profile_attributes = self.h5['profile_attributes']
        
                
    def time_utc(self):
        
        # day
        utc_date = self.profile_attributes['Profile_UTC_Date'][:]
        # fraction of day
        utc_time = self.profile_attributes['Profile_UTC_Time'][:]
        dt = datetime_from_utc(utc_date, utc_time)

        return dt
        
        
    def atb532(self):
        
        atb = self.h5['channel_532']['FFOV']['Total_Attenuated_Backscatter532_Fore_FOV'][:]
        #atbmol = self.h5['ancillary_data']['MET']['Molecular_Backscatter532'][:]
        #atb = atb + atbmol
        return atb


    def atb1064(self):
        
        atb = self.h5['channel_1064']['FFOV']['Total_Attenuated_Backscatter1064_Fore_FOV'][:]
        #atbmol = self.h5['ancillary_data']['MET']['Molecular_Backscatter532'][:]
        #atb = atb + atbmol
        return atb
        
        
    def altitude(self):
        
        altitude = self.h5['metadata_parameters']['Bin_Altitude_Array'][:]
        return altitude

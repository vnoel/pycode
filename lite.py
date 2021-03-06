#!/usr/bin/env python
#encoding:utf-8

# Created by VNoel on 2014-08-11

'''
Class to read LITE lidar data.

LITE = Lidar In-space Technology Experiment
http://www-lite.larc.nasa.gov

Based on code from M. Reverdy (ESA/FX-Conseil)

Example usage:

    >>> import lite
    >>> l = lite.LITE('LITE_L1_19940910_164558_164706')
    >>> print l.rawdata['latitude']
    >>> l.describe()

The script can also be called on a LITE data file for quick checks:

    ~> ./lite.py LITE_L1_19940910_164558_164706 2
    Profile :  2
    Version number:  1 0
    Orbit number:  13
    ID number:  5001587
    Date:  253 16 45 58 8
    MetDate:  0 18 23 3 13
    lat, lon:  45.1957 120.065
    Number of profiles:  (572, 3000)

After reading, all contents of the LITE file are available in l.rawdata[].    
The photon counts are in fields l.rawdata['profile355'], l.rawdata['profile532'] and l.rawdata['profile064']

Backscatter fields ("atb") are initially empty. To fill them, they need to be calibrated:

    >>> l.calibrate_atb_to(532, 1e-4, between_altitudes=[20,22])
    
This will calibrate the atb at 532nm by bringing the average value between 20 and 22 km to 1e-4.
Instead of a single value, a vector of calibration values can be passed (length = number of profiles in the file).

'''

import numpy as np
from datetime import datetime, timedelta

header = np.dtype( [ 
                ('syncvalue', '>i2'), 
                ('majorversionnumber', '>i1'),
                ('minorversionnumber', '>i1'),
                ('datatakeid', '>7i1'),
                ('orbitnumber', '>i1'),
                ('idnumber', '>i4'),
                ('gmtday', '>i2'),
                ('gmthour', '>i1'),
                ('gmtmin', '>i1'),
                ('gmtsec', '>i1'),
                ('gmthund', '>i1'),
                ('metday', '>i2'),
                ('methour', '>i1'),
                ('metmin', '>i1'),
                ('metsec', '>i1'),
                ('methund', '>i1'),
                ('latitude', '>f4'),
                ('longitude', '>f4'),
                ('shufflealtitude', '>f4'),
                ('offnadirangle', '>f4'),
                ('digitizerondelay', '>f4'),
                ('datatakemode', '>i1'),
                ('specialopsmode', '>i1'),
                ('profilevalidstatus', '>i1'),
                ('landwaterflag', '>i1'),
                ('surfelevfootprint', '>f4'),
                ('metdataalts', '>18f4'),
                ('mettemps', '>18f4'),
                ('alttropopause', '>f4'),       # in km
                ('temptropopause', '>f4'),      # in Kelvins
                ('laserselected', '>i1'),
                ('baalignmentstatus', '>i1'),
                ('isdbstatus', '>i1'),
                ('baddatastatus', '>i1'),
                ('aodatastatus', '>i1'),
                ('motorinmotion', '>i1'),
                ('aperwheelstatus', '>i1'),
                ('backgroundmongain', '>i1'),
                
                ('surfacemode355', '>i1'),
                ('dbattenuation355', '>i1'),
                ('numbersatabovesurf355', '>i2'),
                ('highestsatsample355', '>f4'),
                ('numberunderflow355', '>i2'),
                ('filterstatus355', '>i1'),
                ('calibrationstatus355', '>i1'),
                ('calibrationfactor355', '>f4'),
                ('baselinerippleremvd355', '>i1'),
                ('oscillationremoved355', '>i1'),
                ('backgroundvalue355', '>i1'),
                ('highvoltage355enabled', '>i1'),
                ('highvoltage355', '>f4'),
                ('energymonitor355', '>f4'),
                ('pmtgain355', '>f4'),
                ('baselinesubmethod355', '>i1'),
                ('subregionunderflow355', '>i1'),
                ('anomalousprof355', '>i1'),
                ('fillbyte1', 'i1'),

                ('surfacemode532', '>i1'),
                ('dbattenuation532', '>i1'),
                ('numbersatabovesurf532', '>i2'),
                ('highestsatsample532', '>f4'),
                ('numberunderflow532', '>i2'),
                ('filterstatus532', '>i1'),
                ('calibrationstatus532', '>i1'),
                ('calibrationfactor532', '>f4'),
                ('baselinerippleremvd532', '>i1'),
                ('oscillationremoved532', '>i1'),
                ('backgroundvalue532', '>i1'),
                ('highvoltage532enabled', '>i1'),
                ('highvoltage532', '>f4'),
                ('energymonitor532', '>f4'),
                ('pmtgain532', '>f4'),
                ('baselinesubmethod532', '>i1'),
                ('subregionunderflow532', '>i1'),
                ('anomalousprof532', '>i1'),
                ('fillbyte2', 'i1'),
                
                ('surfacemode064', '>i1'),
                ('dbattenuation064', '>i1'),
                ('numbersatabovesurf064', '>i2'),
                ('highestsatsample064', '>f4'),
                ('numberunderflow064', '>i2'),
                ('filterstatus064', '>i1'),
                ('calibrationstatus064', '>i1'),
                ('calibrationfactor064', '>f4'),
                ('baselinerippleremvd064', '>i1'),
                ('oscillationremoved064', '>i1'),
                ('backgroundvalue064', '>i1'),
                ('highvoltage064enabled', '>i1'),
                ('highvoltage064', '>f4'),
                ('energymonitor064', '>f4'),
                ('pmtgain064', '>f4'),
                ('baselinesubmethod064', '>i1'),
                ('subregionunderflow064', '>i1'),
                ('anomalousprof064', '>i1'),
                ('fillbyte3', 'i1'),
                
                ('timeedsinthour', '>i1'),
                ('timeedsintmin', '>i1'),
                ('timeedsintsec', '>i1'),
                ('timeedsinthund', '>i1'),
                ('level0fileidnumber', '>i1'),
                ('level0fileidletter', '>i1'),
                ('reserved', '>6i1'),
                ('highvoltage355cmd', '>f4'),
                ('highvoltage532cmd', '>f4'),
                ('reserved2', '>i4'),
                ('b0_355', '>f4'),
                ('b0_532', '>f4'),
                ('b0_064', '>f4'),
                
                ('outofrng355abv40', '>i1'),
                ('outofrng532abv40', '>i1'),
                ('outofrng064abv40', '>i1'),
                ('outofrange355', '>375i1'),
                ('outofrange532', '>375i1'),
                ('outofrange064', '>375i1'),
                
                ('top355', '>i2'),
                ('bottom355', '>i2'),
                ('top532', '>i2'),
                ('bottom532', '>i2'),
                ('top064', '>i2'),
                ('bottom064', '>i2'),
                
                ('profile355', '>3000f4'),
                ('profile532', '>3000f4'),
                ('profile064', '>3000f4')
                
                ] )
    
wv_keys = {355:'profile355', 532:'profile532', 1064:'profile064'}

class LITE(object):
    
    def __init__(self, filename, altitude_bottom_up=True):
        
        fid = open(filename, 'r')
        self.rawdata = np.fromfile(fid, dtype=header, count=-1)
        fid.close()
        self.altitude = np.linspace(40, -4.985, 3000)
        self.nprof = self.rawdata['profile355'].shape[0]
        self.latitude = self.rawdata['latitude']
        self.longitude = self.rawdata['longitude']

        if altitude_bottom_up:
            self.altitude = self.altitude[::-1]
            for field in 'profile355', 'profile532', 'profile064':
                self.rawdata[field] = self.rawdata[field][:,::-1]

        self.datetimes = []
        for i in np.r_[0:self.nprof]:
            date = datetime(1994,1,1) + timedelta(days=self.rawdata['gmtday'][i]-1)
            fulldate = datetime(date.year, date.month, date.day, self.rawdata['gmthour'][i], self.rawdata['gmtmin'][i], self.rawdata['gmtsec'][i], self.rawdata['gmthund'][i]*10000)
            self.datetimes.append(fulldate)

        self.atb = dict()
        self.atb[355] = None
        self.atb[532] = None
        self.atb[1064] = None            
    
    def calibrate_atb_to(self, wv, mol_atb, between_altitudes=[20,24], horiz_avg=10):

        '''
        calibrate the photon counts as ATB, by bringing the average beween altcal[0] and altcal[1] to the values in mol_atb
        mol_atb is either a single value or a vector with length nprof
        
        horiz_avg is the optional window size for the rolling average to smooth out the calibration factor.
        
        the function returns the calibration factor.
        '''
        
        idx = (self.altitude >= between_altitudes[0]) & (self.altitude < between_altitudes[1])
        reference = np.mean(self.rawdata[wv_keys[wv]][:,idx], axis=1)
        if horiz_avg > 1:
            reference = _rolling_average_padded(reference, horiz_avg)
        cal_factor = mol_atb / reference
        self.atb[wv] = (self.rawdata[wv_keys[wv]].T * cal_factor).T
        
        return cal_factor
        
    def vertical_downsample(self, new_z_step):
        '''
        vertically downsample the backscatter profiles. E.g go from the initial resolution of 30m to 100m.
        Does not upsample (e.g. resolution of 10m)
        Correctly takes into account partial boundaries between origin and target binsizes.
        
        new_z_step in km
        
        backscatter needs to be calibrated (using calibrate_atb_to) before being downsampled
        '''
        
        assert new_z_step > (self.altitude[1] - self.altitude[0]), 'new_z_step cannot be smaller than 30m'
        z2 = np.r_[-1.:40:new_z_step]
        for wv in self.atb:
            if self.atb[wv] is None:
                continue
            data2 = np.zeros([self.atb[wv].shape[0], z2.shape[0]])
            for i in np.r_[0:data2.shape[0]]:
                data2[i,:] = _downsample_profile(self.altitude, self.atb[wv][i,:], z2)
            self.atb[wv] = data2
        self.altitude = z2    
                
    def describe(self, prof=0):
    
        print 'Number of profiles in file: ', self.nprof
        print 'Profile : ', prof
        print '\tVersion number: ', self.rawdata['majorversionnumber'][prof], self.rawdata['minorversionnumber'][prof]
        print '\tOrbit number: ', self.rawdata['orbitnumber'][prof]
        print '\tID number: ', self.rawdata['idnumber'][prof]
        print '\tRaw Date: ', self.rawdata['gmtday'][prof], self.rawdata['gmthour'][prof], self.rawdata['gmtmin'][prof], self.rawdata['gmtsec'][prof], self.rawdata['gmthund'][prof]
        print '\tDatetime: ', self.datetimes[prof]
        print '\tMetDate: ', self.rawdata['metday'][prof], self.rawdata['methour'][prof], self.rawdata['metmin'][prof], self.rawdata['metsec'][prof], self.rawdata['methund'][prof]
        print '\tlat, lon: ', self.rawdata['latitude'][prof], self.rawdata['longitude'][prof]
        
    def plot_photon_profiles(self, wv=355):
        
        import matplotlib.pyplot as plt
        import matplotlib.dates as mdates
        
        ndates = mdates.date2num(self.datetimes)
        fig = plt.figure(figsize=[15,5])
        ax = plt.gca()
        plt.pcolormesh(ndates, self.altitude, self.rawdata[wv_keys[wv]].T)
        ax.xaxis.axis_date()
        plt.ylabel('Altitude')
        fig.autofmt_xdate()
        
        cb = plt.colorbar()
        cb.set_label('Photon counts')
        plt.show()


def main(f='LITE_L1_19940910_164558_164706', iprof=0):

    iprof = int(iprof)
    
    l = LITE(f)
    l.describe(iprof)
    l.plot_profiles(355)


if __name__ == '__main__':
    import plac
    plac.call(main)


# utility functions


def _downsample_profile(z, prof, z2):
    prof2 = np.zeros_like(z2)
    for i in np.r_[0:z2.shape[0]-1]:
        idx = np.where((z >= z2[i]) & (z < z2[i+1]))[0]
        weights = np.ones_like(idx)
        # if the last point is not fully inside the new bin, weight it down
        if z[idx[-1]] < z2[i+1]:
            # print z[idx[-1]], ' < ',  z2[i+1], ' -> weighting down the last'
            weights[-1] = (z2[i+1] - z[idx[-1]]) / (z[idx[-1]] - z[idx[-2]])
        if z[idx[0]] > z2[i]:
            # print z[idx[0]], ' > ' , z2[i], ' -> adding a bit of below'
            idx = np.append(idx, idx[0]-1)
            weights = np.append(weights, (z[idx[0]] - z2[i]) / z[idx[0]-1])
        prof2[i] = np.average(prof[idx], weights=weights)
    return prof2


def _rolling_average_padded(vector, window):
    import pandas as pd
    d = pd.Series(vector)
    d = pd.rolling_mean(d, window)
    vector = d.values
    vector[0:window] = vector[window]
    vector[-window:] = vector[-window]
    return vector
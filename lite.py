#!/usr/bin/env python
#encoding:utf-8

# Created by VNoel on 2014-08-11

'''
Class to read LITE lidar data.

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
    

class LITE(object):
    
    def __init__(self, filename):
        
        self.rawdata = np.fromfile(filename, dtype=header, count=-1)
        self.altitude = np.linspace(40, -4.985, 3000)
        self.nprof = self.rawdata['profile355'].shape[0]
        self.photons = dict()
        self.latitude = self.rawdata['latitude']
        self.longitude = self.rawdata['longitude']
        self.photons[355] = self.rawdata['profile355']
        self.photons[532] = self.rawdata['profile532']
        self.photons[1064] = self.rawdata['profile064']
        
        self.datetimes = []
        for i in np.r_[0:self.nprof]:
            date = datetime(1994,1,1) + timedelta(days=self.rawdata['gmtday'][i]-1)
            fulldate = datetime(date.year, date.month, date.day, self.rawdata['gmthour'][i], self.rawdata['gmtmin'][i], self.rawdata['gmtsec'][i], self.rawdata['gmthund'][i]*10000)
            self.datetimes.append(fulldate)
                
    def describe(self, prof=0):
    
        print 'Profile : ', prof
        print 'Version number: ', self.rawdata['majorversionnumber'][prof], self.rawdata['minorversionnumber'][prof]
        print 'Orbit number: ', self.rawdata['orbitnumber'][prof]
        print 'ID number: ', self.rawdata['idnumber'][prof]
        print 'Raw Date: ', self.rawdata['gmtday'][prof], self.rawdata['gmthour'][prof], self.rawdata['gmtmin'][prof], self.rawdata['gmtsec'][prof], self.rawdata['gmthund'][prof]
        print 'Datetime: ', self.datetimes[prof]
        print 'MetDate: ', self.rawdata['metday'][prof], self.rawdata['methour'][prof], self.rawdata['metmin'][prof], self.rawdata['metsec'][prof], self.rawdata['methund'][prof]
        print 'lat, lon: ', self.rawdata['latitude'][prof], self.rawdata['longitude'][prof]
        print 'Number of profiles: ', self.nprof
        
    def plot_profiles(self, wv=355):
        
        import matplotlib.pyplot as plt
        import matplotlib.dates as mdates
        
        ndates = mdates.date2num(self.datetimes)
        fig = plt.figure(figsize=[15,5])
        ax = plt.gca()
        plt.pcolormesh(ndates, self.altitude, self.photons[wv].T)
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



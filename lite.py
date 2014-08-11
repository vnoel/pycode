#!/usr/bin/env python
#encoding:utf-8

# Created by VNoel on 2014-08-11

import numpy as np

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
        
        fid = open(filename, 'r')
        self.rawdata = np.fromfile(fid, dtype=header, count=-1)
        fid.close()
        
        self.altitude = np.linspace(40, -4.985, 3000)
        self.nprof = self.rawdata['profile355'].shape[0]
                

def main(f = 'LITE_L1_19940910_164558_164706'):

    l = LITE(f)
    print 'First profile : '
    print 'Version number: ', l.rawdata['majorversionnumber'][0], l.rawdata['minorversionnumber'][0]
    print 'Orbit number: ', l.rawdata['orbitnumber'][0]
    print 'ID number: ', l.rawdata['idnumber'][0]
    print 'Date: ', l.rawdata['gmtday'][0], l.rawdata['gmthour'][0], l.rawdata['gmtmin'][0], l.rawdata['gmtsec'][0], l.rawdata['gmthund'][0]
    print 'MetDate: ', l.rawdata['metday'][0], l.rawdata['methour'][0], l.rawdata['metmin'][0], l.rawdata['metsec'][0], l.rawdata['methund'][0]
    print 'lat, lon: ', l.rawdata['latitude'][0], l.rawdata['longitude'][0]
    print 'Number of profiles: ', l.rawdata['profile355'].shape
    # import matplotlib.pyplot as plt
    # plt.pcolormesh(l.rawdata['latitude'], l.altitude, l.rawdata['profile355'].T)
    # plt.show()


if __name__ == '__main__':
    import plac
    plac.call(main)



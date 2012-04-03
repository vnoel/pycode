#!/usr/bin/env python

'''
Module to read Radiosonde measurement files in the PCCORA proprietary binary format.
See http://badc.nerc.ac.uk/data/ukmo-rad-hires/pc-coradata.html

Invalid values are masked in the output

This module can be called to print information out of a PCCORA file.
VNoel 2012 

dependencies: numpy, plac to run the module standalone (easy to remove)
'''

import numpy as np
from datetime import datetime, timedelta

i = 'i2'
pccoraheader = np.dtype( [ ('copyright', 'a20'), ('lenident', i), ('lensyspar', i), 
                            ('nrecdata', i), ('numsdtlevels', i), ('datatype', i), ('lendata', i),
                            ('readyflag', 'b1'),
                            ('resbytarr', '17b1')
                            ] )
pccoraident = np.dtype( [ ('stationtype', i), ('region', i), ('block', i), ('station', i),
                            ('lat', i), ('long', i), ('alt', i), ('windunit', i), ('teleheading', i),
                            ('reserved', i), ('soundingtype', i), ('startmode', i),
                            ('ascenttime', i), ('pturate', i),
                            ('cardnum', 'i4'),
                            ('year', i), ('month', i), ('day', i), ('jday', i), ('hour', i), ('min', i),
                            ('messyear', i), ('messmonth', i), ('messday', i), ('messhour', i),
                            ('cloudgroup', 'a6'), ('weathergroup', 'a6'), ('napp', 'a6'),
                            ('surfpress', i), ('surftemp', i), ('surfhum', i), 
                            ('surfwinddir', i), ('surfwindspeed', i),
                            ('sondenum', 'a10'), ('soundingnum', 'a10'),
                            ('presscorr', i), ('tempcorr', i), ('humcorr', i), ('succsig', i),
                            ('pressacc', i), ('pressrep', i), ('pressrej', i),
                            ('tempacc', i), ('temprep', i), ('temprej', i),
                            ('humacc', i), ('humrep', i), ('humrej', i),
                            ('totomega', i), ('termination', i),
                            ('omegacount', '22b1'),
                            ('windcompmode', i), ('windmode', i), ('omegastation', i),
                            ('loraninfo', '8b1'),
                            ('phaseint', '6i2'),
                            ('phaseintlev', '6i2'),
                            ('refpress', i), ('reftemp', i), ('refhum', i) ] )
                            
syspar =  np.dtype( [ ('junk', '8087b1') ] )

pccoradata = np.dtype( [ ('time', 'f4'), ('logpress', i), ('temp', i), ('hum', i), ('u', i), ('v', i),
                        ('alt', i), ('press', i), ('dew', i), ('mix', i), ('winddir', i), ('windspeed', i),
                        ('az', i), ('distance', i), ('lat', i), ('long', i), 
                        ('sigkey', '2b1'), ('editedsigkey', '2b1'), 
                        ('radarheight', i) ] )

rsfilepath = '/users/vnoel/Projects/rs/'
             
def find_pccora_file(date):
    import glob
    basepath = rsfilepath + '/%04d/%02d/' % (date.year, date.month)
    files = glob.glob(basepath + 'radiosonde_rothera_%04d%02d%02d*.edt' % (date.year, date.month, date.day))
    if not files:
        return None
    if len(files) > 1:
        print 'More than 1 file for that date, using first'
    files = files[0]
    return files
             
                      
def convert_object_to_dict(data, names):

    # convert the data in a regular dictionary
    # to permit type conversion
    data2 = {}

    for name in names:
        if data[name].dtype == 'i2':
            data2[name] = data[name].astype(np.float)
        else:
            data2[name] = data[name]

    return data2
  

def mask_missing_values(data):
    # check for missing values, mask them
    for dataname in data:
        idx = (data[dataname] == -32768.0)
        if np.sum(idx) == 0:
            continue
        data[dataname] = np.ma.masked_where(idx, data[dataname])
    # also mask points with no time data if relevant
    if 'time' in data:
        for dataname in data:
            if data[dataname].dtype == np.float:
                data[dataname] = np.ma.masked_where(data['time']==0, data[dataname])
        data['time'] = np.ma.masked_where(data['time']==0, data['time'])
    return data


def convert_units(data):

    data['alt'] += 30000.
    data['temp'] = data['temp'].astype(np.float) * 0.1  # Kelvin
    data['spress'] = np.exp(data['logpress']/4096.)
    data['press'] *= 0.1
    data['long'] *= 0.01
    data['lat'] *= 0.01
    data['windspeed'] *= 0.1
  
    return data
    

def create_data_datetimes(launchtime, seconds):

    datetimes = [launchtime+timedelta(seconds=int(x)) if (x is not np.ma.masked) else 0 for x in seconds]
    return datetimes
  
  
def pccora_read(file):

    ## HEADER
    
    fid = open(file, 'r')
    
    head = np.fromfile(fid, dtype=pccoraheader, count=1)
    
    ## DATA IDENTIFIER
    
    ident = np.fromfile(fid, dtype=pccoraident, count=1)
    ident = convert_object_to_dict(ident, pccoraident.names)
    ident = mask_missing_values(ident)
    ident['long'] *= 0.01
    ident['lat' ] *= 0.01
    ident['surfpress'] *= 0.1
    ident['datetime'] = datetime(ident['year'], ident['month'], ident['day'], ident['hour'], ident['min'])
    ident['launchtime'] = ident['datetime'] + timedelta(seconds=ident['ascenttime'])
    
    junk = np.fromfile(fid, dtype=syspar, count=1)
    
    data = np.fromfile(fid, dtype=pccoradata, count=25)
    data = convert_object_to_dict(data, pccoradata.names)
    data = mask_missing_values(data)
    data = convert_units(data)
    data['datetime'] = create_data_datetimes(ident['launchtime'], data['time'])

    hires = np.fromfile(fid, dtype=pccoradata, count=head['nrecdata'])
    hires = convert_object_to_dict(hires, pccoradata.names)
    hires = mask_missing_values(hires)
    hires = convert_units(hires)
    hires['datetime'] = create_data_datetimes(ident['launchtime'], data['time'])
    
    fid.close()
    
    return head, ident, data, hires
    
    
def main(file='/users/noel/Projects/blue5/rs/2008/07/radiosonde_rothera_2008070111.edt'):
    
    head, ident, data, hires = pccora_read(file)
    print '## HEADER'
    ndata = head['nrecdata']
    
    print 'Copyright', head['copyright']    
    print 'len ident', head['lenident']
    print 'datatype:',head['datatype']
    print 'ndata:', ndata
    print

    print '## DATA IDENTIFIER'
    
    print 'Station type: ',ident['stationtype']
    print 'Coords: ', ident['lat'], ident['long']
    print 'Date: ', ident['year'], ident['month'], ident['day'], ident['hour'], ident['min']
    print 'Message data: ', ident['messyear'], ident['messmonth'], ident['messday'], ident['messhour']
    print 'Cloud group: ', ident['cloudgroup']
    print 'Weather group: ', ident['weathergroup']
    print 'Napp: ', ident['napp']
    print 'Surface pressure [hPa]: ', ident['surfpress'] 
    print 'Surface temperature [K]: ', 0.1 * ident['surftemp']
    print 'Surface humidity [% RH]: ', ident['surfhum']
    print 'Surface winddir [deg]:', ident['surfwinddir']
    print 'Surface windspeed [m/s]:', ident['surfwindspeed']*0.1
    print 'Temperature correction: ', ident['tempcorr']
    print 'Sonde number: ', ident['sondenum']
    print 'Sounding number: ', ident['soundingnum']
    print 'Wind computing mode: ', ident['windcompmode']
    print 'Wind mode: ', ident['windmode']
    print 'Reference pressure, temp, hum:', ident['refpress'], ident['reftemp'], ident['refhum']
    
    print '## DATA'
    
    for i in np.r_[0:25]:
        print i
        print 'Elapsed time since release [s]:', data['time'][i]
        print 'Time: ', data['datetime'][i]
        print 'Scaled P and P [hPa]:', data['spress'][i], data['press'][i]
        print 'Temperature [K, C]: ', data['temp'][i], data['temp'][i]-273
        print 'Wind direction[deg]: ', data['winddir'][i]
        print 'Wind speed[m/s]: ', data['windspeed'][i]
        print 'Coords: ', data['long'][i], data['lat'][i]


    
if __name__ == '__main__':
    import plac
    plac.call(main)

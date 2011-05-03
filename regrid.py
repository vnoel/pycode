#!/usr/bin/env python
#encoding: utf-8

'''
    Regridding utilities
    Noel V. - LMD/CNRS
'''

import numpy as np
from scipy.interpolate import interp1d

def regrid_profiles(pvec, pprof, vprof):
    '''
    regrid_profiles(pvec, pprof, vprof)
        regrids profiles of variable vprof on the pressure vector pvec
        pprof, vprof [nprof, n1] - profiles of pressure and variable.  Pressure must be decreasing.
        pvec [n2]
        output : vprof2 [nprof, n2]
    '''

    assert pvec.ndim == 1, 'pvec must be a vector'
    assert (pprof.ndim==2) & (vprof.ndim==2), 'Wrong dimensions in pprof or vprof'
    assert np.abs(np.max(pvec)-np.max(pprof)) < 800, 'Check your pressures (possible linear/log confusion)'
    
    nprof = vprof.shape[0]    
    vprof2 = np.empty([nprof, pvec.size])
    for i in np.arange(nprof):
        if not np.all(np.diff(pprof[i,:]) < 0):
            print 'Warning : Pressure values are not all decreasing'
        vprof2[i,:] = np.interp(pvec, pprof[i,::-1], vprof[i,::-1])
        
    return vprof2
    
def regrid_array(pvec, parr, varr, kind='linear'):
    '''
    regrid_array(pvec, parr, varr)
        regrids an array of variable varr on a pressure vector pvec
        parr, varr [nz, nx, ny] - profiles of pressure and variable. Pressure must be decreasing.
        pvec [nz2]
        output : varr2 [nz2, nx, ny]
    '''
    
    assert pvec.ndim == 1, 'pvec must be a vector'
    assert (parr.ndim==3) & (varr.ndim==3), 'parr and varr must have 3 dimensions'
    assert np.abs(np.max(pvec)-np.max(parr)) < 800, 'Check your pressures (possible linear/log confusion)'
    
    if parr.shape[0] == (varr.shape[0]-1):
        print 'Warning: varr and parr not same vertical length. Trying to fix it...'
        varr = varr[:-1,:,:]
    
    assert parr.shape[0] == varr.shape[0], 'Pressure and variable must have same vertical length'
    
    n1, n2, n3 = varr.shape
    varr2 = np.empty([pvec.size, n2, n3])
    for i in np.arange(n2):
        for j in np.arange(n3):
            if not np.all(np.diff(parr[:,i,j]) < 0):
                print 'Warning : Pressure values are not all decreasing'
            if kind is 'linear':
                varr2[:,i,j] = np.interp(pvec, parr[::-1,i,j], varr[::-1,i,j])
            elif kind in ('nearest', 'cubic', 'zero', 'slinear', 'quadratic'):
                f = interp1d(parr[::-1,i,j], varr[::-1,i,j], kind=kind, bounds_error=False)
                varr2[:,i,j] = f(pvec)
            else:
                print 'Unknown interpolation method'
                
    return varr2

def _pvec_log(n=1000):
    npres = n
    pmaxlog = np.log10(1000)
    pminlog = np.log10(1)
    psteplog = (pmaxlog - pminlog) / npres
    pveclog = np.r_[pmaxlog:pminlog:-psteplog]
    return pveclog

def regrid_profiles_on_logp(pprof, vprof):
    '''
    Creates a linear log-step pressure vector and call regrid_profiles on it.
        output: pvec, vprof2
    '''
    pvec = _pvec_log(n=100)
    vprof2 = regrid_profiles(pvec, pprof, vprof)
    return pvec, vprof2
    
def regrid_array_on_logp(parr, varr, n=1000, kind='linear'):
    '''
    Creates a linear log-step pressure vector and call regrid_array on it.
        output: pvec, varr2
    '''
    pvec = _pvec_log(n=n)
    varr2 = regrid_array(pvec, parr, varr, kind=kind)
    return pvec, varr2
    
wrftestfile = '/users/noel/Projects/blue_penguin/analysis2006/wrfout_d01_2006-06-27_00:00:00'

def _pcolor_w(ax, lat, pres, w, logscale=False):
    print lat.shape, pres.shape, w.shape
    pc = ax.pcolormesh(lat, pres, w, vmin=-2, vmax=2)
    ax.set_xlim(-75, -64)
    if logscale:
        ax.set_yscale('log')
        ax.set_ylim(1000, 10)
    else:
        ax.set_ylim(3, 1)


def test_array_regridding():
    import wrf
    import matplotlib.pyplot as plt
    import numpy as np
    
    print 'Reading WRF 3.2.1 data'
    it = 1      # 03:00 on 2006-06-27
    wrffile = wrf.wrf(wrftestfile)
    wrftime = wrffile.time(it=it)
    wrfptop = wrffile.p_top(it=it)
    wrfw = wrffile.w(it=it)
    wrfpres = wrffile.pressure(it=it)
    wrffile.close()
    
    print wrftime, wrfptop
    wrfw = wrfw[:-1,:,:]
    
    print wrfw.shape, wrfpres.shape
    
    # 1st test : regridding on a linear-step pressure vector
    pvec = np.r_[1000:1:-1.]
    wrfw2 = regrid_array(pvec, wrfpres, wrfw)
    
    # show a random profile
    fig = plt.figure()
    plt.plot(wrfw[:,50,50], wrfpres[:,50,50], label='On WRF Pressure')
    plt.plot(wrfw2[:,50,50], pvec, label='regridding on linear pressure')
    plt.ylim(1000, 10)
    plt.legend()
    
    fakelat = np.r_[-75:-64:(75-64)/99.]
    
    fig = plt.figure()
    ax = plt.subplot(2,1,1)
    _pcolor_w(ax, fakelat, wrfpres[:,:,50], wrfw[:,:,50], logscale=True)
    ax.set_title('On WRF Pressure grid, logscale')
    
    ax = plt.subplot(2,1,2)
    _pcolor_w(ax, fakelat, pvec, wrfw2[:,:,50], logscale=True)
    ax.set_title('Regridded on constant linear-step pressure')
    
    # 2nd test: regridding on a constant log-step pressure vector
    pveclog = _pvec_log()
    wrfw2log = regrid_array(pveclog, np.log10(wrfpres), wrfw)
    
    # show a random profile
    fig = plt.figure()
    plt.plot(wrfw[:,50,50], np.log10(wrfpres[:,50,50]), label='On WRF pressure')
    plt.plot(wrfw2log[:,50,50], pveclog, label='regridded on linear log-step pressure')
    plt.ylim(3,1)
    plt.legend()
    
    fig = plt.figure()
    ax = plt.subplot(2,1,1)
    _pcolor_w(ax, fakelat, np.log10(wrfpres[:,:,50]), wrfw[:,:,50])
    ax.set_title('On WRF Pressure grid')
    ax.set_ylabel('Log10(P)')
    ax = plt.subplot(2,1,2)
    _pcolor_w(ax, fakelat, pveclog, wrfw2log[:,:,50])
    ax.set_title('regridded on constant linear log-step pressure')
    ax.set_ylabel('log10(P)')
    
    plt.show()
    

def test_profiles_regridding():
    import calipso
    import wrf
    import matplotlib.pyplot as plt
    
    n1=1300
    n2=1625
    
    print 'Reading CALIOP data'
    calfile = calipso.Cal1('/bdd/CALIPSO/Lidar_L1/CAL_LID_L1.v2.01/2006/2006_06_27/CAL_LID_L1-Prov-V2-01.2006-06-27T03-28-22ZN.hdf')
    lon, lat = calfile.coords(navg=30)
    lon, lat = lon[n1:n2], lat[n1:n2]
    calfile.close()

    print 'Reading WRF 3.2.1 data'
    it = 1      # 03:00 on 2006-06-27
    wrffile = wrf.wrf(wrftestfile)
    wrftime = wrffile.time(it=it)
    wrfptop = wrffile.p_top(it=it)
    wrfw = wrffile.w(it=it, on_orbit=(lon, lat))
    wrfpres = wrffile.pressure(it=it, on_orbit=(lon, lat))
    wrffile.close()

    print wrftime, wrfptop
    wrfw = wrfw[:,:-1]
    
    # 1st try : regridding on a linear-step pressure vector
    pvec = np.r_[1000:1:-1.]
    print 'Interpolating stuff'
    wrfw2 = regrid_profiles(pvec, wrfpres, wrfw)

    # showing a random profile
    fig = plt.figure()
    plt.plot(wrfpres[100,:], wrfw[100,:], label='on WRF pressure')
    plt.plot(pvec, wrfw2[100,:], label='regridded on linear pressure')
    plt.legend()

    # showing arrays
    # the first plot - pressure-based grid (nonlinear)
    fig = plt.figure()
    ax = plt.subplot(2,1,1)
    _pcolor_w(ax, lat, wrfpres.T, wrfw.T, logscale=True)
    ax.set_title('On WRF pressure grid, logscale')
    
    # second plot - linear pressure grid
    ax = plt.subplot(2,1,2)
    _pcolor_w(ax, lat, pvec, wrfw2.T, logscale=True)
    ax.set_title('Regridded on constant linear-step pressure')

    # 2nd try : regridding on a log-step pressure vector
    pveclog = pvec_log()
    wrfw2log = regrid_profiles(pveclog, np.log10(wrfpres), wrfw)

    # showing a random profile
    fig = plt.figure()
    plt.plot(np.log10(wrfpres[100,:]), wrfw[100,:], label='On WRF pressure')
    plt.plot(pveclog, wrfw2log[100,:], label='Regridded on log pressure')
    plt.legend()

    # showing arrays
    fig = plt.figure()
    ax = plt.subplot(2,1,1)
    _pcolor_w(ax, lat, np.log10(wrfpres.T), wrfw.T)
    ax.set_title('On WRF log10(pressure) grid')
    
    ax = plt.subplot(2,1,2)
    _pcolor_w(ax, lat, pveclog, wrfw2log.T)
    ax.set_title('Regridded on constant log-step pressure')

    plt.show()
    

if __name__ == '__main__':
    #test_profiles_regridding()
    test_array_regridding()
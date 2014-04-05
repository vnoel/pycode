# thermo.py

# useful functions for water vapor stuff...

from scipy.optimize import bisect, brentq
import numpy as np

def hno3_partial_pressure(mixing_ratio, pressure):
    '''
    Arguments:
        mixing ratio in ppmv
        pressure in Pa
    '''
    
    return (mixing_ratio*1e-6) * pressure
    
def hno3_saturation_pressure_ln(T, wvp):
    '''
    estimates the hno3 saturation pressure as a function of temperature and water vapor partial pressure
        T - temperature in K
        wvp - water vapor partial pressure in Pa
        Hanson and Mauersberger 1988
    '''

    if T > 200:
        log10_hno3sat = 13.622 - 3561.3 / T
        # print 'log10_hno3sat : ', log10_hno3sat
    elif T <= 200:    
        # pressure in Torr
        wvp_torr = wvp / 133.322
        mT = -2.7836 - 0.00088*T
        bT = 38.9855 - (11397. / T) + (0.009179*T)
        log10_hno3sat = mT * np.log10(wvp_torr) + bT
    
    # revert in Pa...

    hno3sat = np.power(10., log10_hno3sat)
    # print 'hno3sat : ', hno3sat, ' Torr'
    hno3sat *= 133.322
    # print 'hno3sat : ', hno3sat, ' Pa'
    ln_hno3sat = np.log(hno3sat)
    
    return ln_hno3sat

def hno3_frost_point(hno3mr, wvmr, pressure):
    '''
    Arguments:
        hno3mr = hno3 mixing ratio (ppmv)
        wvmr = water vapor mixing ratio (ppmv)
        pressure = pressure in Pa
    '''
    
    if hno3mr.ndim > 0:
        
        frost_point = np.zeros_like(hno3mr)
        for i in np.r_[0:hno3mr.shape[0]]:            

            hno3p = hno3_partial_pressure(hno3mr[i], pressure[i])
            wvp = water_vapor_partial_pressure_over_ice(wvmr[i], pressure[i])
            ln_hno3p = np.log(hno3p)

            f = lambda T: (ln_hno3p - hno3_saturation_pressure_ln(T, wvp))

            frost_point[i] = brentq(f, 170., 210.)
            if frost_point[i] < 180:
                print frost_point[i]
                print hno3mr[i]
                print wvmr[i]
                print pressure[i]
            
    else:
        
        hno3p = hno3_partial_pressure(hno3mr, pressure)
        wvp = water_vapor_partial_pressure_over_ice(wvmr, pressure)
        ln_hno3p = np.log(hno3p)

        f = lambda T: (ln_hno3p - hno3_saturation_pressure_ln(T, wvp))

        frost_point = brentq(f, 170., 210.)
        

    return frost_point
    # if np.size(pressure) is 1:
    #     # find the temperature for which ln_e = ln_e_sat
    #     frost_point = bisect(f, 170., 233., args=(ln_e))
    # else:
    #     frost_point = np.zeros_like(ln_e)
    #     ln_e_flat = ln_e.flatten()
    #     for i in np.arange(ln_e_flat.size):
    #         indices = np.unravel_index(i, ln_e.shape)
    #         frost_point[indices] = brentq(f, 170., 233., args=(ln_e[indices]))
            
    

def water_vapor_partial_pressure_over_ice(mixing_ratio, pressure):
    '''
    Arguments:
        mixing ratio in ppmv
        pressure in Pa
    '''
    return (mixing_ratio*1e-6) * pressure

def water_vapor_saturation_pressure_over_ice_ln(T):
    '''
        estimates the water vapor saturation pressure over ice as a function of temperature
        from Murphy and Koop 2005
        see http://cires.colorado.edu/~voemel/vp.html
        argument :
            temperature in Kelvins
        returns :
            water vapor saturation pressure over ice in Pa
    '''
    
    v = [9.550426, 5723.265, 3.53068, 0.00728332]
    ln_esat = v[0] - v[1]/T + v[2]*np.log(T) - v[3]*T
    return ln_esat

# compute look-up tables : ln_esat = f(ft)
ft = np.r_[170:233:0.1]
ln_esat_as_f_t = water_vapor_saturation_pressure_over_ice_ln(ft)
ln_esat_min = np.min(ln_esat_as_f_t)
ln_esat_max = np.max(ln_esat_as_f_t)
ln_esat_range = ln_esat_max - ln_esat_min
# look-up table : x = ln_esat_lut (esat, fixed step) -> fn_as_f_ln_esat = f(x) (frost point, interpolated)
n_esat = 1000.
ln_esat_lut = np.linspace(ln_esat_min, ln_esat_max, num=n_esat)
ft_as_f_ln_esat = np.interp(ln_esat_lut, ln_esat_as_f_t, ft)

def water_vapor_saturation_pressure_over_ice_ln_lut(T):
    i = np.argmin(np.abs(T-ft_lut))
    ln_esat = ln_esat_lut[i]
    return ln_esat

def frost_point_temperature(mixing_ratio, pressure):
    '''
        Find the frost point temperature using lookup tables of Tsat = f(esat)
        
        arguments :
            mixing ratio (in ppmv) and a pressure (in Pa)
            arguments must be either both scalar or both arrays
        returns :
            frost point temperature in K
            
        for arrays of size [100,100,100], 90.7 s per loop. Still slow.
    '''
    
    ln_e = np.log(water_vapor_partial_pressure_over_ice(mixing_ratio, pressure))
    
    # simple, using lookup tables
    
    if np.size(pressure) is 1:
        # find the temperature for which ln_e = ln_e_sat
        ii = np.floor((ln_e - ln_esat_min) / ln_esat_range * n_esat)
        frost_point = ft_as_f_ln_esat[ii]
    else:
        frost_point = np.zeros_like(ln_e)
        ln_e_flat = ln_e.flatten()
        ii = np.floor((ln_e - ln_esat_min) / ln_esat_range * n_esat)
        for i in np.arange(ln_e_flat.size):
            indices = np.unravel_index(i, ln_e.shape)
            frost_point[indices] = ft_as_f_ln_esat[ii[indices]]
            
    return frost_point
    
def frost_point_temperature_slow(mixing_ratio, pressure):
    '''
    Comme frost_point_temperature... sans look-up tables.
    plus correct mais plus lent.
    
    '''
    
    ln_e = np.log(water_vapor_partial_pressure_over_ice(mixing_ratio, pressure))
    f = lambda x, ln_e: (ln_e - water_vapor_saturation_pressure_over_ice_ln(x))
     
    if np.size(pressure) is 1:
        # find the temperature for which ln_e = ln_e_sat
        frost_point = bisect(f, 170., 233., args=(ln_e))
    else:
        frost_point = np.zeros_like(ln_e)
        ln_e_flat = ln_e.flatten()
        for i in np.arange(ln_e_flat.size):
            indices = np.unravel_index(i, ln_e.shape)
            frost_point[indices] = brentq(f, 170., 233., args=(ln_e[indices]))
            
    return frost_point
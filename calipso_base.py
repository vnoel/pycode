#!/usr/bin/env python
#encoding:utf-8

# Generic calipso file class
# Use Cal1 and Cal2 classes instead

# need pyhdf for hdf4
from pyhdf.SD import SD, SDC
import warnings
import datetime
import numpy as np


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
        self.id = filename[-25:-4]
        self.date = datetime.datetime(self.year, self.month, self.day,
                                      self.hour, self.minutes, self.seconds)

    def __repr__(self):
        return self.filename

    def close(self):
        self.hdf.end()
        self.hdf = None


# Useful maths

def _vector_average(v0, navg, missing=None, valid=None):
    """
    v = _vector_average (v0, navg)
    moyenne le vector v0 tous les navg points.
    """

    v0 = v0.squeeze()

    assert v0.ndim == 1, 'in _vector_average, v0 should be a vector'
    if navg == 1:
        return v0

    n = np.floor(1. * np.size(v0, 0) / navg)
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
        print("_array_average: navg is even, turning weights off")

    # create triangle-shaped weights
    if weighted:
        w = np.zeros(navg)
        w[:navg / 2. + 1] = np.r_[1:navg / 2. + 1]
        w[navg / 2. + 1:] = np.r_[int(navg / 2.):0:-1]
    else:
        w = None

    # create averaged array a with number of averaged profiles n
    n = np.floor(1. * np.size(a0, 0) / navg)
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
            interpole les donnees du tableau z0 sur un nouveau y.
            utile pour regridder les donnees meteo genre temp
    """

    z = np.zeros([np.size(z0, 0), np.size(y)])
    for i in np.arange(np.size(z, 0)):
        z[i, :] = np.interp(y, np.flipud(y0), np.flipud(z0[i, :]))
    return z

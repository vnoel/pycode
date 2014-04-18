#!/usr/bin/env python
# encoding: utf-8

"""
arraydict.py

My own little module to read and write numpy arrays.

Created by Vincent Noel - LMD/CNRS on 2011-11-24.
"""

import numpy as np


class ArrayDict(dict):
    """
    ArrayDict objects are dictionaries containing named numpy arrays
    """

    def __init__(self, from_file=None, **kwargs):
        """
        an ArrayDict can be created empty
            x = ArrayDict()
        or filled with content from a npz file
            x = ArrayDict('file.npz')
        or merge content from several npz files
        (if there's a * or ? in the filename)
            x = ArrayDict('dir/*.npz')
        or from variables
            x = ArrayDict(lon=lon, lat=lat)
        """
        dict.__init__(self)

        if from_file:

            if ('?' in from_file) or ('*' in from_file):
                # from_file looks like a file pattern

                import glob

                filelist = glob.glob(from_file)
                if not filelist:
                    print('no file matching pattern', from_file)
                else:
                    filelist.sort()
                    print('Aggregating data from %d files' % len(filelist))

                    for f in filelist:
                        try:
                            this_array = ArrayDict(from_file=f)
                        except:
                            continue
                        self.append(this_array)

            else:

                npz = np.load(from_file)
                fill_value = None
                if 'fill_value' in npz:
                    fill_value = npz['fill_value']
                for f in npz.files:
                    if f is 'fill_value':
                        continue
                    self[f] = npz[f]
                    if fill_value is not None:
                        idx = (self[f] == fill_value)
                        if idx.sum() > 0:
                            self[f] = np.ma.masked_where(self[f] == fill_value, self[f])
                npz.close()

        if kwargs:
            for key in kwargs:
                self[key] = kwargs[key]

    def append(self, ad, axis=0):
        """
        Appends numpy arrays contained in another arraydict with those present in self.
        0-d arrays are ignored.
        """

        arrnames = list(ad.keys())
        if not arrnames:
            return

        for arrname in arrnames:
            if np.shape(ad[arrname]) is ():
                continue
            if arrname in list(self.keys()):
                self[arrname] = np.concatenate((self[arrname], ad[arrname]), axis=axis)
            else:
                self[arrname] = ad[arrname]


    def list(self):
        """
        display the list of arrays contained in self and their shape
        """

        for arrname in list(self.keys()):
            print(arrname + ':', self[arrname].shape)


    def save(self, filename, verbose=True, fill_value=-99999.):
        """
        save the arrays in a numpy file
        :type fill_value: float
        """
        if verbose:
            print('Saving', filename)

        # special-case masked arrays
        masked = False
        for key in self:
            if np.ma.isMaskedArray(self[key]):
                masked = True
        if masked:
            print('Saving masked arrays with fill_value=%f' % fill_value)
            copy = arraydict(**self)
            for key in copy:
                if np.ma.isMaskedArray(copy[key]):
                    copy[key] = np.ma.filled(self[key], fill_value=-99999.)
            np.savez_compressed(filename, fill_value=fill_value, **copy)
        else:
            np.savez_compressed(filename, **self)


    def dump(self, filename):
        """
        save the arrays in a numpy file
        """
        np.savez(filename, **self)


    def get_vars(self, varnamelist):
        """
        input: a list containing names of variables
        returns: a list containing the associated variables
        """
        varlist = []
        for varname in varnamelist:
            varlist.append(self[varname])
        return varlist


    def subset(self, idx):
        """
        filters out the contained variables along their first dimension according to an index vector
        the index vector must have the same number of items in the first dimension as every variable in the arraydict
        """

        for arrname in list(self.keys()):
            self[arrname] = self[arrname][idx, ...]


class arraydict(ArrayDict):
    pass


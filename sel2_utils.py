#!/usr/bin/env python
# encoding: utf-8
"""
sel2_utils.py
Created by Vincent Noel on 2010-12-02, LMD/CNRS.
"""

import numpy as np

def create_cloud_id(base, top):
    
    '''
    Create cloud IDs and attribute them to layers.
    Cloud IDs begin at 0.
    '''
    
    nprof = base.shape[0]
    nl = base.shape[1]

    cloud_id = np.zeros_like(base)

    def cloud_id_for_layer(iprof, ilayer):

        b, t = base[iprof,ilayer], top[iprof,ilayer]
        if b < 0 or t < 0:
            return -9999.

        if cloud_id[iprof,ilayer] > 0:
            # if the layer already has an ID, return it
            return cloud_id[iprof,ilayer]

        if (iprof+1) >= nprof:
            # if the layer is the last profile, new ID
            return np.max(cloud_id)+1
        
        for j in np.arange(nl):
            # otherwise look right to find an overlapping layer with an ID
            if (base[iprof+1,j] > t) or (top[iprof+1,j] < b):
                continue
                
            # if we get here the layer [iprof+1,j] overlaps
            if cloud_id[iprof,ilayer] == 0:
                # current layer has no ID yet
                # take the ID of the overlapping layer
                cloud_id[iprof,ilayer] = cloud_id_for_layer(iprof+1,j)
            else:
                # current layer has an ID
                # push the ID to the overlapping layer if not ID yet
                if cloud_id[iprof+1,j] == 0:
                    push_cloud_id_to_layers(iprof+1,j,cloud_id[iprof,ilayer])
                # else:
                #     if cloud_id[iprof+1,j] != cloud_id[iprof,ilayer]:
                        # print 'Problem: current ID %d, overlapping layer with ID %d' % (cloud_id[iprof,ilayer], cloud_id[iprof+1,j])
                
                    
        if cloud_id[iprof,ilayer] == 0:
            # no overlapping layers, create new ID
            cloud_id[iprof,ilayer] = np.max(cloud_id) + 1
                        
        return cloud_id[iprof,ilayer]
        

    def push_cloud_id_to_layers(iprof,ilayer,cid):
        
        cloud_id[iprof,ilayer] = cid
        if iprof+1 >= nprof:
            return
            
        for j in np.arange(nl):
            if (base[iprof+1,j] > top[iprof,ilayer]) or (top[iprof+1,j] < base[iprof,ilayer]):
                continue
                
            if cloud_id[iprof+1,j] > 0:
                # if cloud_id[iprof+1,j] != cid:
                    # print 'Problem: right layer has a cloud ID %d and present %d' % (cloud_id[iprof+1,j], cid)
                continue
                
            push_cloud_id_to_layers(iprof+1,j,cid)   


    # loop on profiles and layers

    for i in np.arange(nprof):
        for j in np.arange(nl):            
            if cloud_id[i,j]==0:
                cloud_id[i,j] = cloud_id_for_layer(i, j)
                           
    # make IDs begin at 0 
    cloud_id = cloud_id - 1
    cloud_id[cloud_id<0] = -9999.
                            
    return cloud_id
    
    
def find_cloud_horizontal_extension(cloud_id, hstep=5):
    '''
    this function returns a vector of length n_id (number of individual cloud layers)
    containing the horizontal extent of the associated cloud layer (indexed through its ID)

    cloud_id is a layer array containing the cloud id for each layer
    
    hstep is the horizontal distance between two profiles, 5 km for averaging 15 profiles
    '''

    n_id = np.max(cloud_id) + 1
    hext = np.zeros(n_id)

    for cid in np.arange(n_id):
        # find the layers with with cloud id
        idx = (cloud_id == cid)
        # find profiles with this cloud id
        idx_prof = idx.sum(axis=1)
        # number of profiles where this id shows up = length in profiles of the cloud
        hext[cid] = np.sum(idx_prof > 0) * hstep
        
    return hext

def find_cloud_maximum_profile_extension(cloud_id, base, top):
    '''
    this function returns a vector of length n_id (number of individual cloud layers)
    containing the maximum vertical extent of the associated cloud layer (indexed through its ID)

    cloud_id is a layer array containing the cloud id for each layer    
    '''
    
    n_id = np.max(cloud_id) + 1
    vext = np.zeros(n_id)
    for cid in np.arange(n_id):
        idx = (cloud_id == cid)
        if np.sum(idx) == 0:
            continue
        cz = top[idx] - base[idx]
        vext[cid] = np.max(cz)
        
    return vext

def find_cloud_minimum_profile_extension(cloud_id, base, top):
    '''
    this function returns a vector of length n_id (number of individual cloud layers)
    containing the maximum vertical extent of the associated cloud layer (indexed through its ID)

    cloud_id is a layer array containing the cloud id for each layer    
    '''

    n_id = np.max(cloud_id) + 1
    vext = np.zeros(n_id)
    for cid in np.arange(n_id):
        idx = (cloud_id == cid)
        if np.sum(idx) == 0:
            continue
        cz = top[idx] - base[idx]
        vext[cid] = np.min(cz)

    return vext

def find_cloud_profile_range(cloud_id, base, top):
    '''
    this function returns two vectors of length n_id (number of individual cloud layers)
    containing the maximum and minimum vertical profile extents of the associated cloud layer (indexed through its ID)

    cloud_id is a layer array containing the cloud id for each layer    
    '''

    n_id = np.max(cloud_id) + 1
    vextmin = np.zeros(n_id)
    vextmax = np.zeros(n_id)
    for cid in np.arange(n_id):
        idx = (cloud_id == cid)
        cz = top[idx] - base[idx]
        vextmax[cid] = np.max(cz)
        vextmin[cid] = np.min(cz)

    return vextmin, vextmax
    
def find_cloud_maximum_vertical_range(cloud_id, base, top):
    '''
    this function returns two vectors of length n_id (number of individual cloud layers)
    containing the maximum and minimum vertical extents of the associated cloud layer (indexed through its ID)
    
    This function returns the difference between the highest top and the lowest base for all cloud layers

    cloud_id is a layer array containing the cloud id for each layer    
    '''

    n_id = np.max(cloud_id) + 1
    vextmax = np.zeros(n_id)
    for cid in np.arange(n_id):
        idx = (cloud_id == cid)
        if np.sum(idx) == 0:
            continue
        vextmax[cid] = np.max(top[idx]) - np.min(base[idx])

    return vextmax
    
def find_cloud_base_variation(cloud_id, base):
    '''
    this function returns two vectors of length n_id (number of individual cloud layers)
    containing the maximum and minimum vertical extents of the associated cloud layer (indexed through its ID)

    cloud_id is a layer array containing the cloud id for each layer    
    '''

    n_id = np.max(cloud_id) + 1
    bvar = np.zeros(n_id)
    for cid in np.arange(n_id):
        idx = (cloud_id == cid)
        bvar[cid] = np.max(base[idx]) - np.min(base[idx])

    return bvar

def layer_average_in_cloud(data, cloud_id):
    '''
    Average a given property (data) over layers with the same cloud ID
    '''
    
    n_id = np.max(cloud_id) + 1
    avg = np.zeros(n_id)
    
    for cid in np.arange(n_id):
        idx = (cloud_id == cid) & (np.isfinite(data))
        avg[cid] = np.mean(data[idx])
        
    return avg

def layer_max_in_cloud(data, cloud_id):
    '''
    Find the maximum for a given property (data) over layers with the same cloud ID
    '''

    n_id = np.max(cloud_id) + 1
    m = np.zeros(n_id)

    for cid in np.arange(n_id):
        idx = (cloud_id == cid) & (np.isfinite(data))
        m[cid] = np.max(data[idx])

    return m
    
    

def main():
    pass

if __name__ == '__main__':
    main()

#!/usr/bin/env python
#encoding: utf-8

'''
Useful function for studying PSCs
'''

ID_ALL     = 0b11111
ID_STS     = 0b00001
ID_NAT     = 0b00010
ID_ICE     = 0b00100
ID_ENHNAT  = 0b01000
ID_UNKNOWN = 0b10000

# ID_ENHNAT is the "Mix2-enh" from Pitts et al. 2011 (ACP)

psc_id = {'unknown':ID_UNKNOWN, 'all':ID_ALL, 'STS':ID_STS, 'NAT':ID_NAT, 'Ice':ID_ICE, 'ENHNAT':ID_ENHNAT}

def classify_psc_type(r, depol):
    '''
    give a type to a psc layer according to its r and depol coefficients.
    r = ratio between b532 and bmol, with b532 the total volume backscatter and bmol the molecular backscatter
    depol = particle-only depolarization ratio
    (cf Pitts et al 2009)
    '''
    
    sts_max_depol = 0.035
    sts_min_r = 1.25
    ice_min_r = 4.5 # 3.5

    ice_min_depol = sts_max_depol
    nat_max_r = ice_min_r
    
    if (r >= ice_min_r) and (depol >= ice_min_depol) and (depol < 0.8):
        psc_type = psc_id['Ice']
    elif (r >= sts_min_r) and (depol < sts_max_depol) and (depol >= -0.2):
        psc_type = psc_id['STS']
    elif (r > 2) and (r < 5) and (depol > 0.1):
        # mix2 enhanced NAT
        psc_type = psc_id['ENHNAT']
    elif (depol > sts_max_depol) and (depol < 0.8) and (r < nat_max_r) and (r > 1.0):
        psc_type = psc_id['NAT']
    elif (depol < sts_max_depol) and (depol >= 0.) and (r < sts_min_r) and (r > 1.0):
        psc_type = psc_id['NAT']
    else:
        psc_type = psc_id['unknown']
        
    return psc_type

def _plot_type_colors():
    colors = ['#000033', '#0000FF', '#00FF00', '#FF0000']
    return colors

def plot_type_cmap():
    from matplotlib.colors import ListedColormap
    colors = _plot_type_colors()
    cmap = ListedColormap(colors)
    return cmap
    
def plot_type_legend():
    from matplotlib.patches import Rectangle
    import matplotlib.pyplot as plt
    colors = _plot_type_colors()
    liste_rec = [Rectangle((0, 0), 1, 1, fc=color) for color in colors]
    labels = [types[i] for i in types]
    plt.legend(liste_rec, labels)
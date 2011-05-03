
from shapely.geometry import Polygon, Point

# definition du polygone dans lequel on a des CO
dp_l = [0.01, 0.17, 0.08, 0.01, 0.01]
atb_l = [0.05, 0.095, 0.2, 0.2, 0.05]
co_zone = Polygon(zip(dp_l, atb_l))

#definition du polygone dans lequel on a des nuages de glace
dp_i = [0.13, 0.3, 0.55, 0.55, 0.2]
atb_i = [-0.01, 0.07, 0.06, -0.01, -0.01]
izone = Polygon(zip(dp_i, atb_i))

def is_oriented(idp, iatb):
    p = Point(idp, iatb)
    return co_zone.contains(p)

def is_ice(idp, iatb):
    p = Point(idp, iatb)
    return izone.contains(p)

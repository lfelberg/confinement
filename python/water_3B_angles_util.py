import numpy as np
import itertools

from util import d_pbc, translate_pbc
from water_angles_util import angle_between, find_closest

RAD_TO_DEG = 180./np.pi

def cal_ang(o_coords, rng, dim = 3):
    '''Given a list of coords (dimensions: [nat = 1, nwat, ndim = 3]), 
       move the coords into the nearest image of water of interest,
       calculate the distances to 2 nearest neighbors and angle''' 
    d1, d2, t1 = [], [], []

    for i in range(o_coords.shape[0]):
        curr = o_coords[i]; others = np.delete(o_coords,i,0)
        ot_wr = translate_pbc(curr, others, rng) # trans other ox

        cl_idx, cl_d = find_closest(ot_wr, curr, 2)
        v1 = ot_wr[cl_idx[0]] - curr; v2 = ot_wr[cl_idx[1]] - curr
               
        d1.append(cl_d[0]); d2.append(cl_d[1]);
        t1.append(angle_between(v1[np.newaxis,:], v2[np.newaxis,:])[0]*RAD_TO_DEG)
       #print("{0} {1} {2} {3}".format(i,*curr))
       #print("{0} {1} {2} {3}".format(i,*ot_wr[cl_idx[0]]))
       #print("{0} {1} {2} {3}".format(i,*others[cl_idx[0]]))
       #print("{0} {1} {2} {3}".format(i,*ot_wr[cl_idx[1]]))

   #print([t*(180./np.pi) for t in t1])
    return d1, d2, t1


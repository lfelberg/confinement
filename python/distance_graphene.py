import sys
import re
import numpy as np

from xyzfile import XYZFile
from volfile import VolFile
from util    import d_pbc

GRAPHENE = 3

def find_closest(dist):
    '''Given an array of distances, return the smallest dist value
       and the argmin'''
    closest = np.argmin(dist)
    d_sqrt = np.sqrt(dist[closest])
    return closest, d_sqrt

def dist_btw(frm, to, rang, cond == [], pbcs = [1.0,1.0,1.0]):
    '''Given an array of molecules you would like to compute the distance 
       between, find the atom closest in the to list, given pbc range'''
    dt, at_lst = [], []
    for fm in range(len(frm)):
        frm_ar = np.repeat(frm[np.newaxis,fm,:],len(to),axis=0)
        dist = d_pbc(to, frm_ar, rang, pbcs)
        clos_at, dis_at = find_closest(dist)
        if cond != []: at_lst.append(np.where(cond)[0][clos_at]])
        dt.append(dis_at)
    return dt, cond
        
def dev_frm_avg(points):
    '''Given a list of 1D coords, calc the mean and then the signed dev from 
       that'''
    mean = np.mean(points)
    mean_ar = np.repeat(mean, len(points))
    return (points - mean_ar)

def get_all_dist(xyz, dims):
    '''Method to get the distance between graphene walls, and graphene
       and all oxygen atoms'''
    dist_btw_C, dv_g = [], []

    oxys = xyz.get_type_i(1)
    g_less = xyz.get_graph_wall(0); g_grat = xyz.get_graph_wall(1)
    
    for i in range(len(xyz.atom)):
        hlf_c = xyz.atom[i,g_less,:]; oth_c = xyz.atom[i,g_grat,:]
        oxC = xyz.atom[i,oxys,:]
        rng = dims.get_rng_i(i) # pbc range
        dist_c1, idx_c1 = dist_btw(hlf_c, oxC, rng, g_less) 
        dist_c2, idx_c2 = dist_btw(oth_c, oxC, rng, g_grat)
        dist_cs = dist_btw(, oth_c, rng, g_grat, [0.0,1.0,1.0])

        dist_to_C.append(dist_other); dist_btw_C.append(dist_cs)
        dv_g.append([dev_frm_avg(hlf_c[:,0]), dev_frm_avg(oth_c[:,0])])
    return dist_to_C, dist_btw_C, dv_g

def print_dist_to_c(fname, xyz, dists, dists_C):
    '''Print distances to carbon wall in xyz like format'''
    f = open(fname, 'w')
    for i in range(len(dists)):
        ctr, cct = 0, 0
        f.write("{0}\n".format(len(xyz.types)))
        f.write("Atoms. Timestep {0}\n".format(xyz.time[i]))
        for j in range(len(xyz.types)):
            if xyz.types[j] == GRAPHENE and cct < len(dists_C[i]):
                clo = dists_C[i][cct][1]; di = dists_C[i][cct][0]; cct += 1
            elif xyz.types[j] == GRAPHENE: clo = j; di = 0.0
            else:
                clo = dists[i][ctr][1]; di = dists[i][ctr][0]; ctr += 1
            f.write("{0} {1} {2:.4f}\n".format(xyz.types[j], clo, di))
    f.close()

def main():
    ''' usage: '''
    xyzname=sys.argv[1]; sep=sys.argv[2]; ln=sys.argv[3]; itr=sys.argv[4]
    nm = str(sep)+"_"+str(ln)+"_"+str(itr)
    volC = VolFile("run"+nm+".vol") 
    xyz_cl = XYZFile(xyzname, volC)
    dists, dists_C, grp_st = get_dist( xyz_cl, volC)

    print_dist_to_c("run"+nm+".dist",xyz_cl,dists,dists_C)

if __name__=="__main__":
    main()

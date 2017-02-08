import sys
import re
import numpy as np
from xyzfile import XYZFile
from volfile import VolFile

GRAPHENE = 3

def d_pbc(c1, c2, rng):
    '''Compute distance between points with PBCS
       c1 and c2 are an array of coordinates, rng is PBC range in x,y,z'''
    boxl = np.round((c1-c2)/rng) # find what is rounded distance
    d = c1-c2 - boxl*rng # add that rounded distance to actual dist
    ds = (d * d).sum(axis=1)
    return ds

def find_closest(dist):
    '''Given an array of distances, return the smallest dist value
       and the argmin'''
    closest = np.argmin(dist)
    d_sqrt = np.sqrt(dist[closest])
    return closest, d_sqrt

def dist_btw(frm, to, rang, cond ):
    '''Given an array of molecules you would like to compute the distance 
       between, find the atom closest in the to list, given pbc range'''
    dt = []
    for fm in range(len(frm)):
        frm_ar = np.repeat(frm[np.newaxis,fm,:],len(to),axis=0)
        dist = d_pbc(to, frm_ar, rang)
        clos_at, dis_at = find_closest(dist)
        dt.append([dis_at, np.where(cond)[0][clos_at]])
    return dt
        
def dev_frm_avg(points):
    '''Given a list of 1D coords, calc the mean and then the signed dev from 
       that'''
    mean = np.mean(points)
    mean_ar = np.repeat(mean, len(points))
    return (points - mean_ar)

def get_all_dist(xyz, dims):
    '''Method to get the distance between graphene walls, and graphene
       and all other atoms in the system'''
    grap = xyz.types == GRAPHENE; other = xyz.types != GRAPHENE; dist_to_C = []
    dist_btw_C, dv_g = [], []

    # Finding atoms with x values LT/GT half the box
    x_hlf = (dims[0][1]-dims[0][0])/2.0 # Half of box divides 2 walls
    print("In get_all_dist", xyz.atom.shape)
    x_less  = xyz.atom[0,:,0] < x_hlf; x_great = xyz.atom[0,:,0] > x_hlf
    g_less = np.all(np.array([grap, x_less]), axis=0)
    g_grat = np.all(np.array([grap, x_great]), axis=0)
    
    for i in range(len(xyz.atom)):
        # for the graphene wall-to-wall distance
        hlf_c = xyz.atom[i,g_less,:]; oth_c = xyz.atom[i,g_grat,:]
        # for the other to graphene distance
        grpC = xyz.atom[i,grap,:]; othC = xyz.atom[i,other,:]
        rng = np.array([dims[i][1]-dims[i][0], dims[i][3]-dims[i][2],
                        dims[i][5]-dims[i][4]]) # pbc range

        dist_other = dist_btw(othC, grpC, rng, grap)
        dist_cs = dist_btw(hlf_c, oth_c, rng, g_grat)
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
    ''' Given a list of pairs for g(r) calcs, do this as a function as 
        distance from wall, but only 2D and only for your wall side '''
    xyzname=sys.argv[1]; sep=sys.argv[2]; ln=sys.argv[3]; itr=sys.argv[4]
    nm = str(sep)+"_"+str(ln)+"_"+str(itr)
    volC = VolFile("run"+nm+".vol") 
    xyz_cl = XYZFile(xyzname, volC)

    print("Arry shape", xyz_cl.atom.shape)
    dists, dists_C, grp_st = get_all_dist( xyz_cl, volC.dims)

if __name__=="__main__":
    main()

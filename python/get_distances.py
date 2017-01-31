import sys
import re
import numpy as np
from get_xyz import xyz_reader

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
        
def get_all_dist(tim_xyz, coords, types, dims):
    '''Method to get the distance between graphene walls, and graphene
       and all other atoms in the system'''
    grap = types == GRAPHENE; other = types != GRAPHENE; dist_to_C = []
    dist_btw_C = []

    # Finding atoms with x values LT/GT half the box
    x_hlf = (dims[0][1]-dims[0][0])/2.0 # Half of box divides 2 walls
    x_less  = coords[0,:,0] < x_hlf; x_great = coords[0,:,0] > x_hlf
    g_less = np.all(np.array([grap, x_less]), axis=0)
    g_grat = np.all(np.array([grap, x_great]), axis=0)
    
    for i in range(len(tim_xyz)):
        hlf_c = coords[i,g_less,:]; oth_c = coords[i,g_grat,:]
        grpC = coords[i,grap,:]; othC = coords[i,other,:]
        rng = np.array([dims[i][1]-dims[i][0], dims[i][3]-dims[i][2],
                        dims[i][5]-dims[i][4]]) # pbc range

        dist_other = dist_btw(othC, grpC, rng, grap)
        dist_cs = dist_btw(hlf_c, oth_c, rng, g_grat)
        dist_to_C.append(dist_other); dist_btw_C.append(dist_cs)
    return dist_to_C, dist_btw_C

def get_dist_to_c(tim_xyz, coords, types, dims):
    '''Method to get the distances between graphene and all other particles
       from lammps xyz file'''
    grap = types == GRAPHENE; other = types != GRAPHENE; dist_to_C = []
    for i in range(len(tim_xyz)):
        grpC = coords[i,grap,:]; othC = coords[i,other,:]
        rng = np.array([dims[i][1]-dims[i][0], dims[i][3]-dims[i][2],
                        dims[i][5]-dims[i][4]])
        dis = []
        for ot in range(len(othC)): #finding closest graphene
            atomGd = np.repeat(othC[np.newaxis,ot,:],len(grpC),axis=0)
            dist = d_pbc(grpC, atomGd, rng)
            clos_C, dis_c = find_closest(dist)
            dis.append([dis_c, np.where(grap)[0][clos_C]])
        dist_to_C.append(dis)
    return dist_to_C

def print_dist_to_c(fname, time, dists, dists_C, types):
    '''Print distances to carbon wall in xyz like format'''
    f = open(fname, 'w')
    for i in range(len(dists)):
        ctr, cct = 0, 0
        f.write("{0}\n".format(len(types)))
        f.write("Atoms. Timestep {0}\n".format(time[i]))
        for j in range(len(types)):
            if types[j] == GRAPHENE and cct < len(dists_C[i]):
                clo = dists_C[i][cct][0]; di = dists_C[i][cct][1]; cct += 1
            elif types[j] == GRAPHENE: clo = j; di = 0.0
            else:
                clo = dists[i][ctr][1]; di = dists[i][ctr][0]; ctr += 1
            f.write("{0} {1} {2:.4f}\n".format(types[j], clo, di))
    f.close()

def get_box_dim(vol_out):
    '''For all times in run file, get timestamp, x, y, z min and max'''
    f = open(vol_out, "r")
    time, dims = [], []
    for line in f:
        tmp = line.split()
        time.append([int(tmp[0]), int(tmp[1])])
        dims.append([float(tmp[2]),float(tmp[3]),float(tmp[4]),
                     float(tmp[5]),float(tmp[6]),float(tmp[7])])
    f.close()
    return time, dims

def main():
    xyzname = sys.argv[1]
    sep = sys.argv[2]
    itr = sys.argv[3]
    time, dims = get_box_dim("run"+str(sep)+"_"+str(itr)+".vol")
    tim_xyz, coords, types = xyz_reader(xyzname, time, dims)
    dists, dists_C = get_all_dist(tim_xyz, np.array(coords), np.array(types), 
                                  np.array(dims))

    print_dist_to_c("run"+str(sep)+"_"+str(itr)+".dist",
                    tim_xyz,dists, dists_C,types)

if __name__=="__main__":
    main()

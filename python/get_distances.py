import sys
import re
import numpy as np
from get_xyz import xyz_reader

BENZENE = 3

def d_pbc(c1, c2, rng):
    '''Compute distance between points with PBCS'''
    diff = 0.0
    ds = 0.0
    return ds

def get_coords(filename, tim_xyz, coords, types, dims):
    '''Method to get the distances between graphene and all other particles
       from lammps xyz file'''
    grap = types == BENZENE
    print(grap)
    other = types != BENZENE
    for i in range(2) : #len(tim_xyz)):
        grp = np.array(coords[i][grap])
        print(grp.shape)
        for ot in range(2): #range(len(other)): #finding closest graphene
            atomGd = np.repeat(coords[np.newaxis,other[ot],:],len(grap),axis=0)
            print(atomGd.shape)
    
    return time, dims


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
    get_coords("run"+str(sep)+"_"+str(itr)+".dist", tim_xyz, coords,
               np.array(types), dims)


if __name__=="__main__":
    main()

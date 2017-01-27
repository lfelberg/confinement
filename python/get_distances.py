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

def get_dist_to_c(tim_xyz, coords, types, dims):
    '''Method to get the distances between graphene and all other particles
       from lammps xyz file'''
    grap = types == GRAPHENE; other = types != GRAPHENE; dist_to_C = []
    for i in range(len(tim_xyz)):
        grpC = coords[i,grap,:]
        othC = coords[i,other,:]
        rng = np.array([dims[i][1]-dims[i][0], dims[i][3]-dims[i][2],
                        dims[i][5]-dims[i][4]])
        dis = []
        for ot in range(len(othC)): #finding closest graphene
            atomGd = np.repeat(othC[np.newaxis,ot,:],len(grpC),axis=0)
            dist = d_pbc(grpC, atomGd, rng)
            clos_C = np.argmin(dist)
            closestC = np.where(grap)[0][clos_C]
           #print(clos_C, np.min(dist), closestC)
            dis.append([np.sqrt(dist[clos_C]), closestC])
        dist_to_C.append(dis)
    return dist_to_C

def get_dist_btw_cs():
    '''Method to get the distance for each graphene carbon to other wall'''
    return 0.0

def print_dist_to_c(fname, time, dists, types):
    '''Print distances to carbon wall in xyz like format'''
    f = open(fname, 'w')
    for i in range(len(dists)):
        ctr = 0
        f.write("{0}\n".format(len(types)))
        f.write("Atoms. Timestep {0}\n".format(time[i]))
        for j in range(len(types)):
            if types[j] == GRAPHENE:
                clo = j; di = 0.0
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
    dists = get_dist_to_c(tim_xyz, np.array(coords), np.array(types), 
                          np.array(dims))

    print_dist_to_c("run"+str(sep)+"_"+str(itr)+".dist",tim_xyz,dists,types)

if __name__=="__main__":
    main()

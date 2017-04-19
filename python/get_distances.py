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

def dist_btw(frm, to, rang, cond, pbcs = [1.0,1.0,1.0]):
    '''Given an array of molecules you would like to compute the distance 
       between, find the atom closest in the to list, given pbc range'''
    dt = []
    for fm in range(len(frm)):
        frm_ar = np.repeat(frm[np.newaxis,fm,:],len(to),axis=0)
        dist = d_pbc(to, frm_ar, rang, pbcs)
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
    g_less = xyz.get_graph_wall(0); g_grat = xyz.get_graph_wall(1)
    
    for i in range(len(xyz.atom)):
        # for the graphene wall-to-wall distance
        hlf_c = xyz.atom[i,g_less,:]; oth_c = xyz.atom[i,g_grat,:]
        # for the other to graphene distance
        grpC = xyz.atom[i,grap,:]; othC = xyz.atom[i,other,:]
        rng = np.array([dims[i][1]-dims[i][0], dims[i][3]-dims[i][2],
                        dims[i][5]-dims[i][4]]) # pbc range
        dist_other = dist_btw(othC, grpC, rng, grap)
        dist_cs = dist_btw(hlf_c, oth_c, rng, g_grat, [0.0,1.0,1.0])

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

def print_dev_of_c(fname, xyz, dev_g, boxl_x, sep = False):
    '''Print deviation of carbon wall to xyz like file'''
    grap = xyz.types == GRAPHENE;
    
    g_less = xyz.get_graph_wall(0); g_grat = xyz.get_graph_wall(1)
    ct = 1 if sep == True else 2
    
    for f in range(ct):
        fn = open(fname+str(f)+".dat", 'w')
        wall = g_less if f == 1 else g_grat
        for i in range(len(xyz.atom)): # for each timestep
            step = xyz.atom[i, wall, :] # find wall coords
            fn.write("{0}\nAtoms. Timestep {1}\n".format(len(step), i))
            for j in range(len(step)):
                sep_d = dev_g[i][j][0] if sep == True else dev_g[i][f][j]
                fn.write("{0} {1:.3f} {2:.3f} {l:.3f}\n".format(1, 
                         *step[j,1:], l=sep_d)) 
        fn.close()

def main():
    ''' usage: '''
    xyzname=sys.argv[1]; sep=sys.argv[2]; ln=sys.argv[3]; itr=sys.argv[4]
    nm = str(sep)+"_"+str(ln)+"_"+str(itr)
    volC = VolFile("run"+nm+".vol") 
    xyz_cl = XYZFile(xyzname, volC)
    dists, dists_C, grp_st = get_all_dist( xyz_cl, volC.dims)

    print_dist_to_c("run"+nm+".dist",xyz_cl,dists,dists_C)
    print_dev_of_c("run"+nm+"_graph", xyz_cl, grp_st, volC.get_x_rng())
    print_dev_of_c("run"+nm+"_graph_sep", xyz_cl, dists_C, 
                   volC.get_x_rng(), sep = True)

if __name__=="__main__":
    main()

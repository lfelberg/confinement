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

def dist_btw(frm, to, rang, pbcs = [1.0,1.0,1.0]):
    '''Given an array of molecules you would like to compute the distance 
       between, find the atom closest in the to list, given pbc range'''
    dt, cl_crd = [], np.zeros(frm.shape)
    for fm in range(len(frm)):
        frm_ar = np.repeat(frm[np.newaxis,fm,:],len(to),axis=0)
        dist = d_pbc(to, frm_ar, rang, pbcs)
        clos_at, dis_at = find_closest(dist)
        dt.append(dis_at); cl_crd[fm] = to[clos_at]
    return dt, cl_crd
        
def get_dist(xyz, dims):
    '''Method to get the distance between graphene walls, and graphene
       and all oxygen atoms'''
    oxys = xyz.get_type_i(1); 
    water_dists = np.zeros((3,xyz.get_nsnap(),xyz.get_ct_i(1)))
    g_less = xyz.get_graph_wall(0); g_grat = xyz.get_graph_wall(1)
    
    for i in range(len(xyz.atom)): # for each snapshot
        hlf_c = xyz.atom[i,g_less,:]; oth_c = xyz.atom[i,g_grat,:]
        oxC = xyz.atom[i,oxys,:]
        rng = dims.get_rng_i(i) # pbc range
        dist_c1, crd_c1 = dist_btw(oxC, hlf_c, rng) 
        dist_c2, crd_c2 = dist_btw(oxC, oth_c, rng)
        dist_cs = d_pbc(crd_c1, crd_c2, rng)
        water_dists[0,i] = np.array(dist_c1)
        water_dists[1,i] = np.array(dist_c2)
        water_dists[2,i] = dist_cs

    return water_dists

def print_dist_to_c(fname, dists):
    '''Print distances to carbon wall in xyz like format'''
    f = open(fname, 'w')
    st = ""; dis_n = ["_dg1", "_dg2", "_dgg"]
    for i in range(dists.shape[1]):
        for j in range(len(dis_n)): st+="{0}{1},".format(i,dis_n[j])
    f.write(st[:-1]+"\n")
   
    for k in range(dists.shape[2]):
        st = ""
        for i in range(dists.shape[1]):
            for j in range(len(dis_n)): st+="{0:.5f},".format(dists[j][i][k])
        f.write(st[:-1]+"\n")
    f.close()

def main():
    ''' usage: '''
    xyzname=sys.argv[1]; sep=sys.argv[2]; ln=sys.argv[3]; itr=sys.argv[4]
    nm = str(sep)+"_"+str(ln)+"_"+str(itr)
    volC = VolFile("run"+nm+".vol") 
    xyz_cl = XYZFile(xyzname, volC)
    dists = get_dist( xyz_cl, volC)
    print_dist_to_c("run"+nm+".distgg", dists)

if __name__=="__main__":
    main()

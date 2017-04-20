import sys
import itertools 
import numpy as np

from xyzfile import XYZFile
from volfile import VolFile
from util    import d_pbc, translate_pbc

GRAPHENE = 3
SPAC = 2.5

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
    oxys,_ = xyz.get_inner_wat(); w_d1, w_d2 = [], []; cc_d = []
    water_dists = np.zeros((3,xyz.get_nsnap()-1,xyz.get_ct_i(1)))
    c1_x = np.zeros((1,1)); c2_x = np.zeros((1,1))
    g_less = xyz.get_graph_wall(0); g_grat = xyz.get_graph_wall(1)
    yedg = np.arange(0,dims.get_y_len(),SPAC); n_y_ed = len(yedg)
    zedg = np.arange(0,dims.get_z_len(),SPAC); n_z_ed = len(zedg)
    
    for i in range(1,len(xyz.atom)): # for each snapshot
        c1 = xyz.atom[i,g_less,:]; c2 = xyz.atom[i,g_grat,:]
        oxC = xyz.atom[i,oxys,:]; rng = dims.get_rng_i(i) # pbc range
        c1_yy = np.digitize(c1[:,1], yedg);c2_yy = np.digitize(c2[:,1], yedg)
        c1_zz = np.digitize(c1[:,2], zedg);c2_zz = np.digitize(c2[:,2], zedg)
        o_yy  = np.digitize(oxC[:,1], yedg);o_zz  = np.digitize(oxC[:,2], zedg)
        for yy in range(n_y_ed):
            for zz in range(n_z_ed):
                o_p = np.all(np.array([o_yy==yy, o_zz==zz]),axis=0)
                if sum(o_p.astype(int)) > 0:
                    c1_p = np.all(np.array([c1_yy==yy, c1_zz==zz]),axis=0)
                    c1_cs = c1[c1_p,0];
                    c1_xp = translate_pbc(np.zeros((1,1)),c1_cs,rng[0])
                    c1_x[0,0] = np.mean(c1_xp)
                    c2_p = np.all(np.array([c2_yy==yy, c2_zz==zz]),axis=0)
                    c2_cs = c2[c2_p,0];
                    c2_xp = translate_pbc(np.array(rng[0]),c2_cs,rng[0])
                    c2_x[0,0] = np.mean(c2_xp)
                    d1 = list(d_pbc(oxC[o_p,0][:,np.newaxis],c1_x,rng[0],[1.]))
                    d2 = list(d_pbc(oxC[o_p,0][:,np.newaxis],c2_x,rng[0],[1.]))
                    w_d1.append(d1); w_d2.append(d2)
                    dgg = len(d1) * list((c2_x - c1_x)[0])
                    cc_d.append(dgg)
    w_d1 = list(itertools.chain(*w_d1)); w_d2 = list(itertools.chain(*w_d2))
    return w_d1, w_d2, list(itertools.chain(*cc_d))

def print_dist_to_c(fname, dists):
    '''Print distances to carbon wall in xyz like format'''
    f = open(fname, 'w')
    st = ""; dis_n = ["dg1", "dg2", "dgg"]
    for j in range(len(dis_n)): st+="{0},".format(dis_n[j])
    f.write(st[:-1]+"\n")
   
    for k in range(len(dists[0])):
        st = ""
        for i in range(len(dis_n)): st+="{0:.5f},".format(dists[i][k])
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

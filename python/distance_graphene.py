import sys
import itertools 
import numpy as np

from xyzfile import XYZFile
from volfile import VolFile
from util    import d_pbc, translate_pbc, translate_1st_im

GRAPHENE = 3
SPAC = 2.5

def get_dist(xyz, dims):
    '''Method to get the distance between graphene walls, and graphene
       and all oxygen atoms'''
    oxys,_ = xyz.get_inner_wat(); w_d1, w_d2 = [], []; cc_d = []
    water_dists = np.zeros((3,xyz.get_nsnap()-1,xyz.get_ct_i(1)))
    g_less = xyz.get_graph_wall(0); g_grat = xyz.get_graph_wall(1)
    yedg = np.arange(0,dims.get_y_len(),SPAC); n_y_ed = len(yedg)
    zedg = np.arange(0,dims.get_z_len(),SPAC); n_z_ed = len(zedg)
    
    for i in range(1,len(xyz.atom)): # for each snapshot
        c1 = xyz.atom[i,g_less,:]; c2 = xyz.atom[i,g_grat,:]
        oxC = xyz.atom[i,oxys,:]; rng = dims.get_rng_i(i) # pbc range
        c1[:,0] = translate_pbc(np.array(0.5),c1[:,0],rng[0]) 
        c2[:,0] = translate_pbc(np.array(rng[0]-0.5),c2[:,0],rng[0]) 
        c1_yy = np.digitize(c1[:,1],yedg); c2_yy = np.digitize(c2[:,1],yedg)
        c1_zz = np.digitize(c1[:,2],zedg); c2_zz = np.digitize(c2[:,2],zedg)
        o_yy  = np.digitize(oxC[:,1],yedg); o_zz  = np.digitize(oxC[:,2],zedg)
        for yy in range(n_y_ed):
            for zz in range(n_z_ed):
                o_p = np.all(np.array([o_yy==yy, o_zz==zz]),axis=0)
                if sum(o_p.astype(int)) > 0:
                    o_w = translate_1st_im(oxC[o_p,0], rng[0])
                    c1_p = c1[np.all(np.array([c1_yy==yy,c1_zz==zz]),axis=0),0]
                    c2_p = c2[np.all(np.array([c2_yy==yy,c2_zz==zz]),axis=0),0]
                    c1_x = np.mean(c1_p); c2_x = np.mean(c2_p)
                    d1 = list(abs(o_w-c1_x));  d2 = list(abs(o_w-c2_x))
                    if c2_x - c1_x < 2.0: print("This is pts and means",  c1_p, c2_p, c1_x, c2_x)
                    w_d1.append(d1); w_d2.append(d2)
                    dgg = len(d1) * [abs(c2_x - c1_x)]
                    cc_d.append(dgg)
                else:
                    c1_p = c1[np.all(np.array([c1_yy==yy,c1_zz==zz]),axis=0),0]
                    c2_p = c2[np.all(np.array([c2_yy==yy,c2_zz==zz]),axis=0),0]
                    w_d1.append(0.0); w_d2.append(0.0)
                    cc_d.append([abs(c2_x - c1_x)])
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

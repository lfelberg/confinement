import sys
import itertools 
import numpy as np

from xyzfile import XYZFile
from volfile import VolFile
from util    import d_pbc, translate_pbc, translate_1st_im

GRAPHENE = 3; OXY = 1
SPAC = 3.5

def get_dist(xyz, dims):
    '''Method to get the distance between graphene walls, and graphene
       and all oxygen atoms'''
    oxys = xyz.get_type_i(OXY); w_d1, w_d2 = [], []; cc_d = []
    water_dists = np.zeros((3,xyz.get_nsnap()-1,xyz.get_ct_i(OXY)))
    g_less = xyz.get_graph_wall(0); g_grat = xyz.get_graph_wall(1)
    yedg = np.arange(0,dims.get_y_len(),SPAC); n_y_ed = len(yedg)
    zedg = np.arange(0,dims.get_z_len(),SPAC); n_z_ed = len(zedg)
    
    onxon = np.zeros((1,1)) # 1X1 array used for translating when array inp
    for i in range(1,len(xyz.atom)): # for each snapshot
        rng=dims.get_rng_i(i)
        c1 = xyz.atom[i,g_less]; c2=xyz.atom[i,g_grat]; oxC=xyz.atom[i,oxys] 

        # for both walls, check if they are wrapped
        mnc = min(c1[:,0]); mxc = max(c1[:,0])
        if mxc - mnc > 2:
            onxon[0,0] = mnc; c1[:,0] = translate_pbc(onxon, c1[:,0], rng[0])

        mnc = min(c2[:,0]); mxc = max(c2[:,0])
        if mxc - mnc > 2:
            onxon[0,0] = mnc; c2[:,0] = translate_pbc(onxon, c2[:,0], rng[0])

        m1 = np.mean(c1[:,0]); m2 = np.mean(c2[:,0])
        if m2 < m1 : continue # skip frame if second wall now comes before first

        #move all so that 1st wall is at x=0
        c1[:,0] -= m1; c2[:,0] -= m1; oxC[:,0] -= m1; onxon[0,0] = 0.0
        oxC[:,0] = translate_pbc(onxon, oxC[:,0], rng[0])

        c1_yy = np.digitize(c1[:,1],yedg); c2_yy = np.digitize(c2[:,1],yedg)
        c1_zz = np.digitize(c1[:,2],zedg); c2_zz = np.digitize(c2[:,2],zedg)
        o_yy  = np.digitize(oxC[:,1],yedg); o_zz  = np.digitize(oxC[:,2],zedg)
        for yy in range(n_y_ed):
            for zz in range(n_z_ed):
                o_p = np.all(np.array([o_yy==yy, o_zz==zz]),axis=0)
                c1_p = sum(np.all(np.array([c1_yy==yy,c1_zz==zz]),axis=0).astype(int))
                c2_p = sum(np.all(np.array([c1_yy==yy,c1_zz==zz]),axis=0).astype(int))
                if sum(o_p.astype(int)) > 0:
                    o_w = oxC[o_p,0]
                    c1_p = c1[np.all(np.array([c1_yy==yy,c1_zz==zz]),axis=0),0]
                    c2_p = c2[np.all(np.array([c2_yy==yy,c2_zz==zz]),axis=0),0]
                    c1_x = np.mean(c1_p); c2_x = np.mean(c2_p)
                    d1 = list(o_w-c1_x); w_d1.append(d1)
                    dgg = len(d1) * [abs(c2_x - c1_x)]; cc_d.append(dgg)
 
                # this part was used ONLY for the 6A flexible case where the gg
                # are in contact
               #elif c1_p > 0 and c2_p > 0:
               #    c1_p = c1[np.all(np.array([c1_yy==yy,c1_zz==zz]),axis=0),0]
               #    c2_p = c2[np.all(np.array([c2_yy==yy,c2_zz==zz]),axis=0),0]
               #    c1_x = np.mean(c1_p); c2_x = np.mean(c2_p)
               #    w_d1.append([0.0]); w_d2.append([0.0])
               #    cc_d.append([abs(c2_x - c1_x)])
    return list(itertools.chain(*w_d1)), list(itertools.chain(*cc_d))

def print_dist_to_c(fname, dists):
    '''Print distances to carbon wall in csv like format'''
    f = open(fname, 'w'); st = ""; dis_n = ["dg1", "dgg"]
    for j in range(len(dis_n)): st+="{0},".format(dis_n[j])
    f.write(st[:-1]+"\n")
   
    for k in range(len(dists[0])):
        st = ""
        for i in range(len(dis_n)): st+="{0:.5f},".format(dists[i][k])
        f.write(st[:-1]+"\n")
    f.close()

def main():
    ''' usage: python distance_graphene.py xyzname sep len iter'''
    xyzname=sys.argv[1]; sep=sys.argv[2]; ln=sys.argv[3]; itr=sys.argv[4]
    nm = str(sep)+"_"+str(ln)+"_"+str(itr)
    volC = VolFile("run"+nm+".vol"); xyz_cl = XYZFile(xyzname, volC)
    dists = get_dist( xyz_cl, volC); print_dist_to_c("run"+nm+".distgg", dists)

if __name__=="__main__":
    main()

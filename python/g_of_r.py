import sys
import re
import math
import numpy as np

from xyzfile import XYZFile
from volfile import VolFile
from util    import d_pbc

#TODO: Find out if doing 2D g(r) is worth it?

GRAPHENE = 3
BIN_GR = 80
LMAX = 12.0

def gr_cal(crd0, crd1, rang, xbins, dists):
    '''Given an array of molecules you would like to compute the distance 
       between, find their distances, histogram them and create hist based
       on distance of crd0 from graphene. Will do 2D or 3D, depending on the
       length of the rang array.'''
    # Here x... is for a histogram of the type 0 particles dist from wall.
    gr, dim, vl = [], len(rang), 1.0; st = 3 - dim; xst = xbins[1]-xbins[0]
    xmx = xbins[-1]+xst; dr = LMAX/float(BIN_GR)
    # Normalization factors for g(r)
    rd_rng = np.arange(0, LMAX, dr); low, upp = rd_rng[:-1], rd_rng[1:]
    if dim == 3: rd_mu = 4.0/3.0*np.pi*(np.power(upp,3.) - np.power(low,3.))
    else: rd_mu = 2.0*np.pi*(np.power(upp,2.) - np.power(low,2.)) 

    # number density normalization
    n0 = float(len(crd0))
    n1 = n0 if crd1 == [] else float(len(crd1))
    for d in range(st, 3): vl *= rang[d]
    num_dens = n0*n1 / vl

    # For storing the data
    at_ct = np.zeros((len(xbins))); gr_ar = np.zeros((len(xbins), BIN_GR))
    for fm in range(len(crd0)): # for each of type 0, find dist to each type 1
        to_ar = crd1[:,st:] if crd1 != [] else np.delete(crd0[:,st:],fm,0)
        frm_ar = np.repeat(crd0[np.newaxis,fm,st:],len(to_ar),axis=0)
        dist = d_pbc(to_ar, frm_ar, rang)
        his_den, benz = np.histogram(dist, BIN_GR, range=(0, LMAX))
        x_loc = math.floor((dists[fm]/xmx)*float(len(xbins))) # find idx of cor0[fm]
        at_ct[x_loc] += 1.0
        gr_ar[x_loc] += (his_den.astype(float)/rd_mu)

    at_ct[at_ct==0] = 1.0 # if no at in bin, set to 1, for no div/0
    at_ar = np.repeat(at_ct[:,np.newaxis],gr_ar.shape[1],axis=1)
    return gr_ar/at_ar
        
def get_gr(xyz, disC, dists, volC, grPair):
    '''Method to get the g(r) for two atom types'''
    # find atoms of type 0 and type 1
    ty0 = xyz.get_type_i(grPair[0]); ty1 = xyz.get_type_i(grPair[1]); bnz = 10
    g_r_3, g_r_2_m, rng_m = [], [], np.zeros((3))

    x_hlf = volC.get_x_max()/2.0 # Half of box divides 2 walls
    # for x dir, use bins 0 -> x_hlf. ATs will be max L/2 from wall
    x_binz = np.arange(2.5, x_hlf, x_hlf/float(bnz)); 

    # for the g(r) pairs, need to get groups on the same wall side
    ty00 = np.all(np.array([ty0, xyz.get_inner_ats()]), axis=0) # btw 2 walls
    ty10 = np.all(np.array([ty1, xyz.get_inner_ats()]), axis=0)

    ty01 = np.all(np.array([ty0, xyz.get_outer_ats()]), axis=0) # out 2 walls
    ty11 = np.all(np.array([ty1, xyz.get_outer_ats()]), axis=0)

   #for i in range(1,len(xyz.atom)):
    for i in range(len(xyz.atom)):
        rng = np.array([volC.get_x_rng_i(i), volC.get_y_rng_i(i),
                        volC.get_z_rng_i(i)]) # pbc range
        rng_m += rng # take average of range for printing

        c00 = xyz.atom[i,ty00,:]; # coords for them
        if grPair[0] == grPair[1]: c10 = [] 
        else: c10 = xyz.atom[i,ty10,:]
        g_r_3.append(gr_cal(c00, c10, rng, x_binz, dists[i,ty00,-1]))
       #print(c00.shape, c10.shape, x_binz.shape, dists[i,ty00,-1].shape)

        c01 = xyz.atom[i,ty01,:]; 
        if grPair[0] == grPair[1]: c11 = [] 
        else: c11 = xyz.atom[i,ty10,:]
        g_r_3.append(gr_cal(c01, c11, rng, x_binz, dists[i,ty01,-1]))

    # Time averages of histograms
    g_r_3 = np.array((g_r_3))
    g_r_3_m = np.mean(g_r_3, axis = 0)
    return g_r_3_m, g_r_2_m, x_binz, rng_m/float(len(xyz.atom))

def print_gr(x, grs, fname):
    '''Print distances to carbon wall in xyz like format'''
    f = open(fname, 'w'); 
    f.write("Bin,"); st = ""
    dimbins = np.arange(0, LMAX, LMAX/float(BIN_GR))
    for i in range(len(x)): st += "{0:.5f},".format(x[i])
    f.write("{0}\n".format(st[:-1]))
    for i in range(BIN_GR):
        st = ""
        for j in range(len(grs[i])): st += "{0:.5f},".format(grs[i][j])
        f.write("{0:.4f},{1}\n".format(dimbins[i], st[:-1]))
    f.close()

def main():
    ''' Given a list of pairs for g(r) calcs, do this as a function as 
        distance from wall, but only 2D and only for your wall side '''
    xyzname=sys.argv[1]; sep=sys.argv[2]; ln=sys.argv[3]; itr=sys.argv[4]
    n_gr = int(sys.argv[5]); grS, grPr = 6, []
    for i in range(grS, grS + n_gr*2, 2): 
        grPr.append([int(sys.argv[i]),int(sys.argv[i+1])])

    nm = str(sep)+"_"+str(ln)+"_"+str(itr)
    volC = VolFile("run"+nm+".vol") 
    disC = XYZFile("run"+nm+".dist", VolFile("")) 
    xyz_cl = XYZFile(xyzname, volC)
    print("Arry shape", xyz_cl.atom.shape, disC.atom.shape)

    for i in range(n_gr):
        prNm = '_'+str(grPr[i][0])+'_'+str(grPr[i][1])+'.csv'
        grs3D, grs2D, xrng, rng_m = get_gr(xyz_cl, disC, disC.atom, volC, grPr[i])
        print_gr(xrng, grs3D.T, 'g_r_3D_'+nm+prNm)

if __name__=="__main__":
    main()

import sys
import re
import math
import numpy as np
from xyzfile import XYZFile
from volfile import VolFile

GRAPHENE = 3
BIN_GR = 80
LMAX = 12.0

def d_pbc(c1, c2, rng):
    '''Compute distance between points with PBCS
       c1 and c2 are an array of coordinates, rng is PBC range in x,y,z'''
    boxl = np.round((c1-c2)/rng) # find what is rounded distance
    d = c1-c2 - boxl*rng # add that rounded distance to actual dist
    ds = (d * d).sum(axis=1)
    return np.sqrt(ds)

def gr_cal(crd0, crd1, rang, xbins, dists):
    '''Given an array of molecules you would like to compute the distance 
       between, find their distances, histogram them and create hist based
       on distance of crd0 from graphene. Will do 2D or 3D, depending on the
       length of the rang array.'''
    # Here x... is for a histogram of the type 0 particles dist from wall.
    gr, dim = [], len(rang); st = 3 - dim; xst = xbins[1]-xbins[0]
    xmx = xbins[-1]+xst; dr = LMAX/float(BIN_GR)
    # Normalization factors for g(r)
    rd_rng = np.arange(0, LMAX, dr)
    rd_mu = 4.*math.pi*rd_rng*rd_rng*dr if dim == 3 else 2.*math.pi*rd_rng*dr
    rd_mu[0] = 1.0 # set norm of first pos to zero so no runtime errors

    # For storing the data
    at_ct = np.zeros((len(xbins))); gr_ar = np.zeros((len(xbins), BIN_GR))
    for fm in range(len(crd0)): # for each of type 0, find dist to each type 1
        frm_ar = np.repeat(crd0[np.newaxis,fm,st:],len(crd1),axis=0)
        dist = d_pbc(crd1[:,st:], frm_ar, rang)
        his_den, benz = np.histogram(dist, BIN_GR, range=(0, LMAX))
        x_loc = math.floor((dists[fm]/xmx)*float(len(xbins))) # find idx of cor0[fm]
        at_ct[x_loc] += 1.0
        gr_ar[x_loc] += (his_den.astype(float)/rd_mu)

    at_ct[at_ct==0] = 1.0 # if no at in bin, set to 1, for no div/0
    at_ar = np.repeat(at_ct[:,np.newaxis],gr_ar.shape[1],axis=1)
   #print("This is bins", benz)
    return gr_ar/at_ar
        
def get_gr(xyz, disC, dists, volC, grPair):
    '''Method to get the g(r) for two atom types'''
    # find atoms of type 0 and type 1 and type graphene
    ty0 = xyz.types == grPair[0]; ty1 = xyz.types == grPair[1]; bnz = 5
    g_r_3, g_r_2, rng_m = [], [], np.zeros((3)); grap = xyz.types == GRAPHENE

    x_hlf = volC.get_x_max()/2.0 # Half of box divides 2 walls
    # for x dir, use bins 1.4 -> x_quarter. ATs will be max L/4 from wall
    x_qu = x_hlf/2.0; x_binz = np.arange(2.5, x_qu, x_qu/float(bnz)); 
    # Finding graphene with x values LT/GT half box, will delimit g_r
    xl  = xyz.atom[0,:,0] < x_hlf; xgr = xyz.atom[0,:,0] > x_hlf
    g_less = np.all(np.array([grap, xl]), axis=0)
    g_grat = np.all(np.array([grap, xgr]), axis=0)

    for i in range(len(xyz.atom)):
        # for locations of graphene walls find average X value
        w0 = np.mean(xyz.atom[i,g_less,0]); w1 = np.mean(xyz.atom[i,g_grat,0])

        rng = np.array([volC.get_x_rng_i(i), volC.get_y_rng_i(i),
                        volC.get_z_rng_i(i)]) # pbc range
        rng_m += rng # take average of range for printing

        # for the g(r) pairs, need to get groups on the same wall side
        x_m0  = xyz.atom[i,:,0] > w0; x_m1 = xyz.atom[i,:,0] < w1
        wall = np.all(np.array([x_m0, x_m1]), axis=0); nt_wl = wall == False
        ty00 = np.all(np.array([ty0, wall]), axis=0) # between 2 walls
        ty10 = np.all(np.array([ty1, wall]), axis=0)
        c00 = xyz.atom[i,ty00,:]; c10 = xyz.atom[i,ty10,:] # coords for them

        print(c00.shape, c01.shape, x_binz.shape, dists[i,ty00,-1].shape)


        g_r_3.append(gr_cal(c00, c10, rng, x_binz, dists[i,ty00,-1]))
       #g_r_2.append(gr_cal(c00, c10, rng[1:], x_binz, dists[i,ty00, -1]))

        ty01 = np.all(np.array([ty0, nt_wl]), axis=0) # outside of 2 walls
        ty11 = np.all(np.array([ty1, nt_wl]), axis=0)
        c01 = xyz.atom[i,ty01,:]; c11 = xyz.atom[i,ty11,:] # coords for them
        g_r_3.append(gr_cal(c01, c11, rng, x_binz, dists[i,ty01,-1]))

    # Time averages of histograms
    g_r_2 = np.array((g_r_2)); g_r_3 = np.array((g_r_3))
    g_r_2_m = g_r_2; #np.mean(g_r_2, axis = 0); 
   #g_r_2_m = np.mean(g_r_2, axis = 0); 
    g_r_3_m = np.mean(g_r_3, axis = 0)

    print(g_r_2.shape, g_r_3.shape)
    return g_r_3_m, g_r_2_m, x_binz, rng_m/float(len(xyz.atom))

def print_gr(x, grs, fname):
    '''Print distances to carbon wall in xyz like format'''
    print(grs.shape)
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
    print("Gr pairs", grPr)

    for i in range(n_gr):
        prNm = '_'+str(grPr[i][0])+'_'+str(grPr[i][1])+'.csv'
        grs3D, grs2D, xrng, rng_m = get_gr(xyz_cl, disC, disC.atom, volC, grPr[i])
       #print_gr(xrng, grs2D.T, 'g_r_2D_'+nm+prNm)
        print_gr(xrng, grs3D.T, 'g_r_3D_'+nm+prNm)

if __name__=="__main__":
    main()

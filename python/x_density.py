import sys
import csv
import numpy as np

from util import translate_pbc

from xyzfile import XYZFile
from volfile import VolFile
from csvfile import CSVFile

GRAPHENE = 3
OXY = 1; HYD = 2
A3_TO_CM3 = 1.e-24     # 1 A^3 = 1.0e-24  cm^3
AMU_TO_GM = 1.6605e-24 # 1 amu = 1.66e-24 grams

type_wt = {
            1: 16.00,
            2:  1.01,
            3: 12.00,
            4: 12.00,
            5:  1.01,
          }

def histog_dist(xyzC, volC):
    ''' Get a histogram of distances and transform into distance from plate'''
    xm = np.ceil(volC.get_x_max()/2)
    h_rng = (-xm,xm); bns = int((2.*xm)/0.1)
   #bns = 100
    kk = 0

    # use wall atoms as reference points
    w0 = xyzC.get_graph_wall(0); w1 = xyzC.get_graph_wall(1)
    ox = xyzC.get_type_i(OXY); hy = xyzC.get_type_i(HYD)
    nox = xyzC.get_ct_i(OXY); nhy = xyzC.get_ct_i(HYD)
    
    w0_wrap = np.zeros((1,1)); c_wl = np.zeros((3,1))

    # histogram of density for each type
    dens = np.zeros((4, len(xyzC.atom), bns))
    for j in range(1, len(xyzC.atom)):
        coords = xyzC.atom[j,:,0]; xj = volC.get_x_rng_i(j)
        rngj = volC.get_rng_i(j)
        # for each snap, find wall x loc, make sure wall isn't wrapped
        w0_m=np.mean(coords[w0]); w0_s=np.std(coords[w0])
        w1_m=np.mean(coords[w1]); w1_s=np.std(coords[w1])

       #print(" w1 mean: {0:.2f}, w1 std: {1:.2f} and w2 m: {2:.2f}, w2 s: {3:.3f}".
       #      format(w0_m, w0_s, w1_m, w1_s))

        if w0_s > 0.25 or w1_s > 0.25: continue #wall is right at edge, so skip

        coords = coords - w0_m # recenter so w0 is at x = 0
        w_o = translate_pbc(w0_wrap, coords[ox], xj)
        dens[0,kk],be=np.histogram(w_o, bins=bns, range=h_rng); dx=be[1]-be[0]
        dens[0,kk] /= float(len(w_o)*dx)
        w_h = translate_pbc(w0_wrap, coords[hy], xj)
        dens[1,kk], be = np.histogram(w_h, bins=bns, range=h_rng)
        dens[1,kk] /= float(len(w_h)*dx)

        c_sol,_ = xyzC.get_sol_crd_i(j,rngj,False) 
        sl = translate_pbc(w0_wrap, c_sol[:,0] - w0_m, xj)
        dens[2,kk], be = np.histogram(sl, bins=bns, range=h_rng)
        dens[2,kk] /= float(len(sl)*dx)

        c_wl[1] = w1_m - w0_m
        if w1_m > w0_m: c_wl[2] = c_wl[1] - xj
        else:           c_wl[2] = c_wl[1] + xj
        dens[3,kk], be = np.histogram(c_wl, bins=bns, range=h_rng)
        dens[3,kk] /= float(3.0*dx)
        kk += 1

   #dens *= (AMU_TO_GM/A3_TO_CM3)
    print(dens[:,:kk].shape)
    dens_mn = np.mean(dens[:,:kk], axis = 1)
    
    return be[1:]-(dx/2), dens_mn

def print_hist(hist_out,dens_mn, benz):
    f, headr = open(hist_out, 'w'), 'Bin'
    ty_ln = [ 'w_o', 'w_h', 'sol', 'g' ]
    for i in range(len(ty_ln)): headr += (","+ty_ln[i])
    f.write(headr+"\n")
    for i in range(len(benz)):
        st = "{0:8f},".format(benz[i])
        for j in range(len(ty_ln)): st += "{0:.5f},".format(dens_mn[j,i])
        f.write(st[:-1]+"\n")
    f.close()

def main():                                                                        
    '''Uses volume file and XYZ or CSV (ext: *dens_hist.csv) file'''
    xyzname=sys.argv[1]; sep=sys.argv[2]; ln=sys.argv[3]; itr=sys.argv[4]                    
    sol_nm = sys.argv[5]; sol_ct = int(sys.argv[6])
    nm = str(sep)+"_"+str(ln)+"_"+itr
    volC = VolFile("run"+nm+".vol")
    xyzC = XYZFile(xyzname, volC)
    xyzC.sol_ty = sol_nm
    xyzC.nsol = sol_ct
    nm = str(sep)+"_"+str(ln)+"_"+str(sol_ct/2)+"_"+itr+"_"
    bb, dm = histog_dist(xyzC, volC)
    print_hist("run"+nm+"dens_hist.csv", dm, bb)

if __name__=="__main__":
    main() 

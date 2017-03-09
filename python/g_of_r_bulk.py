import sys
import re
import math
import numpy as np

from xyzfile import XYZFile
from volfile import VolFile
from util    import d_pbc

DR = 0.1
BNS = np.arange( 0., 20., DR)
RADI = DR/2.+BNS[:-1] # center of each histogram bin

def gr_cal(rang, crd0, crd1 = []):
    '''Given an array of molecules you would like to compute the distance 
       between, find their distances, histogram them 
       If crd1 is empty, then we are G of R with the same type
    '''
    gr, dim, vl, n0 = [], len(rang), 1.0, float(len(crd0)); st = 3 - dim;
    # Normalization factors for g(r)
    for d in range(st, 3): vl *= rang[d]
    n1 = n0 if crd1 == [] else float(len(crd1))
    num_dens = n0*n1 / vl; upp, low = BNS[1:], BNS[:-1]
    
    if dim == 3: nrm = 4.0/3.0*np.pi*(np.power(upp, 3.) - np.power(low,3.))
    else:        nrm = 2.0*np.pi*(np.power(upp, 2.) - np.power(low,2.))
    print(nrm, vl, rang)

    # For storing the data
    his_all = np.zeros(len(RADI))
    end = len(crd0) if crd1 != [] else len(crd0)-1
    for fm in range(end): # for each of type 0, find dist to each type 1
        to_arr = crd1[:,st:] if crd1 != [] else np.delete(crd0[:,st:],fm,0)
        frm_ar = np.repeat(crd0[np.newaxis,fm,st:],len(to_arr),axis=0)
        dist = d_pbc(to_arr, frm_ar, rang)
        his_den, benz = np.histogram(dist, bins=BNS)
        his_all += his_den
       #print(len(dist), to_arr.shape, frm_ar.shape, benz)
    print(num_dens, BNS, len(crd0))

    return his_all/ nrm / num_dens
        
def get_gr(xyz, volC, grPair):
    '''Method to get the g(r) for two atom types'''
    # find atoms of type 0 and type 1
    ty0 = xyz.get_type_i(grPair[0]); ty1 = xyz.get_type_i(grPair[1])
    g_r_3, g_r_2_m = [], []
    
    for i in range(len(xyz.atom)):
        rng = np.array([volC.get_x_rng_i(i), volC.get_y_rng_i(i),
                        volC.get_z_rng_i(i)]) # pbc range

        c0 = xyz.atom[i,ty0,:]; # coords for each atom pair
        if grPair[0] == grPair[1]: c1 = []
        else: c1 = xyz.atom[i,ty1,:]
        g_r_3.append(gr_cal(rng, c0, c1))
        
        print(len(c0), len(c1), )

    # Time averages of histograms
    g_r_3 = np.array((g_r_3))
    g_r_3_m = np.mean(g_r_3, axis = 0)
    print(g_r_3.shape,g_r_3_m.shape)
    return g_r_3_m, 

def print_gr(grs, fname):
    '''Print time averaged g-r'''
    f = open(fname, 'w'); 
    f.write("Bin,GR\n")
    for i in range(len(RADI)):
        f.write("{0:.4f},{1:.4f}\n".format(RADI[i], grs[i]))
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
    xyz_cl = XYZFile(xyzname, volC)
    print("Arry shape", xyz_cl.atom.shape)

    for i in range(n_gr):
        prNm = '_'+str(grPr[i][0])+'_'+str(grPr[i][1])+'.csv'
        grs3D,  = get_gr(xyz_cl, volC, grPr[i])
        print_gr(grs3D, 'g_r_3D_'+nm+prNm)

if __name__=="__main__":
    main()

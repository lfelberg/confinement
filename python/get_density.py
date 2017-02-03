import sys
import numpy as np
import matplotlib.pyplot as plt  

from xyzfile import XYZFile
from volfile import VolFile

GRAPHENE = 3
A3_TO_CM3 = 1.e-24     # 1 A^3 = 1.0e-24  cm^3
AMU_TO_GM = 1.6605e-24 # 1 amu = 1.66e-24 grams

type_wt = {
            1: 16.00,
            2:  1.01,
            3: 12.00,
            4: 12.00,
            5:  1.01,
          }


def histog_dist(hist_out, xyzC, volC):
    ''' Get a histogram of distances and transform into distance from plate'''
    xr, bns = volC.get_x_len(), 100
    type_lst = type_wt.keys()
    ty_ln = len(type_lst)

    # histogram of density for each type
    dens = np.zeros((ty_ln, len(xyzC.atom), bns))
    for j in range(len(xyzC.atom)):
        ty_ct = 0
        vol = (xr/bns)*volC.get_y_rng_i(j)*volC.get_z_rng_i(j)
        for t in range(ty_ln):
            sub_ar = xyzC.types == type_lst[t]
            # getting x coords of atoms type ty
            cords = xyzC.atom[j,sub_ar,0]
            his_den, benz = np.histogram(cords, bns, range=(0, xr))
            dens[ty_ct][j] = his_den.astype(float)/vol
            dens[ty_ct][j] *= type_wt[type_lst[t]]
            ty_ct += 1
    dens *= (AMU_TO_GM/A3_TO_CM3)
    dens_mn = np.mean(dens, axis = 1)

    f = open(hist_out, 'w')
    f.write("Bin,type1,type2,type3,type4,type5\n")
    for i in range(bns):
        f.write("{0:>5.3f},{1:.4f},{2:.4f},{3:.4f},{4:.4f}\n".format(
                 benz[i], *tuple(dens_mn[:,i])))
    f.close()


def plot_density_hist(fname):
    d_sums = np.zeros((3,dens_mn.shape[1]))
    d_sums[0] = dens_mn[0]+dens_mn[1]
    d_sums[1] = dens_mn[2]
    d_sums[2] = dens_mn[3]+dens_mn[4]

    for t in range(ty_ln):
        f = plt.figure(1, figsize = (3.0, 3.0))                                    
        ax = f.add_subplot(111)    
        ax.plot(benz[:-1], dens_mn[t])
        plt.savefig('dens'+str(type_lst[t])+'.png', format='png',                            
                            bbox_inches = 'tight', dpi=300) 
        plt.close()

    f = plt.figure(1, figsize = (3.0, 3.0))                                    
    ax = f.add_subplot(111)    
    ax.plot(benz[:-1], d_sums[0])
    plt.savefig('dens_wat'+'.png', format='png',                            
                            bbox_inches = 'tight', dpi=300) 
    plt.close()

def main():                                                                        
    xyzname = sys.argv[1]; sep = sys.argv[2]; itr = sys.argv[3]                    
    volC = VolFile("run"+str(sep)+"_"+str(itr)+".vol")                             
    xyzC = XYZFile(xyzname, volC)

    histog_dist(xyzname[:-3]+"dens_hist", xyzC, volC)
                                                                                   
if __name__=="__main__":                                                           
    main() 

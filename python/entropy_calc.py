import sys
import re
import math
import numpy as np
import itertools
import matplotlib.pyplot as plt 

from scipy.integrate import simps

from csvfile import CSVFile
from volfile import VolFile
from util    import d_pbc

NANGLES = 5
KB = 1.38064852e-23   # joules/kelvin
MOL = 6.022140857e23  # avogadro's number
J_TO_CAL = 1.0/4.1868 # 1 calorie/4.1868 Joules 
KB_CAL_MOL = KB * J_TO_CAL * MOL

def g_of_r(bns, dist, density = 1.0):
    ''' Given an array of distances, histogram to compute g(r) '''
    hist, bins = np.histogram(dist, bins = bns)
    upp, low = bns[1:], bns[:-1]
    ndens = float(len(dist)) / density
    gr = hist/(4.0/3.0*np.pi*(np.power(upp, 3.) - np.power(low,3.)))/ndens
    return gr, ndens

def g_of_r_multi(bns, dist, density = 1.0):
    ''' Given an array of distances, histogram to compute g(r) '''
    hist, bins = np.histogramdd(dist, bins = bns)
    upp, low = bns[1:], bns[:-1]
    ndens = float(len(dist)) / density
    gr = hist/(4.0/3.0*np.pi*(np.power(upp, 3.) - np.power(low,3.)))/ndens
    return gr, ndens

def trans_entropy(dists, vl):
    '''Compute the translational entropy from a list of pairwise distances
    '''
    ntimes, npairs = dists.shape[0], dists.shape[1]
    binsiz_r, grs, ndens = 0.53, [], np.zeros(ntimes)
    bns_r = np.arange( 2.5, 15., binsiz_r)
    radi_r  = binsiz_r/2.+bns_r[:-1] # center of each histogram bin

    # calculating g(R) for each time snapshot
    grs = np.zeros((ntimes, len(bns_r)-1))
    for t in range(ntimes):
        grs[t], ndens[t] = g_of_r(bns_r, dists[t], vl.get_vol_i(t))

    grs_avg = np.mean(grs, axis = 0)
    s_t_integrand = grs_avg*np.log(grs_avg) - grs_avg + 1.0
    ent_t = simps(s_t_integrand, radi_r)
    return -0.5 * ent_t * KB_CAL_MOL * np.mean(ndens)

def orien_entropy(dists, angles, vl):
    '''Compute the translational entropy from a list of pairwise distances
    '''
    ntimes, npairs = dists.shape[0], dists.shape[1]
    binsiz_ra, binsiz_a, grs, nden = 0.53, 0.17, [], np.zeros(ntimes)
    bns_ra = np.arange(2.5,20.,binsiz_ra); bns_a = np.arange(0,np.pi,binsiz_a)
    radi_a  = binsiz_a/2.+bns_a[:-1] # center of each histogram bin
    radi_ra = binsiz_ra/2.+bns_ra[:-1] # center of each histogram bin
    nb_a, nb_ra = len(radi_a), len(radi_ra) # Number of bins for each hist
    f = plt.figure(1, figsize = (3.0, 3.0))
    ax = f.add_subplot(111)

    # calculating g(R) for each time snapshot
    grs = np.zeros((ntimes, nb_ra, NANGLES, nb_a))
    grs_r = np.zeros((ntimes, nb_ra))
    sh_shell = np.zeros((nb_ra, NANGLES))
    print(len(bns_ra), nb_ra)
    for t in range(ntimes):
        r_bins = np.digitize(dists[t], bns_ra[:-1])
        grs_r[t], nden[t] = g_of_r(bns_ra, dists[t], vl.get_vol_i(t))
       #for each r bin, compute the g(angle)
        for bn in range():
            this_r = r_bins == bn
            for an in range(NANGLES):
                grs[t][bn][an], _ = g_of_r(bns_a, angles[t*NANGLES + an][this_r])
            gr_mean = np.mean(grs[:,bn,:,:], axis = 0)
            gr_mean[gr_mean==0] = 1.0
            sh_shell[bn] = simps(gr_mean*np.log(gr_mean), radi_a)
            ax.plot(radi_a, grs[t][bn][0])
    plt.show()

    grs_avg = np.mean(grs_r, axis = 0)
    s_o_integrand = grs_avg*np.prod(sh_shell,axis = 1)
    ent_o = simps(s_o_integrand, radi_ra)
    print(ent_o)
    return -0.5 * ent_o * KB_CAL_MOL * np.mean(nden)

def main():
    ''' Given a csv file with 5 angles, oxygen pair distance and dist from wall
        calculate translational and rotational entropy
    '''
    angname=sys.argv[1]; sep=sys.argv[2]; ln=sys.argv[3]; itr=sys.argv[4]
    nm = str(sep)+"_"+str(ln)+"_"+str(itr)
    angC = CSVFile(angname)
    volC = VolFile("run"+nm+".vol")
    
    dis_loc = angC.find_keyword("dis")
    other_loc = angC.find_not_keyword("dis")
    ent_t  = trans_entropy(angC.dat[dis_loc], volC)
    print("This is the translational entropy: {0:.4f}".format(ent_t))
    ent_or = orien_entropy(angC.dat[dis_loc], angC.dat[other_loc], volC)

if __name__=="__main__":
    main()

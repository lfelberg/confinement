import sys
import re
import math
import numpy as np
import itertools
import matplotlib.pyplot as plt 

from scipy.integrate import simps

from csvfile import CSVFile
from util    import d_pbc

def g_of_r(bns, dist, ndensity = 1.0):
    ''' Given an array of distances, histogram to compute g(r) '''
    hist, bins = np.histogram(dist, bins = bns)
    upp, low = bns[1:], bns[:-1]
    gr = hist/(4.0/3.0*np.pi*(np.power(upp, 3.) - np.power(low,3.)))/ndensity
    return gr

def trans_entropy(dists):
    '''Compute the translational entropy from a list of pairwise distances
    '''
    binsiz, grs = 0.03, []
    ent_t, bns = 0.0, np.arange( 2.5, 30., binsiz)
    radi = binsiz/2.+bns[:-1] # center of each histogram bin
    ntimes, npairs = dists.shape[0], dists.shape[1]
    f = plt.figure(1, figsize = (3.0, 3.0))
    ax = f.add_subplot(111)

    # calculating g(R) for each time snapshot
    grs = np.zeros((ntimes, len(bns)-1))
    for t in range(ntimes):
        print(min(dists[t]))
        grs[t] = g_of_r(bns, dists[t])
        ax.plot(radi, grs[-1])
    plt.show()

    grs_avg = np.mean(grs, axis = 0)
    s_t_integrand = grs_avg*np.log(grs_avg) - grs_avg + 1.0
    print(grs_avg.shape, s_t_integrand.shape, radi.shape)
    ent_t = simps(s_t_integrand, radi)
    print(ent_t)
    return ent_t

def main():
    ''' Given a csv file with 5 angles, oxygen pair distance and dist from wall
        calculate translational and rotational entropy
    '''
    angname=sys.argv[1]; sep=sys.argv[2]; ln=sys.argv[3]; itr=sys.argv[4]
    nm = str(sep)+"_"+str(ln)+"_"+str(itr)
    angC = CSVFile(angname)
    
    ent_t = trans_entropy(angC.dat[angC.find_keyword("dis")])

if __name__=="__main__":
    main()

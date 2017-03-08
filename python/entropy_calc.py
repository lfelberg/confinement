import sys
import re
import math
import numpy as np
import itertools
import matplotlib.pyplot as plt 

from csvfile import CSVFile
from util    import d_pbc


def trans_entropy(dists):
    '''Compute the translational entropy from a list of pairwise distances
    '''
    binsiz = 0.5
    ent_t, bns = 0.0, np.arange( 0., 30., binsiz)
    radi = binsiz/2.+bns[:-1] # center of each histogram bin
    upp, low = bns[1:], bns[:-1]
    ntimes, npairs = dists.shape[0], dists.shape[1]
    f = plt.figure(1, figsize = (3.0, 3.0))
    ax = f.add_subplot(111)

    for t in range(ntimes):
        hist, bins = np.histogram(dists[t], bins = bns)
        gr = hist/(4.0/3.0*np.pi*(np.power(upp, 3.) - np.power(low,3.)))
        ax.plot(radi, gr)

    print(bins, radi)
    plt.show()

    return ent_t

def main():
    ''' Given a csv file with 5 angles, oxygen pair distance and dist from wall
        calculate translational and rotational entropy
    '''
    angname=sys.argv[1]; sep=sys.argv[2]; ln=sys.argv[3]; itr=sys.argv[4]

    nm = str(sep)+"_"+str(ln)+"_"+str(itr)
    angC = CSVFile(angname)
    
    print(angC.find_keyword("dis"))
    ent_t = trans_entropy(angC.dat[angC.find_keyword("dis")])

if __name__=="__main__":
    main()

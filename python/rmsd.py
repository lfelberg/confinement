import sys
import re
import numpy as np

from xyzfile import XYZFile

NWALL = 2

def rmsd(dists):
    '''Method to compute the rmsd for one wall'''
    rms = np.zeros(len(dists.atom))
    for i in range(len(dists.atom)):
        rms[i] = np.sum(dists.atom[i,:,-1]*dists.atom[i,:,-1])
    print(len(rms))
    return rms

def print_rms(fname, time, rms):
    '''Print rmsd of carbon walls'''
    f = open(fname, 'w')
    f.write("atime,rmsd\n")
    for i in range(len(rms)):
        f.write("{0},{1:.3f}\n".format(time[i], rms[i]))
    f.close()

def main():
    ''' usage: python rmsd.py graphsep_prefix sep len iter'''
    ntime = 0
    disname=sys.argv[1]; sep=sys.argv[2]; ln=sys.argv[3]; itr=sys.argv[4]
    nm = "run"+str(sep)+"_"+str(ln)+"_"+str(itr)+"_graph"
    for wl in range(NWALL): 
        disC = XYZFile(nm+str(wl)+".dat")
        if ntime == 0: 
            ntime = len(disC.atom)
            msd = np.zeros(ntime)
        msd += rmsd(disC)

    print_rms(nm+".rmsd", disC.time, msd)

if __name__=="__main__":
    main()

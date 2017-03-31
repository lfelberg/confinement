import sys
import numpy as np
import matplotlib.pyplot as plt

from xyzfile import XYZFile
from volfile import VolFile


def get_water_dev(xyz, volC):
    bnz = 5
    oi,hi = xyz.get_inner_wat(); oou,hou = xyz.get_outer_wat() # outside walls
    print(xyz.atom.shape)
    bns = np.arange(-3.,3.5, 0.1)
    hs_t = np.zeros(len(bns)-1)
    for i in range(1,len(xyz.atom)): # for each time snapshot, except first
        x_mn = np.mean(xyz.atom[i,oi,0])
        hist, bins = np.histogram(xyz.atom[i,oi,0]-x_mn, bins = bns)
        hs_t += hist
        

    f = plt.figure(1, figsize = (3.0, 3.0))
    ax, ct, leg = f.add_subplot(111), 0, []
    ax.plot(bins[:-1]+0.05, hs_t)
    plt.show()

    return bns

def print_angles(angls, fname):
    '''Print file of data in csv-like format, angls is a list of values:
       angls = [ [thet1],[thet2], [chi1], [chi2], [phi], [rs] ]'''
    f = open(fname, 'w'); 
    nsnap, stn = len(angls[0]), ''
    vals = ['the1_', 'the2_','chi1_', 'chi2_','phi_', 'dis_', 'walld_']
    for i in range(nsnap):  # writing header
        for j in range(len(vals)):
            stn += "{0}{1},".format(vals[j],i)
    f.write(stn[:-1]+'\n')

    for k in range(len(angls[0][0])):
        st = ''
        for j in range(nsnap):
            for i in range(len(vals)):
                st += "{0:.5f},".format(angls[i][j][k])
        f.write("{0}\n".format(st[:-1]))
    f.close()

def main():
    ''' Given an xyz file, sort all waters, calculate 5 angles between each
        pair on each side of the wall '''
    xyzname=sys.argv[1]; sep=sys.argv[2]; ln=sys.argv[3]; itr=sys.argv[4]

    nm = str(sep)+"_"+str(ln)+"_"+str(itr)
    volC = VolFile("run"+nm+".vol") 
    xyz_cl = XYZFile(xyzname, volC)

    angs = get_water_dev(xyz_cl, volC)
   #print_angles(angs, "run"+nm+"_angles.csv")

if __name__=="__main__":
    main()
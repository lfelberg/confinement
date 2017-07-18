import sys
import numpy as np

from util import d_pbc
from xyzfile import XYZFile
from volfile import VolFile

WOXY = 1
CUT = 8.50

def get_angles_wats(xyz, volC):
    '''Method to get angles between y-z plane and benzene plane'''
    oo = xyz.get_type_i(WOXY); 

    stats = np.zeros((2,len(xyz.atom),xyz.nsol))
    for i in range(1,len(xyz.atom)): # for each time snapshot, except first
        rng = np.array([volC.get_x_rng_i(i), volC.get_y_rng_i(i),
                        volC.get_z_rng_i(i)]) # pbc range

        coms, angs = xyz.get_sol_crd_i(i, rng)
        if angs != []: stats[0,i] = angs
        for j in range(xyz.nsol):
            dst = d_pbc(coms[j][np.newaxis,:],xyz.atom[i,oo], rng)
            stats[1,i] = sum((dst < CUT).astype(int))
    return stats[:,1:]

def print_stats(angls, fname):
    '''Print file of data in csv-like format, angls is a list of values:
       angls = [ [sol1_ang], [sol1_wneig], [sol2_ang], [sol2_wneig]... ]'''
    f = open(fname, 'w'); 
    nsnap, stn = len(angls[0]), ''
    vals = ['ang','neigh']
    for sl in range(len(angls[0][0])):
        for j in range(len(vals)): stn += "sol{0}_{1},".format(sl,vals[j])
    f.write("time,"+stn[:-1]+'\n')

    for j in range(nsnap):
        st = '{0},'.format(j)
        for k in range(len(angls[0][j])): # the number of solutes
            for i in range(len(vals)):
                st += "{0:.5f},".format(angls[i][j][k])
        f.write("{0}\n".format(st[:-1]))
    f.close()

def main():
    xyzname=sys.argv[1]; sep=sys.argv[2]; ln=sys.argv[3]; itr=sys.argv[4]
    sol_nm = sys.argv[5]; sol_ct = sys.argv[6]
    nm = str(sep)+"_"+str(ln)+"_"+itr; volC = VolFile("run"+nm+".vol") 
    xyz_cl = XYZFile(xyzname, volC)
    xyz_cl.sol_ty = sol_nm
    xyz_cl.nsol = int(sol_ct)

    stats = get_angles_wats(xyz_cl, volC)
    print_stats(stats, "run"+nm+"_benz_stats.csv")

if __name__=="__main__":
    main()

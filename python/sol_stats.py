import sys
import numpy as np

from util import d_pbc
from xyzfile import XYZFile
from volfile import VolFile

WOXY = 1
CUT = 6.20

def get_angles_wats(xyz, volC):
    '''Method to get angles between y-z plane and benzene plane'''
    oo = xyz.get_type_i(WOXY); ct = 0

    stats = np.zeros((3,len(xyz.atom),xyz.nsol))
    sol = np.zeros((1,1)); w1 = np.zeros((1,1))
    for i in range(1,len(xyz.atom)): # for each time snapshot, except first
        w1[0,0] = np.mean(xyz.atom[i, xyz.get_graph_wall(0),0])
        w_s = np.std(xyz.atom[i, xyz.get_graph_wall(0),0])
        if w_s > 0.35: continue
        rng = np.array([volC.get_x_rng_i(i), volC.get_y_rng_i(i),
                        volC.get_z_rng_i(i)]) # pbc range

        coms, angs = xyz.get_sol_crd_i(i, rng)

        # computing distance from wall to each center   
        stats[2,ct] = d_pbc(coms[:,0][:,np.newaxis], w1, rng[0], [1.0])

        # saving angle distribtion
        if angs != []: stats[0,ct] = angs
        for j in range(xyz.nsol):
            dst = d_pbc(coms[j][np.newaxis,:],xyz.atom[i,oo], rng)
            sol = coms[j][np.newaxis,0]
            stats[1,ct,j] = sum((dst < CUT).astype(int))
        ct += 1
    return stats[:,:ct]

def print_stats(angls, fname):
    '''Print file of data in csv-like format, angls is a list of values:
       angls = [ [sol1_ang], [sol1_wneig], [sol2_ang], [sol2_wneig]... ]'''
    f = open(fname, 'w'); 
    nsnap, stn = len(angls[0]), ''
    vals = ['ang','neigh', 'dist']
    for sl in range(len(angls[0][0])):
        for j in range(len(vals)): stn += "sol{0}_{1},".format(sl,vals[j])
    f.write("time,"+stn[:-1]+'\n')

    print(angls.shape)
    for j in range(nsnap):
        st = '{0},'.format(j)
        for k in range(len(angls[0][j])): # the number of solutes
            for i in range(len(vals)):
                st += "{0:.2f},".format(angls[i][j][k])
        f.write("{0}\n".format(st[:-1]))
    f.close()

def main():
    xyzname=sys.argv[1]; sep=sys.argv[2]; ln=sys.argv[3]; itr=sys.argv[4]
    sol_nm = sys.argv[5]; sol_ct = sys.argv[6]
    nm = str(sep)+"_"+str(ln)+"_"+itr; volC = VolFile("run"+nm+".vol") 
    nm = str(sep)+"_"+str(ln)+"_"+sol_ct+"_"+itr
    xyz_cl = XYZFile(xyzname, volC)
    xyz_cl.sol_ty = sol_nm
    xyz_cl.nsol = int(sol_ct)

    stats = get_angles_wats(xyz_cl, volC)
    print_stats(stats, "run"+nm+"_sol_stats.csv")

if __name__=="__main__":
    main()

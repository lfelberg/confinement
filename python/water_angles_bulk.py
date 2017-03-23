import sys
import numpy as np
import itertools

from water_angles_util import translate_pbc,unit_vector,angle_between,cal_ang
from water_angles_util import plane_eq,calc_dip,project_plane,find_closest
from xyzfile import XYZFile
from volfile import VolFile

WOXY = 1; WHYD = 2; GRAPHENE = 3

def get_angles(xyz, volC):
    '''Method to get various angles between two waters'''
    # find water oxys, hyds
    oo = xyz.get_type_i(WOXY); hh = xyz.get_type_i(WHYD); 
    t1s, t2s, c1s, c2s, phs, rs, vls, wat_angles = [],[],[],[],[],[],[],[]

    n_w = sum(oo.astype(int)); wat = np.zeros((3, n_w, 3));
    h=np.zeros((2,n_w), dtype=int); h_ct, hidx = 0, 0
    for i in range(len(oo)): #looping over all atoms, make a list of H1 and H2
        if hh[i] == True:
            if h_ct == 0:
                h[0][hidx] = i; h_ct = 1
            elif h_ct == 1:
                h[1][hidx] = i; h_ct = 0; hidx +=1

    for i in range(len(xyz.atom)): # for each time snapshot, except first
        rng = np.array([volC.get_x_rng_i(i), volC.get_y_rng_i(i),
                        volC.get_z_rng_i(i)]) # pbc range

        wat[0] = xyz.atom[i,oo,:]; wat[1] = xyz.atom[i,h[0],:]
        wat[2] = xyz.atom[i,h[1],:]
        t1, t2, c1, c2, ph, r, _ = cal_ang(wat, rng)
        t1s += [t1];t2s += [t2];c1s += [c1];c2s += [c2];phs += [ph];
        rs += [r]; vls += [[np.prod(rng) for i in range(len(r))]]
    return list([t1s, t2s, c1s, c2s, phs, rs, vls]) 

def print_angles(angls, fname):
    '''Print file of data in csv-like format, angls is a list of values:
       angls = [ [thet1],[thet2], [chi1], [chi2], [phi], [rs] ]'''
    f = open(fname, 'w'); 
    nsnap, stn = len(angls[0]), ''
    vals = ['the1_', 'the2_','chi1_', 'chi2_','phi_', 'dis_','vol_']
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

    angs = get_angles(xyz_cl, volC)
    print_angles(angs, "run"+nm+"_angles.csv")

if __name__=="__main__":
    main()

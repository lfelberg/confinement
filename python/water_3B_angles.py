import sys
import numpy as np

from water_3B_angles_util import cal_ang
from xyzfile import XYZFile
from volfile import VolFile

WOXY = 1

def get_angles(xyz, volC):
    '''Method to get various angles between three waters, in confined space'''
    # find water oxys, hyds
    oo = xyz.get_type_i(WOXY); bnz = 5
    oi,_ = xyz.get_inner_wat(); oou,_ = xyz.get_outer_wat() # outside walls
    d1s, d2s, t1s, = [],[],[]

    n_w_in = sum(oi.astype(int)); n_w_ou = sum(oou.astype(int))
    in_wat = np.zeros((1, n_w_in, 3)); ou_wat = np.zeros((1, n_w_ou, 3))

    for i in range(1,3): #len(xyz.atom)): # for each time snapshot, except first
        di1, di2, th1 = [],[],[]
        rng = volC.get_rng_i(i) # pbc range
        in_wat = xyz.atom[i,oi,:]; ou_wat = xyz.atom[i,oou,:];

        d1,d2,t1=cal_ang(in_wat,rng)
        di1+=d1;di2+=d2;th1+=t1;
        d1,d2,t1=cal_ang(in_wat,rng)
        di1+=d1;di2+=d2;th1+=t1;

        d1s += [di1];d2s += [di2];t1s += [th1];
    return list([d1s, d2s, t1s]) 

def print_angles(angls, fname):
    '''Print file of data in csv-like format, angls is a list of values:
       angls = [ [d1],[d2],[thet1]]'''
    f = open(fname, 'w'); 
    nsnap, stn = len(angls[0]), ''
    vals = ['dis1_','dis2_','the1_']
    for i in range(nsnap):
        for j in range(len(vals)):
            stn += "{0}{1},".format(vals[j],i)
    f.write(stn[:-1]+'\n')

    for k in range(len(angls[0][0])): # the number of pairs will change
        st = ""
        for j in range(nsnap):
            for i in range(len(vals)-1):
                st += "{0:.5f},".format(angls[i][j][k])
        f.write("{0}\n".format(st[:-1]))
    f.close()

def main():
    ''' Given an xyz file, sort all waters, calculate 5 angles between each
        pair on each side of the wall '''
    xyzname=sys.argv[1]; sep=sys.argv[2]; ln=sys.argv[3]; itr=sys.argv[4]

    nm = str(sep)+"_"+str(ln)+"_"+str(itr)
    volC = VolFile("run"+nm+".vol"); xyz_cl = XYZFile(xyzname, volC)

    angs = get_angles(xyz_cl, volC)
    print_angles(angs, "run"+nm+"_3B_angles.csv")

if __name__=="__main__":
    main()
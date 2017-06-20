import sys
import numpy as np

from water_angles_util import cal_ang_3B
from util import translate_pbc
from xyzfile import XYZFile
from volfile import VolFile

WOXY = 1

def get_angles(xyz, volC):
    '''Method to get various angles between three waters, in confined space'''
    # find water oxys, hyds
    oo = xyz.get_type_i(WOXY); bnz = 5
    oi,_ = xyz.get_inner_wat(); oou,_ = xyz.get_outer_wat() # outside walls
    d1s, d2s, t1s, = [],[],[]

    for i in range(1,len(xyz.atom)): # for each time snapshot, except first
        di1, di2, th1 = [],[],[]
        rng = volC.get_rng_i(i) # pbc range
        in_wat = xyz.atom[i,oi,:]; ou_wat = xyz.atom[i,oou,:];
        in_wat[1:,0] = translate_pbc(in_wat[0,0], in_wat[1:,0], rng[0])
        ou_wat[1:,0] = translate_pbc(ou_wat[0,0], ou_wat[1:,0], rng[0])

        nm = xyz.get_spacing_for_interlayer()
        b_in = np.digitize(in_wat[:,0], min(in_wat[:,0])+nm)
        b_ou = np.digitize(ou_wat[:,0], min(ou_wat[:,0])+nm)
        
        for j in range(len(nm)):
            if sum((b_in==j).astype(int))>10:
               d1,d2,t1=cal_ang_3B(in_wat[b_in == j],rng)
               di1+=d1;di2+=d2;th1+=t1;
            if sum((b_ou==j).astype(int))>10:
               d1,d2,t1=cal_ang_3B(ou_wat[b_ou == j],rng)
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

    for k in range(len(angls[0][0])):
        st = ""
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
    volC = VolFile("run"+nm+".vol"); xyz_cl = XYZFile(xyzname, volC)

    angs = get_angles(xyz_cl, volC)
    print_angles(angs, "run"+nm+"_3B_angles_layers.csv")

if __name__=="__main__":
    main()

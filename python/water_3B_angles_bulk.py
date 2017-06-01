import sys
import numpy as np

from water_angles_util import cal_ang_3B
from xyzfile import XYZFile
from volfile import VolFile

WOXY = 1

def get_angles(xyz, volC):
    '''Method to get various angles between three waters, in confined space'''
    # find water oxys, hyds
    oo = xyz.get_type_i(WOXY); bnz = 5
    d1s, d2s, t1s, = [],[],[]
    ocrd = xyz.atom[:,oo]

    for i in range(1,len(xyz.atom)): # for each time snapshot, except first
        di1, di2, th1, nm = [],[],[], 20
        rng = volC.get_rng_i(i) # pbc range

        w_bn = np.linspace(min(ocrd[i,:,0]),max(ocrd[i,:,0]), num=nm)
        b_w = np.digitize(ocrd[i,:,0], w_bn)
        
        for j in range(4,len(w_bn)-3):
            if sum((b_w==j).astype(int))>10:
               d1,d2,t1=cal_ang_3B(ocrd[i,b_w == j],rng)
               di1+=d1;di2+=d2;th1+=t1;
       #d1,d2,t1=cal_ang_3B(xyz.atom[i,oo,:],rng)
       #di1+=d1;di2+=d2;th1+=t1;

       #d1s += [di1];d2s += [di2];t1s += [th1];
        d1s += di1;d2s += di2;t1s += th1;
    print(len(d1s), len(t1s))
    return list([d1s, d2s, t1s]) 

def print_angles(angls, fname):
    '''Print file of data in csv-like format, angls is a list of values:
       angls = [ [d1],[d2],[thet1]]'''
    f = open(fname, 'w'); 
    nsnap, stn = len(angls[0]), ''
    vals = ['dis1','dis2','the1']
    for j in range(len(vals)): stn += "{0},".format(vals[j])
    f.write(stn[:-1]+'\n')

    for j in range(nsnap):
        st = ""
        for i in range(len(vals)):
            st += "{0:.5f},".format(angls[i][j])
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

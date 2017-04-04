import sys
import numpy as np

from water_angles_util import translate_pbc,cal_ang
from xyzfile import XYZFile
from volfile import VolFile

WOXY = 1; WHYD = 2; GRAPHENE = 3

def get_angles(xyz, volC):
    '''Method to get various angles between two waters, in confined space'''
    # find water oxys, hyds
    oo = xyz.get_type_i(WOXY); hh = xyz.get_type_i(WHYD); bnz = 5
    t1s, t2s, c1s, c2s, phs, rs, ws, wat_angles = [],[],[],[],[],[],[],[]

    wat = np.zeros((3, sum(oo.astype(int)), 3))
    hs=np.zeros((2,sum(oo.astype(int))), dtype=int)

    for i in range(1,len(xyz.atom)): # for each time snapshot, except first
        th1, th2, ch1, ch2, phi, rr, wd = [],[],[],[],[],[],[]
        rng = np.array([volC.get_x_rng_i(i), volC.get_y_rng_i(i),
                        volC.get_z_rng_i(i)]) # pbc range

        wat[0] = xyz.atom[i,oo,:]
        for h in range(2): wat[h+1]=xyz.atom[i,hs[h],:]

        #binning waters by distribution of x positions
        x = np.mean(wat[0,:,0])
        x_bn = np.arange(0, x*2.0, 0.7); 
        bn = np.digitize(wat[0,:,0]-x, x_bn)
        for j in range(len(x_bn)):
            if sum((bn == j).astype(int)) > 1: #if there is > 1 water in bin
                b_arr = bn == j
                t1,t2,c1,c2,ph,r,w=cal_ang(wat[:,b_arr],rng)
                w = [rng[1]*rng[2]] * len(r)
                th1+=t1;th2+=t2;ch1+=c1;ch2+=c2;phi+=ph;rr+=r;wd+=w

        t1s += [th1];t2s += [th2];c1s += [ch1];c2s += [ch2];phs += [phi];
        rs += [rr];ws += [wd]
    return list([t1s, t2s, c1s, c2s, phs, rs, ws]) 

def print_angles(angls, fname):
    '''Print file of data in csv-like format, angls is a list of values:
       angls = [ [thet1],[thet2], [chi1], [chi2], [phi], [rs] ]'''
    f = open(fname, 'w'); 
    nsnap, stn = len(angls[0]), ''
    vals = ['atim','the1', 'the2','chi1', 'chi2','phi', 'dis', 'vol']
    for j in range(len(vals)):
        stn += "{0},".format(vals[j])
    f.write(stn[:-1]+'\n')

    for j in range(nsnap):
        for k in range(len(angls[0][j])): # the number of pairs will change
            st = '{0},'.format(j)
            for i in range(len(vals)-1):
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

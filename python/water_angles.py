import sys
import numpy as np

from water_angles_util import translate_pbc,cal_ang
from xyzfile import XYZFile
from volfile import VolFile

WOXY = 1; WHYD = 2; GRAPHENE = 3

def get_angles(xyz, disC, volC):
    '''Method to get various angles between two waters, in confined space'''
    # find water oxys, hyds
    oo = xyz.get_type_i(WOXY); hh = xyz.get_type_i(WHYD); bnz = 5
    oi,hi = xyz.get_inner_wat(); oou,hou = xyz.get_outer_wat() # outside walls
    t1s, t2s, c1s, c2s, phs, rs, ws, wat_angles = [],[],[],[],[],[],[],[]

    n_w_in = sum(oi.astype(int)); n_w_ou = sum(oou.astype(int))
    in_wat = np.zeros((3, n_w_in, 3)); ou_wat = np.zeros((3, n_w_ou, 3))
    hs_i=np.zeros((2,n_w_in), dtype=int); hs_o=np.zeros((2,n_w_ou), dtype=int)
    ho_ct, hoidx, hi_ct, hiidx = 0, 0, 0, 0
    for i in range(len(oi)): #looping over all atoms, make a list of H1 and H2
        if hi[i] == True:
            if hi_ct == 0:
                hs_i[0][hiidx] = i; hi_ct = 1
            elif hi_ct == 1:
                hs_i[1][hiidx] = i; hi_ct = 0; hiidx +=1
        if hou[i] == True:
            if ho_ct == 0:
                hs_o[0][hoidx] = i; ho_ct = 1
            elif ho_ct == 1:
                hs_o[1][hoidx] = i; ho_ct = 0; hoidx +=1

    rg=7.0; x_bn = np.arange(-rg, rg, (2.*rg)/9); 
    for i in range(1,len(xyz.atom)): # for each time snapshot, except first
        th1, th2, ch1, ch2, phi, rr, wd = [],[],[],[],[],[],[]
        rng = np.array([volC.get_x_rng_i(i), volC.get_y_rng_i(i),
                        volC.get_z_rng_i(i)]) # pbc range

        in_wat[0] = xyz.atom[i,oi,:]; ou_wat[0] = xyz.atom[i,oou,:];
        for h in range(2):
            in_wat[h+1]=xyz.atom[i,hs_i[h],:];ou_wat[h+1]=xyz.atom[i,hs_o[h],:]

        #binning waters by distribution of x positions
        #For the outer waters, need to wrap with PBCS before computing
        # to make sure the mean is not in the middle of the box
        ou_sft = np.zeros((ou_wat.shape[1:])); ou_sft[0] = ou_wat[0,0]
        ou_sft[1:] = translate_pbc(ou_wat[0,0], ou_wat[0,1:], rng)
        x_in = np.mean(in_wat[0,:,0]);x_ou = np.mean(ou_sft[:,0]);
        b_in = np.digitize(in_wat[0,:,0]-x_in, x_bn)
        b_ou = np.digitize(ou_sft[:,0]-x_ou, x_bn)
        for j in range(len(x_bn)):
            if sum((b_in == j).astype(int)) > 1: #if there is > 1 water in bin
                b_arr = b_in == j
                t1,t2,c1,c2,ph,r,w=cal_ang(in_wat[:,b_arr],rng)
                w = [rng[1]*rng[2]] * len(r)
                th1+=t1;th2+=t2;ch1+=c1;ch2+=c2;phi+=ph;rr+=r;wd+=w
            if sum((b_ou == j).astype(int)) > 1:
                b_arr = b_ou == j
                t1,t2,c1,c2,ph,r,w=cal_ang(ou_wat[:,b_arr],rng)
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
    disC = XYZFile("run"+nm+".dist", VolFile("")) 
    xyz_cl = XYZFile(xyzname, volC)

    angs = get_angles(xyz_cl, disC, volC)
    print_angles(angs, "run"+nm+"_angles.csv")

if __name__=="__main__":
    main()

import sys
import re
import math
import numpy as np

from xyzfile import XYZFile
from volfile import VolFile
from util    import d_pbc, translate_pbc

dr = 0.05
LMAX = 15.0+dr
HIS = np.arange(0, LMAX, dr);

def g_of_r(dists, rang):
    '''Given an array of distances, histogram and calc g(r)'''
    # Normalization factors for g(r)
    gr, vl, dim = [], 1.0, len(rang)
    rd_rng = np.arange(0, LMAX, dr); low, upp = rd_rng[:-1], rd_rng[1:]
    bin_gr = len(rd_rng)
    if dim == 3: rd_mu = 4.0/3.0*np.pi*(np.power(upp,3.) - np.power(low,3.))
    else: rd_mu = 2.0*np.pi*(np.power(upp,2.) - np.power(low,2.)) 

    # number density normalization
    for d in range(dim): vl *= rang[d]
    his_den, benz = np.histogram(dists, bins=rd_rng)
    return his_den*vl*2./rd_mu, len(dists) #DIVIDE BY TWO FOR NORM SCHEME

def gr_cal(crd1,crd0,rang,same):
    '''Given an array of molecules you would like to compute the distance 
       on distance of crd0 from graphene. Will do 2D or 3D, depending on the
       length of the rang array.'''
    dist,dim=[],len(rang); st=3-dim; pbc=[1.,1.] if dim==2 else [1.,1.,1.]
    for fm in range(len(crd0)): # for each of type 0, find dist to each type 1
        to_ar = crd1[:,st:] if same == False else np.delete(crd0[:,st:],fm,0)
        frm_ar = np.repeat(crd0[np.newaxis,fm,st:],len(to_ar),axis=0)
        dist.append(list(d_pbc(to_ar, frm_ar, rang, pbc)))

    his, ct = g_of_r(np.array([it for subl in dist for it in subl]), rang)
    return his, ct
        
def get_gr(xyz, volC, grPair):
    '''Method to get the g(r) for two atom types'''
    # find atoms of type 0 and type 1
    ty0 = xyz.get_type_i(grPair[0]); ty1 = xyz.get_type_i(grPair[1]); bnz = 2
    g_r_3, g_r_2_m, rng_m = [], [], np.zeros((3))

    # depending on the system, there will be varying # interlayers
    if ("_12_" in xyz.xyzfname): nm = 3
    elif ("_13_" in xyz.xyzfname or "_14_" in xyz.xyzfname): 
         nm = np.array([0,0.35,1.24,3.24,3.90,5.9,6.7,10.0])
    elif ("_16_" in xyz.xyzfname):#nm = 7
         nm = np.array([0,0.55,1.24,3.24,3.70,5.4,6.1,8.1,8.7,14.0])
    xb_ct=np.zeros((nm,2)); xb_his=np.zeros((nm,2,len(HIS)-1))

    # for the g(r) pairs, need to get groups on the same wall side
    ty00 = np.all(np.array([ty0, xyz.get_inner_ats()]), axis=0) # btw 2 walls
    ty10 = np.all(np.array([ty1, xyz.get_inner_ats()]), axis=0)

    ty01 = np.all(np.array([ty0, xyz.get_outer_ats()]), axis=0) # out 2 walls
    ty11 = np.all(np.array([ty1, xyz.get_outer_ats()]), axis=0)

    for i in range(1,len(xyz.atom)):
        rng = np.array([volC.get_x_rng_i(i), volC.get_y_rng_i(i),
                        volC.get_z_rng_i(i)]) # pbc range
        rng_m += rng # take average of range for printing

        # get cords for each type, and wrap x coords 
        cd00 = xyz.atom[i,ty00]; cd01 = xyz.atom[i,ty01]
        cd10 = xyz.atom[i,ty10]; cd11 = xyz.atom[i,ty11]

        c00[1:,0] = translate_pbc(c00[0,0],c00[1:,0],rng[0])
        c01[1:,0] = translate_pbc(c01[0,0],c01[1:,0],rng[0])
        c10[:,0] = translate_pbc(c00[0,0],c10[:,0],rng[0])
        c11[:,0] = translate_pbc(c01[0,0],c11[:,0],rng[0])

        in_bn = min(c00[:,0])+nm
        ou_bn = min(c10[:,0])+nm
        b00 = np.digitize(c00[:,0], in_bn)
        b01 = np.digitize(c01[:,0], ou_bn)
        if grPair[0] == grPair[1]: 
            c10=np.copy(c00); c11=np.copy(c01); b10=b00; b11=b01; sm=True
        else: 
            sm = False
            b10 = np.digitize(c10[:,0], in_bn)
            b11 = np.digitize(c11[:,0], ou_bn)
        
        for j in range(1,len(in_bn)):
            if sum((b00==j).astype(int))>50 and sum((b10==j).astype(int))>50:
               gr2, pr_ct2 = gr_cal(c00[b00==j],c10[b10==j], rng[1:], sm)
               xb_his[j,1,:] += gr2;xb_ct[j,1] += pr_ct2;
               
            if sum((b01==j).astype(int))>50 and sum((b11==j).astype(int))>50:
               gr2, pr_ct2 = gr_cal(c01[b01==j],c11[b11==j], rng[1:], sm)
               xb_his[j,1,:] += gr2;xb_ct[j,1] += pr_ct2;


    # Time averages of histograms
    xb_ct[xb_ct == 0.] = 1.0
    g_r_3_m = xb_his[:,0,:]/xb_ct[:,0, np.newaxis].astype(int);
    g_r_2_m = xb_his[:,1,:]/xb_ct[:,1, np.newaxis].astype(int);
    return g_r_3_m, g_r_2_m, list(range(len(g_r_2_m)))

def print_gr(x, grs, fname):
    '''Print distances to carbon wall in xyz like format'''
    print(grs.shape)
    f = open(fname, 'w'); 
    f.write("Bin,"); st = ""
    dimbins = HIS[:-1]+dr/2.
    for i in range(len(x)): st += "x{0:.3f},".format(x[i])
    f.write("{0}\n".format(st[:-1]))
    for i in range(len(dimbins)-1):
        st = ""
        for j in range(len(grs[i])): st += "{0:.5f},".format(grs[i][j])
        f.write("{0:.4f},{1}\n".format(dimbins[i], st[:-1]))
    f.close()

def main():
    ''' Given a list of pairs for g(r) calcs, do this as a function as 
        distance from wall, but only 2D and only for your wall side '''
    xyzname=sys.argv[1]; sep=sys.argv[2]; ln=sys.argv[3]; itr=sys.argv[4]
    n_gr = int(sys.argv[5]); grS, grPr = 6, []
    for i in range(grS, grS + n_gr*2, 2): 
        grPr.append([int(sys.argv[i]),int(sys.argv[i+1])])

    nm = str(sep)+"_"+str(ln)+"_"+str(itr)
    volC = VolFile("run"+nm+".vol") 
    xyz_cl = XYZFile(xyzname, volC)
    print("Arry shape", xyz_cl.atom.shape)

    for i in range(n_gr):
        prNm = '_'+str(grPr[i][0])+'_'+str(grPr[i][1])+'defined.csv'
        grs3D, grs2D, xrng = get_gr(xyz_cl, volC, grPr[i])
        print_gr(xrng, grs2D.T, 'g_r_2D_'+nm+prNm)

if __name__=="__main__":
    main()

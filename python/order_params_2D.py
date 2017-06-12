import sys
import numpy as np
import scipy 
from scipy.special import *

from water_angles_util import find_in_cutoff, project_plane,angle_between
from util import translate_pbc, cart_to_sph
from xyzfile import XYZFile
from volfile import VolFile

WOXY = 1

def theta_cal(oxy, nncrds):
    '''Given a list of vectors (from 1 atom to ct nearest neighbors)
       calculate the angle between each vector and Zaxis'''
    x_pos = np.mean(nncrds[:,0]); o_o_vecs = nncrds - oxy
    oo_v_proj = project_plane(o_o_vecs, np.array((1,0,0))[np.newaxis,:])
    # creating vector that is z axis
    ax_proj = np.array([0.0,0.0,1.0+oxy[0,2]])[np.newaxis,:]
    return angle_between(ax_proj, oo_v_proj)

def phi_n_cal(oxy, nns, n):
    '''Given a list of vectors in cartesian coords, calculate the SH for 
       a range of l, -l <= m <= l'''
    thes = theta_cal(oxy, nns); phi_n = np.sum(np.exp(float(n)*1j*thes))
    nb = float(len(nns))
    return phi_n * (1./nb)   

def cal_op(oxygens, rng, cutoff, n):
    '''Given a list of oxygens, find the nearest neighbors for each in the 
       list ( r < 3.7 ) and compute phi_n OP for each '''
    nn = [] # nn is a list of all nearest neighs for each oxygen
    phin = np.zeros(len(oxygens),dtype=complex)

    for ox in range(len(oxygens)):
        curr = oxygens[ox][np.newaxis,:]; others = np.delete(oxygens,ox,0)
        ot_wr = translate_pbc(curr, others, rng) # trans other ox
        clo = find_in_cutoff(curr, ot_wr, cutoff)[0]
        nn_ox = np.where(clo == True)[0]

        if len(nn_ox) == 0: 
            nn += [[]]; continue
        phin[ox] = phi_n_cal(curr, ot_wr[clo], n); nn += [nn_ox]
    print("mean: {0}".format(np.mean(phin.real)))
    return phin

def get_2D_order(xyz, volC, co = 3.7, n = 6):
    '''Method to get various angles between three waters, in confined space'''
    # find water oxys
    oo = xyz.get_type_i(WOXY);phi_all, phi_b_all = [],[]
    oi,_ = xyz.get_inner_wat(); oou,_ = xyz.get_outer_wat() # outside walls

    # depending on the system, there will be varying # interlayers
    nm = xyz.get_spacing_for_interlayer()

    for i in range(1,len(xyz.atom)): # for each time snapshot, except first
        rng = volC.get_rng_i(i) # pbc range
        in_wat = xyz.atom[i,oi,:]; ou_wat = xyz.atom[i,oou,:];
        in_wat[1:,0] = translate_pbc(in_wat[0,0], in_wat[1:,0], rng[0])
        ou_wat[1:,0] = translate_pbc(ou_wat[0,0], ou_wat[1:,0], rng[0])

        in_bn = min(in_wat[:,0])+nm; b_in = np.digitize(in_wat[:,0], in_bn)
        ou_bn = min(ou_wat[:,0])+nm; b_ou = np.digitize(ou_wat[:,0], ou_bn)
        
        for j in range(1,len(in_bn)):
            phi = cal_op(in_wat[b_in == j],rng,co,n); phi_all += list(phi)
            phi = cal_op(ou_wat[b_ou == j],rng,co,n); phi_all += list(phi)

    return list([phi_all]) #, phi_b_all]) 

def print_ops(op, fname):
    '''Print file of data in csv-like format, ops:
       op = [phi, phi_b]'''
    f = open(fname, 'w'); vals = ['phi']
    stn = vals[0]; f.write(stn+'\n')

    for k in range(len(op[0])):
        st = ""
        for j in range(len(op)):
            st += "{0:.7f},".format(op[j][k].real)
        f.write("{0}\n".format(st[:-1]))
    f.close()

def main():
    ''' Given an xyz file, sort all waters, calculate 5 angles between each
        pair on each side of the wall '''
    xyzname=sys.argv[1]; sep=sys.argv[2]; ln=sys.argv[3]; itr=sys.argv[4]
    cutoff, n = float(sys.argv[5]), int(sys.argv[6])

    nm = sep+"_"+ln+"_"+itr
    volC = VolFile("run"+nm+".vol"); xyz_cl = XYZFile(xyzname, volC)

    phi = get_2D_order(xyz_cl, volC, cutoff, n)
    print_ops(phi, "run"+nm+"_op_2D.csv")

if __name__=="__main__":
    main()

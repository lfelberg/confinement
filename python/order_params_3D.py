import sys
import numpy as np
import scipy 
from scipy.special import *

from water_angles_util import find_in_cutoff
from util import translate_pbc, cart_to_sph
from xyzfile import XYZFile
from volfile import VolFile

WOXY = 1

def sph_harm(r, l, m):
    '''Calculate the spherical harmonics for a given vector r = [x,y,z]'''
    r_sph = cart_to_sph(r)
    theta = r_sph[:,1][np.newaxis,:] # polar [0, pi]
    phi   = r_sph[:,2][np.newaxis,:] # azimuth [0, 2pi]
    m = np.repeat(m[:,np.newaxis],phi.shape[1],axis=1)

    Ylm = scipy.special.sph_harm(m, l, phi, theta)  # SH in python has polar/az reversed
    return Ylm

def q_lm_cal(r_cart, l):
    '''Given a list of vectors in cartesian coords, calculate the SH for 
       a range of l, -l <= m <= l'''
    ylm = np.sum(sph_harm(r_cart, l, np.arange(-l,l+1)),axis=1)
    nb = float(len(r_cart))
    return ylm * (1./nb)   

def q_lm_bar_cal(i, qlm, nn):
    '''Given a list of neareast neighbors and a list of qs, calc qlm bar'''
    qlm_b = np.copy(qlm[i])
    for j in range(len(nn)): 
        qlm_b += qlm[nn[j]]  
    return qlm_b * ( 1. / (float(len(nn))+1.))

def q_l_bar_cal(qlm_bar,l):
    '''Given a list of qlm_bar (one for each oxygen, and len(qlm_b) = 2*l + 1
       Calc ql_bar = sqrt((4pi/2l+1)*sum_{-l<=m<=l} |qlm_bar|^2 '''
    pref = (4.0 * np.pi ) / ( 2.*l + 1. )
    sm = np.sum(qlm_bar.real*qlm_bar.real + qlm_bar.imag*qlm_bar.imag, axis=1)
    return np.sqrt( pref * sm )

def cal_op(oxygens, rng):
    '''Given a list of oxygens, find the nearest neighbors for each in the 
       list ( r < 3.7 ) and compute q OPs for each '''
    nn = [] # nn is a list of all nearest neighs for each oxygen
    q4m = np.zeros((len(oxygens),len(np.arange(-4,5))),dtype=complex)
    q6m = np.zeros((len(oxygens),len(np.arange(-6,7))),dtype=complex)

    for ox in range(len(oxygens)):
        curr = oxygens[ox][np.newaxis,:]; others = np.delete(oxygens,ox,0)
        ot_wr = translate_pbc(curr, others, rng) # trans other ox
        clo = find_in_cutoff(curr, ot_wr, 3.7)[0]
        nn_ox = np.where(clo == True)[0]

        if len(nn_ox) == 0: 
            nn += [[]]
            continue
        q4m[ox] = q_lm_cal(ot_wr[clo] - curr, 4)
        q6m[ox] = q_lm_cal(ot_wr[clo] - curr, 6)
        nn += [nn_ox]

    q4m_b = np.zeros((len(oxygens),len(np.arange(-4,5))),dtype=complex)
    q6m_b = np.zeros((len(oxygens),len(np.arange(-6,7))),dtype=complex)

    for ox in range(len(oxygens)):
        q4m_b[ox] = q_lm_bar_cal(ox,q4m,nn[ox])
        q6m_b[ox] = q_lm_bar_cal(ox,q6m,nn[ox])

   #print(min(q4m_b.real.flatten()),max(q4m_b.real.flatten()),min(q4m_b.imag.flatten()),max(q4m_b.imag.flatten()))
   #print(min(q6m_b.real.flatten()),max(q6m_b.real.flatten()),min(q6m_b.imag.flatten()),max(q6m_b.imag.flatten()))
    q4_bar = q_l_bar_cal(q4m_b,4.0)
    q6_bar = q_l_bar_cal(q6m_b,6.0)
    return q4_bar, q6_bar

def get_order(xyz, volC):
    '''Method to get various angles between three waters, in confined space'''
    # find water oxys, hyds
    oo = xyz.get_type_i(WOXY); bnz = 5
    oi,_ = xyz.get_inner_wat(); oou,_ = xyz.get_outer_wat() # outside walls
    q4_b_all, q6_b_all = [], []

    nm = xyz.get_spacing_for_interlayer()

    for i in range(1,len(xyz.atom)): # for each time snapshot, except first
        rng = volC.get_rng_i(i) # pbc range
        in_wat = xyz.atom[i,oi,:]; ou_wat = xyz.atom[i,oou,:];
        in_wat[1:,0] = translate_pbc(in_wat[0,0], in_wat[1:,0], rng[0])
        ou_wat[1:,0] = translate_pbc(ou_wat[0,0], ou_wat[1:,0], rng[0])

        b_in = np.digitize(in_wat[:,0], min(in_wat[:,0])+nm)
        b_ou = np.digitize(ou_wat[:,0], min(ou_wat[:,0])+nm)
        
        for j in range(1,len(in_bn)):
            q4_b, q6_b = cal_op(in_wat[b_in == j],rng)
            q4_b_all += list(q4_b); q6_b_all += list(q6_b)
            q4_b, q6_b = cal_op(ou_wat[b_ou == j],rng)
            q4_b_all += list(q4_b); q6_b_all += list(q6_b)

    return list([q4_b_all, q6_b_all]) 

def print_ops(ops, fname):
    '''Print file of data in csv-like format, ops:
       ops = [q4_bar, q6_bar]'''
    f = open(fname, 'w'); vals = ['q4_bar','q6_bar']
    stn = vals[0]+","+vals[1]
    f.write(stn+'\n')

    for k in range(len(ops[0])):
        st = ""
        for j in range(len(ops)):
            st += "{0:.7f},".format(ops[j][k])
        f.write("{0}\n".format(st[:-1]))
    f.close()

def main():
    ''' Given an xyz file, sort all waters, calculate 5 angles between each
        pair on each side of the wall '''
    xyzname=sys.argv[1]; sep=sys.argv[2]; ln=sys.argv[3]; itr=sys.argv[4]

    nm = str(sep)+"_"+str(ln)+"_"+str(itr)
    volC = VolFile("run"+nm+".vol"); xyz_cl = XYZFile(xyzname, volC)

    angs = get_order(xyz_cl, volC)
    print_ops(angs, "run"+nm+"_order_params.csv")

if __name__=="__main__":
    main()

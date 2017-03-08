import sys
import re
import math
import numpy as np
import itertools

from xyzfile import XYZFile
from volfile import VolFile
from util    import d_pbc

WOXY = 1; WHYD = 2; GRAPHENE = 3

QOXY = -1.0484
QHYD =  0.5242

def translate_pbc(c1, c2, rng):
    '''Translate c2 to distance sub rng/2 wrt c1
       c1 and c2 are arrays of coordinates, rng is PBC range in x,y,z'''
    boxl = np.round((c1-c2)/rng) # find what is rounded distance
    return c2 + boxl*rng

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.sqrt(np.sum(vector*vector,axis=1)[:,np.newaxis])

def angle_between(v1, v2, v1nrm = False, v2nrm = False):
    """ Returns the angle in radians between vectors 'v1' and 'v2'
        v1 and v2 have the dimensions = [dim, nwat]::
        >>> angle_between(np.array((1,0,0))[np.newaxis,:],
                          np.array((0,1,0))[np.newaxis,:])
        1.5707963267948966
        >>> angle_between(np.array((1, 0, 0))[np.newaxis,:],
                          np.array((1, 0, 0))[np.newaxis,:])
        0.0
        >>> angle_between(np.array((1, 0, 0))[np.newaxis,:],
                          np.array((-1, 0, 0))[np.newaxis,:])
        3.141592653589793
    """
    if v1nrm == True: v1_u = v1   #Normalize if needed
    else:             v1_u = unit_vector(v1)
    if v2nrm == True: v2_u = v2
    else:             v2_u = unit_vector(v2)

    return np.arccos(np.clip(np.sum(v1_u*v2_u, axis = 1), -1.0, 1.0))

def plane_eq(c1, c2, c3):
    '''Given three points, calculate a unit normal to the plane made by those 
       points: c1 = c2 = c3 = [dim, nwat] = cp, v1, v2, return
       Math: norm = (c2-c1) x (c3-c1) / || (c2-c1) x (c3-c1) ||
             where x=cross prod(a,b) = [a[1]b[2]-a[2]b[1], a[2]b[0]-a[0]b[2],
                                        a[0]b[1]-a[1]b[0]] ::
       >>> plane_eq(np.array(1,2,3)[:,np.newaxis],np.array(4,6,9)[:,np.newaxis],
                    np.array(12,11,9)[:,np.newaxis])
       np.array([-0.5076004073, 0.8121606517, -0.2876402308])
       >>> plane_eq(np.array(0,0,0)[:,np.newaxis],np.array(1,0,0)[:,np.newaxis],
                    np.array(0,11,0)[:,np.newaxis])
       np.array([0.0, 0.0, -1.0])
    '''
    v1 = c2 - c1; v2 = c3 - c1
    cp = np.array([v1[1]*v2[2]-v1[2]*v2[1], 
                   -1*(v1[0]*v2[2]-v1[2]*v2[0]), 
                    v1[0]*v2[1]-v1[1]*v2[0]])

    return cp / np.sqrt(np.sum(cp*cp, axis = 0))

def calc_dip(crd):
    '''Given an array of coords, dim = [nwat, dim=3], cal the dipole moment
        mu = sum(vh1_o*(qh)+vh2_o(qh)), where vh1_o = r_hi - r_o'''
    dip = (crd[1]-crd[0] + crd[2]-crd[0]) * QHYD
    return dip # vec frm org = oxy, e.g. dip has oxy position subtracted frm it

def project_plane(vec, norm):
    '''Method to project given vector onto plane described by normal
       vec = [nwat, dim], norm = [nwat, dim]
       Math: proj(u onto n) = u - [ (u * n)/ ||n||^2 ] n 
             where u, n are vectors, * is dot product::
       >>> project_plane(np.array((10.,11.,23.))[np.newaxis,:], 
                         np.array((0,1,0))[np.newaxis,:])
       array([[ 10.,   0.,  23.]])
       >>> project_plane(np.array((10.,11.,23.))[np.newaxis,:], 
                         np.array((1,1,1))[np.newaxis,:])
       array([[-4.66666667, -3.66666667,  8.33333333]])
    '''
    return (vec-(np.sum(vec*norm,axis=1)/
                 np.sum(norm*norm,axis=1))[:,np.newaxis]*norm)

def cal_ang(w_coords, rng):
    '''Given a list of coords (dimensions: [nat = 3, nwat, ndim = 3]), 
       move the coords into the nearest image of water of interest,
       calculate the dipole moments and angles''' 
    t1, t2, c1, c2, ph, r = [], [], [], [], [], []
    # mve hyds to be w/in L/2 of their oxy
    w_coords[1] = translate_pbc(w_coords[0],w_coords[1],rng)
    w_coords[2] = translate_pbc(w_coords[0],w_coords[2],rng)

    for i in range(w_coords.shape[1]-1):
        curr = w_coords[:,i]; 
        others = w_coords[:,i+1:]; ot_wr = np.zeros(others.shape)
        cur_ar = np.repeat(curr[np.newaxis,0,:], others.shape[1], axis = 0) 
        curr = curr[:,np.newaxis,:]
        ot_wr[0] = translate_pbc(cur_ar, others[0], rng) # trans other ox
        ot_wr[1] = translate_pbc(ot_wr[0], others[1], rng) # trans hydrogs
        ot_wr[2] = translate_pbc(ot_wr[0], others[2], rng)

        inter_mol_ax = ot_wr[0] - curr[0]
        mu_cur = calc_dip(curr); mu_oth = calc_dip(ot_wr)
        
        the_1 = angle_between(mu_cur, inter_mol_ax)
        the_2 = angle_between(mu_oth, inter_mol_ax)

        w1_nrm = plane_eq(curr[0,:,:].T,curr[1,:,:].T,curr[2,:,:].T).T
        chi1 = angle_between(w1_nrm, inter_mol_ax, True, False)
        
        w2_nrm = plane_eq(ot_wr[0,:,:].T,ot_wr[1,:,:].T,ot_wr[2,:,:].T).T
        chi2 = angle_between(w2_nrm, inter_mol_ax, True, False)
        
        mu_cur_proj = project_plane(mu_cur, inter_mol_ax)
        mu_oth_proj = project_plane(mu_oth, inter_mol_ax)
        phi = angle_between(mu_cur_proj, mu_oth_proj)
        
        dists = d_pbc(curr[0], ot_wr[0], rng) # cal O-O distance

        t1.append(the_1); t2.append(the_2); c1.append(chi1); c2.append(chi2); 
        ph.append(phi); r.append(dists); 

    t1 = list(itertools.chain(*t1));t2=list(itertools.chain(*t2));
    c1 = list(itertools.chain(*c1));c2=list(itertools.chain(*c2));
    ph = list(itertools.chain(*ph));r =list(itertools.chain(*r));
    return t1, t2, c1, c2, ph, r

def get_angles(xyz, volC):
    '''Method to get various angles between two waters'''
    # find water oxys, hyds
    oo = xyz.get_type_i(WOXY); hh = xyz.get_type_i(WHYD); 
    t1s, t2s, c1s, c2s, phs, rs, wat_angles = [],[],[],[],[],[],[]

    n_w = sum(oo.astype(int)); wat = np.zeros((3, n_w, 3));
    h=np.zeros((2,n_w), dtype=int); h_ct, hidx = 0, 0
    for i in range(len(oo)): #looping over all atoms, make a list of H1 and H2
        if hh[i] == True:
            if h_ct == 0:
                h[0][hidx] = i; h_ct = 1
            elif hi_ct == 1:
                h[1][hidx] = i; h_ct = 0; hidx +=1

    for i in range(1,len(xyz.atom)): # for each time snapshot, except first
        rng = np.array([volC.get_x_rng_i(i), volC.get_y_rng_i(i),
                        volC.get_z_rng_i(i)]) # pbc range

        wat[0] = xyz.atom[i,oo,:]; wat[1] = xyz.atom[i,h[0],:]
        wat[2] = xyz.atom[i,h[1],:]
        t1, t2, c1, c2, ph, r, w = cal_ang(wat, rng)
        t1s += [t1];t2s += [t2];c1s += [c1];c2s += [c2];phs += [ph];
        rs += [r];
    return list([t1s, t2s, c1s, c2s, phs, rs]) 

def print_angles(angls, fname):
    '''Print file of data in csv-like format, angls is a list of values:
       angls = [ [thet1],[thet2], [chi1], [chi2], [phi], [rs] ]'''
    f = open(fname, 'w'); 
    nsnap, stn = len(angls[0]), ''
    vals = ['the1_', 'the2_','chi1_', 'chi2_','phi_', 'dis_']
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

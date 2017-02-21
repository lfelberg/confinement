import sys
import re
import math
import numpy as np
from xyzfile import XYZFile
from volfile import VolFile

WOXY = 1; WHYD = 2; GRAPHENE = 3
BIN_GR = 80; LMAX = 12.0

QOXY = -1.0484
QHYD =  0.5242

def translate_pbc(c1, c2, rng):
    '''Translate c2 to distance sub rng/2 wrt c1
       c1 and c2 are arrays of coordinates, rng is PBC range in x,y,z'''
    boxl = np.round((c1-c2)/rng) # find what is rounded distance
    return c2 + boxl*rng

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::
            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.sum(v1_u*v2_u, axis = 1), -1.0, 1.0))
   #return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)) # orig online


def plane_eq(c1, c2, c3):
    '''Given three points, calculate a normal to the plane made by those 
       points'''
    v1 = c2 - c1; v2 = c3 - c1
    print(c2.shape)
    cp = np.array([v1[1]*v2[2]-v1[2]*v2[1], 
                   -1*(v1[0]*v2[2]-v1[2]*v2[0]), 
                    v1[0]*v2[1]-v1[1]*v2[0]])

    print(v1.shape, c1.shape, cp.shape)

    return - (cp[0]*c1 + cp[1]*c2 + cp[2]*c3)


def calc_dip(crd):
    '''Given an array of coords, dim = [nwat, dim=3], cal the dipole moment
        mu = sum(vh1_o*(qo-qh)+vh2_o(qo-qh))'''
    dif_q = QOXY - QHYD
    dip = (crd[1]-crd[0] + crd[2]-crd[0]) * dif_q
   #if crd.shape[1] < 2:
   #    dip_v = dip+crd[0]
   #    print("draw line {{ {0} {1} {2} }} {{ {3} {4} {5} }}".format(crd[0][0][0],
   #           crd[0][0][1], crd[0][0][2], dip_v[0][0],dip_v[0][1], dip_v[0][2]))
    return dip+crd[0]

def cal_ang(w_coords, rng):
    '''Given a list of coords (dimensions: nat=3, nwat, ndim = 3), 
       move the coords into the nearest image of water of interest,
       calculate the dipole moments and angles''' 
       
    # mve hyds to be w/in L/2 of their oxy
    w_coords[1] = translate_pbc(w_coords[0],w_coords[1],rng)
    w_coords[2] = translate_pbc(w_coords[0],w_coords[2],rng)

   #w_rshp = np.zeros((1,w_coords.shape[1]*3, 3)); typ = np.zeros((w_coords.shape[1]*3),dtype=int)
   #w_rshp[0,0:w_coords.shape[1],:] = w_coords[0]
   #w_rshp[0,w_coords.shape[1]:w_coords.shape[1]*2,:] = w_coords[1]
   #w_rshp[0,w_coords.shape[1]*2:w_coords.shape[1]*3,:] = w_coords[2]
   #typ[w_coords.shape[1]:] = 1
   #volC = VolFile("")
   #xyzC = XYZFile("test.xyz", volC, w_rshp, typ)

    for i in range(2): #w_coords.shape[1]-1):
        curr = w_coords[:,i]; 
        others = w_coords[:,i+1:]; ot_wr = np.zeros(others.shape)
        cur_ar = np.repeat(curr[np.newaxis,0,:], others.shape[1], axis = 0) 
        curr = curr[:,np.newaxis,:]
        ot_wr[0] = translate_pbc(cur_ar, others[0], rng) # trans other ox
        ot_wr[1] = translate_pbc(ot_wr[0], others[1], rng) # trans hydrogs
        ot_wr[2] = translate_pbc(ot_wr[0], others[2], rng)

        inter_mol_ax = ot_wr[0] - curr[0]

        mu_cur = calc_dip(curr)
        mu_oth = calc_dip(ot_wr)
        
        the_1 = angle_between(mu_cur, inter_mol_ax)
        the_2 = angle_between(mu_oth, inter_mol_ax)
        print(inter_mol_ax.shape, mu_cur.shape, the_1.shape, the_2.shape)

        w1_nrm = plane_eq(curr[0,:,:].T,curr[1,:,:].T,curr[2,:,:].T).T
        chi1 = angle_between(w1_nrm, inter_mol_ax)
        
        w2_nrm = plane_eq(ot_wr[0,:,:],ot_wr[1,:,:],ot_wr[2,:,:])
        chi2 = angle_between(w2_nrm, inter_mol_ax)
        
        

def get_angles(xyz, disC, dists, volC):
    '''Method to get various angles between two waters'''
    # find water oxys, hyds
    oo = xyz.get_type_i(WOXY); hh = xyz.get_type_i(WHYD); bnz = 5; rng_m = 0
    oi,hi = xyz.get_inner_wat(); oou,hou = xyz.get_outer_wat() # outside walls
    wat_angles = []

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

    # for x dir, use bins 0 -> x_quarter. ATs will be max L/2 from wall
    x_hlf = volC.get_x_max()/2.0 # Half of box divides 2 walls
    x_binz = np.arange(0.5, x_hlf, x_hlf/float(bnz)); 

    for i in range(1,len(xyz.atom)): # for each time snapshot, except first
        rng = np.array([volC.get_x_rng_i(i), volC.get_y_rng_i(i),
                        volC.get_z_rng_i(i)]) # pbc range
        rng_m += rng # take average of range for printing

        in_wat[0] = xyz.atom[i,oi,:]; in_wat[1] = xyz.atom[i,hs_i[0],:]
        in_wat[2] = xyz.atom[i,hs_i[1],:]
        wat_angles.append(cal_ang(in_wat, rng))

        ou_wat[0] = xyz.atom[i,oou,:]; ou_wat[1] = xyz.atom[i,hs_o[0],:]
        ou_wat[2] = xyz.atom[i,hs_o[1],:]
        

    return rng_m/float(len(xyz.atom))

def print_gr(x, grs, fname):
    '''Print distances to carbon wall in xyz like format'''
    print(grs.shape)
    f = open(fname, 'w'); 
    f.write("Bin,"); st = ""
    dimbins = np.arange(0, LMAX, LMAX/float(BIN_GR))
    for i in range(len(x)): st += "{0:.5f},".format(x[i])
    f.write("{0}\n".format(st[:-1]))
    for i in range(BIN_GR):
        st = ""
        for j in range(len(grs[i])): st += "{0:.5f},".format(grs[i][j])
        f.write("{0:.4f},{1}\n".format(dimbins[i], st[:-1]))
    f.close()

def main():
    ''' Given an xyz file, sort all waters, calculate 5 angles between each
        pair on each side of the wall '''
    xyzname=sys.argv[1]; sep=sys.argv[2]; ln=sys.argv[3]; itr=sys.argv[4]

    nm = str(sep)+"_"+str(ln)+"_"+str(itr)
    volC = VolFile("run"+nm+".vol") 
    disC = XYZFile("run"+nm+".dist", VolFile("")) 
    xyz_cl = XYZFile(xyzname, volC)

    print("Arry shape", xyz_cl.atom.shape, disC.atom.shape)
    rng_m = get_angles(xyz_cl, disC, disC.atom, volC)

if __name__=="__main__":
    main()

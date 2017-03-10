import sys
import numpy as np
import itertools

from water_angles_util import translate_pbc, unit_vector, angle_between
from water_angles_util import plane_eq, calc_dip, project_plane
from xyzfile import XYZFile
from volfile import VolFile
from util    import d_pbc

WOXY = 1; WHYD = 2; GRAPHENE = 3

def cal_ang(w_coords, d_to_wall, rng):
    '''Given a list of coords (dimensions: [nat = 3, nwat, ndim = 3]), 
       move the coords into the nearest image of water of interest,
       calculate the dipole moments and angles''' 
    t1, t2, c1, c2, ph, r, wd = [], [], [], [], [], [], []
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
        wd.append(np.repeat(d_to_wall[i][1], ot_wr.shape[1]))

    t1 = list(itertools.chain(*t1));t2=list(itertools.chain(*t2));
    c1 = list(itertools.chain(*c1));c2=list(itertools.chain(*c2));
    ph = list(itertools.chain(*ph));r =list(itertools.chain(*r));
    wd = list(itertools.chain(*wd)); 
    return t1, t2, c1, c2, ph, r, wd

def get_angles(xyz, disC, volC):
    '''Method to get various angles between two waters'''
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

    # for x dir, use bins 0 -> x_quarter. ATs will be max L/2 from wall
    x_hlf = volC.get_x_max()/2.0 # Half of box divides 2 walls
    x_binz = np.arange(0.5, x_hlf, x_hlf/float(bnz)); 

    for i in range(1,len(xyz.atom)): # for each time snapshot, except first
        th1, th2, ch1, ch2, phi, rr, wd = [],[],[],[],[],[],[]
        rng = np.array([volC.get_x_rng_i(i), volC.get_y_rng_i(i),
                        volC.get_z_rng_i(i)]) # pbc range

        in_wat[0] = xyz.atom[i,oi,:]; in_wat[1] = xyz.atom[i,hs_i[0],:]
        in_wat[2] = xyz.atom[i,hs_i[1],:]
        t1, t2, c1, c2, ph, r, w = cal_ang(in_wat, disC.atom[i,oi,:], rng)
        th1+=t1; th2+=t2; ch1+=c1; ch2+=c2; phi+=ph; rr+=r; wd+=w

        ou_wat[0] = xyz.atom[i,oou,:]; ou_wat[1] = xyz.atom[i,hs_o[0],:]
        ou_wat[2] = xyz.atom[i,hs_o[1],:]
        t1, t2, c1, c2, ph, r, w = cal_ang(ou_wat, disC.atom[i,oou,:], rng)
        th1+=t1; th2+=t2; ch1+=c1; ch2+=c2; phi+=ph; rr+=r; wd+=w

        t1s += [th1];t2s += [th2];c1s += [ch1];c2s += [ch2];phs += [phi];
        rs += [rr];ws += [wd]
    return list([t1s, t2s, c1s, c2s, phs, rs, ws]) 

def print_angles(angls, fname):
    '''Print file of data in csv-like format, angls is a list of values:
       angls = [ [thet1],[thet2], [chi1], [chi2], [phi], [rs] ]'''
    f = open(fname, 'w'); 
    nsnap, stn = len(angls[0]), ''
    vals = ['the1_', 'the2_','chi1_', 'chi2_','phi_', 'dis_', 'walld_']
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
    disC = XYZFile("run"+nm+".dist", VolFile("")) 
    xyz_cl = XYZFile(xyzname, volC)

    angs = get_angles(xyz_cl, disC, volC)
    print_angles(angs, "run"+nm+"_angles.csv")

if __name__=="__main__":
    main()

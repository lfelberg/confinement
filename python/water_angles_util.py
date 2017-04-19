import numpy as np
import itertools
import operator

from util import d_pbc, translate_pbc

QOXY = -1.0484
QHYD =  0.5242

def find_closest(frm, to, n_close = 1):
    '''Given a vector of points dim frm = [nvals, nsamp, ndim = 3],
       and to = [nsamp, ndim], return a vector of closest = [size = nsamp]
       where each position has the index (0->nvals) of point closest '''
    closest = np.zeros((1, frm.shape[1])); nsamp = frm.shape[0]
   #to_ar = np.repeat(to[np.newaxis,:,:], nsamp, axis = 0) 
    to_ar = np.repeat(to[np.newaxis,:], nsamp, axis = 0) 
    d2 = np.sum((frm - to_ar)*(frm - to_ar), axis = -1)
    if n_close == 1: return np.argmin(d2, axis =  0)
    else: 
        d2 = list(d2); d2_l = list(range(len(d2)))
        d2 = list(zip(d2,d2_l)); d2.sort(key=operator.itemgetter(0))
        clos, dd = [], []
        for x in range(n_close): 
            clos.append(d2[x][1]); dd.append(np.sqrt(d2[x][0]))
        return clos, dd

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.sqrt(np.sum(vector*vector,axis=-1)[:,np.newaxis])

def angle_between(v1, v2, v1nrm = False, v2nrm = False):
    """ Returns the angle in radians between vectors 'v1' and 'v2'
        v1 and v2 have the dimensions = [dim, nwat]::
        >>> angle_between(np.array((1,0,0))[np.newaxis,:],np.array((0,1,0))[np.newaxis,:])
        array([ 1.57079633])
        >>> angle_between(np.array((1, 0, 0))[np.newaxis,:],np.array((1, 0, 0))[np.newaxis,:])
        array([ 0.])
        >>> angle_between(np.array((1, 0, 0))[np.newaxis,:],np.array((-1, 0, 0))[np.newaxis,:])
        array([ 3.14159265])
    """
    if v1nrm == True: v1_u = v1   #Normalize if needed
    else:             v1_u = unit_vector(v1)
    if v2nrm == True: v2_u = v2
    else:             v2_u = unit_vector(v2)

    return np.arccos(np.clip(np.sum(v1_u*v2_u, axis = -1), -1.0, 1.0))

def plane_eq(c1, c2, c3):
    '''Given three points, calculate a unit normal to the plane made by those 
       points: c1 = c2 = c3 = [dim, nwat] = cp, v1, v2, return
       Math: norm = (c2-c1) x (c3-c1) / || (c2-c1) x (c3-c1) ||
             where x=cross prod(a,b) = [a[1]b[2]-a[2]b[1], a[2]b[0]-a[0]b[2],
                                        a[0]b[1]-a[1]b[0]] ::
       >>> plane_eq(np.array((1,2,3))[:,np.newaxis],np.array((4,6,9))[:,np.newaxis],np.array((12,11,9))[:,np.newaxis])
       array([[-0.50760041],
              [ 0.81216065],
              [-0.28764023]])
       >>> plane_eq(np.array((0,0,0))[:,np.newaxis],np.array((1,0,0))[:,np.newaxis],np.array((0,11,0))[:,np.newaxis])
       array([[ 0.],
              [ 0.],
              [ 1.]])
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
       >>> project_plane(np.array((10.,11.,23.))[np.newaxis,:],np.array((0,1,0))[np.newaxis,:])
       array([[ 10.,   0.,  23.]])
       >>> project_plane(np.array((10.,11.,23.))[np.newaxis,:],np.array((1,1,1))[np.newaxis,:])
       array([[-4.66666667, -3.66666667,  8.33333333]])
    '''
    return (vec-(np.sum(vec*norm,axis=1)/
                 np.sum(norm*norm,axis=1))[:,np.newaxis]*norm)

def cal_ang(w_coords, rng, d_to_wall = [], dim = 3):
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
        cur_ar  = np.repeat(curr[np.newaxis,0,:], others.shape[1], axis = 0) 
        curr = curr[:,np.newaxis,:]
        cur_all = np.repeat(curr, others.shape[1], axis = 1) 
        ot_wr[0] = translate_pbc(cur_ar, others[0], rng) # trans other ox
        ot_wr[1] = translate_pbc(ot_wr[0], others[1], rng) # trans hydrogs
        ot_wr[2] = translate_pbc(ot_wr[0], others[2], rng)

        inter_mol_ax = ot_wr[0] - curr[0]
        mu_cur = calc_dip(curr); mu_oth = calc_dip(ot_wr)
        
        the_1 = angle_between(mu_cur, inter_mol_ax)
        the_2 = angle_between(mu_oth, inter_mol_ax)

        # calculate plane of intermol axis + dipole vector
        inter_o1_nrm = plane_eq(ot_wr[0].T, curr[0].T,(mu_cur+curr[0]).T).T
        inter_o2_nrm = plane_eq(curr[0].T,ot_wr[0].T,(mu_oth+ot_wr[0]).T).T

        # calc the vector between sorted hydrogens
        clos_h1_cur = find_closest(curr[1:], ot_wr[0])
        clos_h1_oth = find_closest(ot_wr[1:], curr[0])
        clos_h1_cur = np.repeat(clos_h1_cur[np.newaxis,:], 3, axis=0)
        clos_h1_oth = np.repeat(clos_h1_oth[np.newaxis,:], 3, axis=0)
        x,y = np.meshgrid(np.arange(cur_all.shape[1]),
                          np.arange(cur_all.shape[2]))
        hyd_cur=(cur_all[2-clos_h1_cur,x,y]-cur_all[clos_h1_cur+1,x,y]).T
        hyd_oth=(ot_wr[2-clos_h1_oth,x,y] - ot_wr[clos_h1_oth+1,x,y]).T

        chi1 = angle_between(hyd_cur, inter_o1_nrm, False, True)
        chi2 = angle_between(hyd_oth, inter_o2_nrm, False, True)
        
        mu_cur_proj = project_plane(mu_cur, inter_mol_ax)
        mu_oth_proj = project_plane(mu_oth, inter_mol_ax)
        phi = angle_between(mu_cur_proj, mu_oth_proj)
        
        if dim == 3: dists=d_pbc(curr[0],ot_wr[0],rng) # cal O-O distance
        else:        dists=d_pbc(curr[0,:,1:],ot_wr[0,:,1:],rng[1:],[1.,1.])

        t1.append(the_1); t2.append(the_2); c1.append(chi1); c2.append(chi2); 
        ph.append(phi); r.append(dists); 
        if d_to_wall!=[]: wd.append(np.repeat(d_to_wall[i][1],ot_wr.shape[1]))

    t1 = list(itertools.chain(*t1));t2=list(itertools.chain(*t2));
    c1 = list(itertools.chain(*c1));c2=list(itertools.chain(*c2));
    ph = list(itertools.chain(*ph));r =list(itertools.chain(*r));
    wd = list(itertools.chain(*wd));
    return t1, t2, c1, c2, ph, r, wd



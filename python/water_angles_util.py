import numpy as np

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

    return np.arccos(np.clip(np.sum(v1_u*v2_u, axis = 1), -1.0, 1.0))

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


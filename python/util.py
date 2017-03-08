import numpy as np


def d_pbc(c1, c2, rng, pbcs = [1.,1.,1.]):
    '''Compute distance between points with PBCS
       c1 and c2 are an array of coordinates, rng is PBC range in x,y,z'''
    pbcs = np.array(pbcs) # use this array to remove pbcs in any dimension
    boxl = np.round((c1-c2)/rng) # find what is rounded distance
    d = c1-c2 - boxl*rng*pbcs # add that rounded distance to actual dist
    ds = (d * d).sum(axis=1)
    return np.sqrt(ds)


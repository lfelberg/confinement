import numpy as np

def d_pbc(c1, c2, rng, pbcs = [1.,1.,1.]):
    '''Compute distance between points with PBCS
       c1 and c2 are an array of coordinates, rng is PBC range in x,y,z'''
    pbcs = np.array(pbcs) # use this array to remove pbcs in any dimension
    boxl = np.round((c1-c2)/rng) # find what is rounded distance
    d = c1-c2 - boxl*rng*pbcs # add that rounded distance to actual dist
    ds = (d * d).sum(axis=1)
    return np.sqrt(ds)

def d3(c1, c2):
    '''Compute distance between points c1 and c2 are an array of coordinates'''
    d = c1-c2
    return d * d

def translate_pbc(c1, c2, rng):
    '''Translate c2 to distance sub rng/2 wrt c1
       c1 and c2 are arrays of coordinates, rng is PBC range in x,y,z'''
    boxl = np.round((c1-c2)/rng) # find what is rounded distance
    return c2 + boxl*rng

def translate_1st_im(c, rng):
    '''Translate coords to first periodic image with box rng: 0 -> edge'''
    return c - np.floor(c/rng)*rng

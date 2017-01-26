import sys
import re
import numpy as np
from numpy import linalg as LA

def main():
    filename=sys.argv[1]
    sigma=float(sys.argv[2])
    ats, typs = atomfinder(filename)
    grids = gridfind(sigma)
    voidfind(ats, typs, grids, sigma)
        
def atomfinder(filename):
    ''' Method to create array of atom coordinates in xyz file'''
    f=open(filename, "r")
    f.readline()
    f.readline()
    print filename
    atoms, atom, typ, types = [], [], [], []
    
    '''First, we read in the coordinates from the trajectory file and 
       save them into atom=[[coords],[coords]...,[coords]]'''
    for line in f:
        #Gets past # and atoms lines
        if line[0]=="A":
            atom+=[atoms]; typ += [types]
            atoms, types = [], []
        # Only grab graphene carbons and water oxygens
        elif len(line.split()) > 1 and line.split()[0] != '4' and \
             line.split()[0] != '5' and line.split()[0] != '2':
            #coords is a list of the coordinates in string formats        
            coords=line.split()
            types.append(int(coords[0]))
            atoms.append([float(coords[1]),float(coords[2]),float(coords[3])])
    #For last timestep, since it doesn't have a line about Atoms after it.
    atom=atom+[atoms]; typ = typ+[types]
    f.close()
    return atom, typ

def gridfind(sigma):
    '''This function returns a list of lists of coordinates. The 
       coordinates are of the gridpoints at distance=0.5*sigma'''
    # Divide up the entire volume into 0.5x0.5x0.5 
    x=np.arange(-6,18,sigma)
    y=np.arange(0, 46.86,sigma)
    z=np.arange(0,49.19,sigma)
    xx,yy,zz = np.meshgrid(x,y,z)
    xx = xx.flatten(); yy = yy.flatten(); zz = zz.flatten()
    grid = np.zeros((xx.shape[0], 3))
    grid[:,0] = xx; grid[:,1] = yy; grid[:,2] = zz;

    return grid

def voidfind(atoms, types, gridpoints, sigma):
    '''Given a grid and a list of atoms from a trajectory, identify the 
       voids'''
    counter=0
    grid_pres = np.zeros((len(atoms),gridpoints.shape[0]))
    print("Gdpts shape {0}, and {1}".format(gridpoints.shape, gridpoints[0]))
    for item in atoms: #each item is one different TIMESTEP
        print("Timestep {0} in voidfind".format(counter))
        counter+=1; coordCt = 0
        for coords in item:
            rad = 3.35 if types[counter-1][coordCt] == 1 else 3.35
            coord = np.array(coords)
            # making an array of atom point gridshape times to find dist to
            # each grid point
            atomGd = np.repeat(coord[np.newaxis,:],gridpoints.shape[0],axis=0)
      
            diff = atomGd - gridpoints
            dists = LA.norm(diff, axis = 1)
            inds = dists < rad
            grid_pres[counter-1][inds] += 1
            coordCt += 1

    voidpoints = []
    for ct in range(len(atoms)):
        vInds = grid_pres[ct] == 0
        voidpoints.append(gridpoints[vInds])

    printVoids(voidpoints, "test.xyz")

def printVoids(voidpoints, fname):
    '''Method to print the locations of void points in xyz file'''
    maxVd = 0
    f = open(fname, 'w')
    for ct in range(len(voidpoints)): # Finding timestep w most voids
        maxVd = len(voidpoints[ct]) if len(voidpoints[ct]) > maxVd else maxVd
    for step in range(len(voidpoints)):
        f.write("{0}\n".format(maxVd))
        f.write("Atoms\n")
        for void in range(len(voidpoints[step])):
            f.write("9 {0} {1} {2}\n".format(
                     *tuple(voidpoints[step][void])))
        if len(voidpoints[step]) < maxVd:
            for vd in range(len(voidpoints[step]), maxVd):
                f.write("23 0.0 0.0 0.0\n")

    f.close()

if __name__=="__main__":
    main()


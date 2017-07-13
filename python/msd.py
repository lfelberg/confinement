import sys
import re
import math
import numpy as np

from xyzfile import XYZFile
from volfile import VolFile
from util    import d3, translate_pbc

OXY = 1

def trans_coords(coords, rng):
    ''' Method to translate coordinates so that each step is within 1 box
        of the previous. This is to account for the total MSD. '''
    cnew = np.zeros(coords.shape); cnew[0] = coords[0]
    st = "{0}\nAtoms".format(coords.shape[1]) # for printing XYZ file

    for i in range(len(coords)-1):
        cnew[i+1] = translate_pbc(cnew[i], coords[i+1], rng[i])
        cnew[i] = cnew[i] - np.mean(cnew[i], axis=0)
       #print(st)
       #for j in range(len(coords[0])): print("{0} {1} {2} {3}".format(8,*cnew[i,j,:]))
    cnew[-1] = cnew[-1] - np.mean(cnew[-1], axis=0)
    return cnew

def get_msd(xyz, volC):
    '''Method to get the mean square displacement of water oxygens and solutes'''
    nsnaps = len(xyz.atom)-1 # number of snaps in XYZ, remove 0th frame
    # find water oxygen atoms
    nOX = xyz.get_ct_i(OXY); nmol = xyz.nsol + nOX
    iOX,_ = xyz.get_inner_wat(); oOX,_ = xyz.get_outer_wat()
    if xyz.nsol > 0:
        iSL = xyz.get_inner_sol(); oSL = xyz.get_outer_sol()
    oicd = trans_coords(xyz.atom[1:,iOX],volC.get_rng()[1:]) # don't use 1st snap
    oocd = trans_coords(xyz.atom[1:,oOX],volC.get_rng()[1:]) # don't use 1st snap
    sicd = trans_coords(xyz.atom[1:,iSL],volC.get_rng()[1:]) # don't use 1st snap
    socd = trans_coords(xyz.atom[1:,oSL],volC.get_rng()[1:]) # don't use 1st snap

    msd = np.zeros((nsnaps,nmol,3)) #save msd for each oxy and solute molecule
    ms_mean = np.arange(nsnaps,0,-1.)[:, np.newaxis, np.newaxis]
    for i in range(nsnaps-1):
        ms = d3(oicd[np.newaxis,i], oicd[i+1:])
        msd[1:nsnaps-i,:sum(iOX.astype(int))] += ms
        if sum(oOX.astype(int)) > 0: # will be zero for bulk water
            ms = d3(oocd[np.newaxis,i], oocd[i+1:])
            msd[1:nsnaps-i,sum(iOX.astype(int)):nOX] += ms
        if xyz.nsol > 0: # If there are solvent molecules
           ms = d3(sicd[np.newaxis,i], sicd[i+1:])
           msd[1:nsnaps-i,nOX:nOX+(xyz.nsol/2)] += ms
           ms = d3(socd[np.newaxis,i], socd[i+1:])
           msd[1:nsnaps-i,nOX+(xyz.nsol/2):] += ms
    return msd/ms_mean, nOX

def print_msd(msd, timstep, fname, nOX):
    '''Print the mean square displacement for for y-z averaged across oxys and for each solute'''
    f = open(fname, 'w'); nsol = msd.shape[1]-nOX
    f.write("ASTP(PS),MSDYZ,MSDYZER")
    for i in range(nsol): f.write(",MSDYZ_sol{0}".format(i))
    f.write("\n")

    tstep = np.arange(0,timstep*msd.shape[0], timstep); d = np.zeros(2)
    for i in range(msd.shape[0]):
        dat = msd[i]; st = ""
        d[0]  = np.mean(np.sum(dat[:nOX,1:], axis = 1))
        d[1] = np.std(np.sum(dat[:nOX,1:], axis = 1))
        for j in range(len(d)): st += "{0:.5f},".format(d[j])
        for s in range(nsol): st += "{0:.5f},".format(np.sum(dat[s,1:]))
        f.write("{0},{1}\n".format(tstep[i]*2./1000., st[:-1]))
    f.close()

def main():
    xyzname=sys.argv[1]; sep=sys.argv[2]; ln=sys.argv[3]; itr=sys.argv[4]
    sol_nm = sys.argv[5]; sol_ct = sys.argv[6]
    nm = str(sep)+"_"+str(ln)+"_"+itr; volC = VolFile("run"+nm+".vol") 
    xyz_cl = XYZFile(xyzname, volC)
    xyz_cl.sol_ty = sol_nm
    xyz_cl.nsol = int(sol_ct)

    msd, nw = get_msd(xyz_cl, volC)
    print_msd(msd, xyz_cl.time[3]-xyz_cl.time[2], 'msd_'+nm+".csv", nw)

if __name__=="__main__":
    main()

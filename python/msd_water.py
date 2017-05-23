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
    cnew = np.zeros((coords.shape[0],1,3)); #cnew = np.zeros(coords.shape); tcnew[0] = coords[0]
    cnew[0] = np.mean(cnew[0], axis=0)# for calc msd of com
   #st = "{0}\nAtom".format(len(coords[0]))
   #for i in range(len(coords)-1):
    for i in range(1,len(coords)):
       #cnew[i+1] = translate_pbc(cnew[i], coords[i+1], rng[i])
        cnew[i] = np.mean(translate_pbc(cnew[i-1], coords[i], rng[i-1]), axis=0) # for calc msd of com
       #cnew[i] = cnew[i] - np.mean(cnew[i], axis=0)
       #print(st)
       #for j in range(len(coords[0])): print("{0} {1} {2} {3}".format(8,*coords[i,j,:]))
   #cnew[-1] = cnew[-1] - np.mean(cnew[-1], axis=0)
   #print(st)
   #for j in range(len(coords[0])): print("{0} {1} {2} {3}".format(8,*coords[-1,j,:]))
    return cnew

def get_msd(xyz, volC):
    '''Method to get the mean square displacement of water oxygens'''
    # find water oxygen atoms
    tOX = xyz.get_type_i(OXY)
    nsnaps = len(xyz.atom)-1 # number of snaps in XYZ, remove 0th frame
    nOX = np.sum(tOX.astype(int))
   #msd = np.zeros((nsnaps,nOX,3)) # saving msd for each oxygen
    msd = np.zeros((nsnaps,1,3)) # msd of com
    ms_mean = np.arange(nsnaps,0,-1.)[:, np.newaxis, np.newaxis]

    ocrd = trans_coords(xyz.atom[1:,tOX],volC.get_rng()[1:]) # don't use 1st snap
    for i in range(nsnaps-1):
        ms = d3(ocrd[np.newaxis,i], ocrd[i+1:])
        print(ocrd[np.newaxis,i].shape, ocrd[i+1:].shape, ms.shape)
        msd[1:nsnaps-i] += ms
    return msd/ms_mean

def print_msd(msd, timstep, fname):
    '''Print the mean square displacement for 3D, for y-z and for x'''
    f = open(fname, 'w'); 
    f.write("ASTP(PS),MSDX,MSDXER,MSDY,MSDYER,MSDZ,MSDZER,MSDYZ,MSDYZER,MSD3D,MSD3DER\n")
    tstep = np.arange(0,timstep*msd.shape[0], timstep)
    d = np.zeros(10)
    for i in range(msd.shape[0]):
        dat = msd[i]; st = ""
        d[0] = np.mean(dat[:,0]); d[1] = np.std(dat[:,0])
        d[2] = np.mean(dat[:,1]); d[3] = np.std(dat[:,1])
        d[4] = np.mean(dat[:,2]); d[5] = np.std(dat[:,2])
        d[6]  = np.mean(np.sum(dat[:,1:], axis = 1))
        d[7] = np.std(np.sum(dat[:,1:], axis = 1))
        d[8] = np.mean(np.sum(dat, axis = 1))
        d[9] = np.std(np.sum(dat, axis = 1))
        for j in range(len(d)): st += "{0:.5f},".format(d[j])
        f.write("{0},{1}\n".format(tstep[i]*2./1000., st[:-1]))
    f.close()

def main():
    ''' Given a list of pairs for g(r) calcs, do this as a function as 
        distance from wall, but only 2D and only for your wall side '''
    xyzname=sys.argv[1]; sep=sys.argv[2]; ln=sys.argv[3]; itr=sys.argv[4]

    nm = str(sep)+"_"+str(ln)+"_"+str(itr)
    volC = VolFile("run"+nm+".vol") 
    xyz_cl = XYZFile(xyzname, volC)

    msd = get_msd(xyz_cl, volC)
    print_msd(msd, xyz_cl.time[3]-xyz_cl.time[2], 'msd_'+nm+".csv")

if __name__=="__main__":
    main()

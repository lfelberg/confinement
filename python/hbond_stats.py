import sys
import re
import numpy as np
from xyzfile import XYZFile
from volfile import VolFile
from csvfile import CSVFile

GRAPHENE = 3
WATER_OXY = 2

''' 
H_HD_HA,frame_1,frame_2,frame_3,frame_4,frame_5,,Average
4817_4816_2584,1,1,1,1,1,,1.0
1938_1936_4123,1,1,1,1,1,,1.0

sorted: 'Average', 'H', HD, HA, 'frame_1', 'frame_2', 'frame_3', 'frame_4', 'frame_5'
so the frames start at 4
'''

def get_hbond_v_dist(xyzC, sep, csvC):
    '''Using output of hbonanza to tabulate hydrogen bonding in graphene'''
    xr, bns, fr_s = sep/2., 100, 4 # where the frame starts

    # histogram of number of hbonds per dist
    dens = np.zeros((len(xyzC.atom), bns))   
    print(xyzC.atom.shape)
    n_wat = sum(xyzC.types == WATER_OXY)
    wats = xyzC.types == WATER_OXY
    for ti in range(len(xyzC.atom)):
        h_ct = np.zeros(len(xyzC.atom[ti]))
        for at in range(len(xyzC.atom[ti])):
            if xyzC.types[at] == WATER_OXY:
                fr = np.zeros((2))
                hd = csvC.dat[1].astype(int) == at
                ha = csvC.dat[2].astype(int) == at
                fr[0] = sum(csvC.dat[ti+fr_s, hd])
                fr[1] = sum(csvC.dat[ti+fr_s, ha])
                h_ct[at] = sum(fr)
        hw, benz = np.histogram(xyzC.atom[ti,:,-1],bins=bns,range=(0.0,xr), weights=h_ct)
        hs, _ = np.histogram(xyzC.atom[ti,:,-1],bins=bns,range=(0.0,xr), weights=wats.astype(float))
        
        hw[hw == 0.0] = 1.0
        print(hw, hs)
        dens[ti] = hs.astype(float)/hw.astype(float)
    dens_mn = np.mean(dens, axis = 0)

    f = open(xyzC.xyzfname[:-3]+"hbond_dat", 'w')
    f.write("Bin,hbond\n")                                 
    for i in range(bns):                                                           
        f.write("{0:>5.3f},{1:.4f}\n".format(              
                 benz[i], dens_mn[i]))                                   
    f.close() 


def main():
    ''' usage: exe *.xyz sep iter '''
    xyzname = sys.argv[1]; sep = sys.argv[2]; itr = sys.argv[3]
    volC = VolFile("run"+str(sep)+"_"+str(itr)+".vol")
    xyzC = XYZFile("run"+str(sep)+"_"+str(itr)+".dist", VolFile(''))
    print("Arry shape", xyzC.atom.shape)
    csvC = CSVFile(xyzname[:-3]+'frame_by_frame_hbonds.csv')
    
    print("This is csvC.dat shape {0}".format(csvC.dat.shape))
    get_hbond_v_dist(xyzC, float(sep), csvC)

if __name__=="__main__":
    main()

import sys
import numpy as np

from volfile import VolFile
from xyzfile import XYZFile

typ_dic = { 
            1: ['O WAT', -1.0484,], 
            2: ['H WAT',  0.5242,], 
            3: ['C GRA',  0.0,], 
            4: ['C BEN', -0.115,], 
            5: ['H BEN',  0.115,],
          }

def write_to_pdb(fname, xyz):
    '''Given coordinates, PBC box dimensions, and type list,
       write out a pdb trajectory file'''
    f = open(fname, 'w')
    for tim in range(len(xyz.atom)):
        res_ct = 0
        for i in range(len(xyz.types)):
            if xyz.types[i] == 1: res_ct += 1
            elif xyz.types[i] == 3 and xyz.types[i-1] == 2: res_ct += 1
            elif xyz.types[i] == 4 and xyz.types[i-1] != xyz.types[i]: res_ct += 1

            f.write("ATOM{0:7d}    {1:s} {2:5d}    {3:8.3f}{4:8.3f}{5:8.3f} {ll:5.2f}  0.00\n".format(
                     i, typ_dic[xyz.types[i]][0], res_ct,
                     *xyz.atom[tim, i, :], 
                     ll=typ_dic[xyz.types[i]][1]))
        f.write("END\n")
    f.close()

def main():
    xyzname=sys.argv[1]; sep=sys.argv[2]; ln=sys.argv[3]; itr=sys.argv[4]
    nm = str(sep)+"_"+str(ln)+"_"+str(itr)
    volC = VolFile("run"+nm+".vol")
    xyzC = XYZFile(xyzname, volC)

    write_to_pdb(xyzname[:-3]+"pdb", xyzC)

if __name__=="__main__":
    main()

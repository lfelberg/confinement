import sys
import re
import numpy as np

from volfile import VolFile
from xyzfile import XYZFile
WOXY = 1; WHYD = 2; GRAPHENE = 3

type_dic = {
             1 : ["WATER", "OW"],
             2 : ["WATER", "HW"],
             3 : ["GRAPH", "C"],
           }

class ZipFile:
    '''A class for combining velocity and xyz files'''

    def __init__(self, fname, volC):
         self.fpref = fname
         self.zip_up(volC)

    def zip_up(self, volC):
        '''A method to take velocity file, xyzfile and volClass and make gro file'''
        self.xyzC = XYZFile(self.fpref+".xyz", volC)
        self.velC = XYZFile(self.fpref+".velxyz", volC)
        self.print_coords(volC)

    def print_coords(self, volC):
        '''If given coordinates, write them to file'''
        f, nat = open("traj_"+volC.volfname[3:-4]+".gro", 'a'), 0

        for ti in range(len(self.xyzC.atom)):
            if nat == 0: nat = len(self.xyzC.atom[ti])
            if nat!=len(self.xyzC.atom[ti]) : print("Different number of atoms")

            resno = 1
            f.write("confined syst, t= {0}\n{1}\n".format(volC.time[ti,1],nat))
            ## "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f"
            for i in range(nat):
                if self.xyzC.types[i]==WOXY or self.xyzC.types[i]==GRAPHENE: resno += 1
                strn = "{0:5}{1:>5}{2:5}{3:5}".format(resno,
                        type_dic[self.xyzC.types[i]][0],
                        type_dic[self.xyzC.types[i]][1], i)

                strn += "{0:8.3f}{1:8.3f}{2:8.3f}".format(*(self.xyzC.atom[ti][i]))
                strn += "{0:8.4f}{1:8.4f}{2:8.4f}".format(*(self.velC.atom[ti][i]))
                strn+="\n"
                f.write(strn)

            f.write("{0:.5f} {1:.5f} {2:.5f}\n".format(*(volC.get_rng_i(ti))))

        f.close()

def main():
    filename=sys.argv[1]; volf=sys.argv[2]
    vC = VolFile(volf)
    zf = ZipFile(filename, vC)

if __name__=="__main__":
    main()


import sys
import re
import numpy as np

'''0. Step 1.Time 2.Temp 3.Press 4.PotEng 5.KinEng 6.TotEng 7.E_vdwl 
   8.E_coul 9.E_pair 10.E_bond 11.E_angle 12.E_dihed 13.Volume 14.Density 
   15.Lx 16.Ly 17.Lz 18.Xlo 19.Xhi 20.Ylo 21.Yhi 22.Zlo 23.Zhi '''


def get_box_dim(filename):
    '''Method to get the dimensions of box for xyz file
       from lammps .out file'''
    f=open(filename, "r")
    time, dims = [], []
    for line in f:
       if (len(line) > 200) and bool(re.search(r'\d', line)) == True:
           tmp = line.split()
           time.append([int(tmp[0]), int(tmp[1])])
           dim = [float(tmp[x]) for x in range(18, 24)]
           dims.append(dim)

    f.close()
    return time, dims


def print_box_dim(time, dims, vol_out):
    '''For all times in run file, print x, y, z min and max'''
    f = open(vol_out, "w")
    for i in range(len(time)):
        f.write("{0} {1} {2} {3} {4} {5} {6} {7}\n".format(
                *tuple(time[i] + dims[i])))
    f.close()

def main():
    filename=sys.argv[1]
    time, dims = get_box_dim(filename)
    print_box_dim(time, dims, filename[:-3]+"vol")


if __name__=="__main__":
    main()

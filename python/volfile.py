import sys
import re
import numpy as np

'''0. Step 1.Time 2.Temp 3.Press 4.PotEng 5.KinEng 6.TotEng 7.E_vdwl 
   8.E_coul 9.E_pair 10.E_bond 11.E_angle 12.E_dihed 13.Volume 14.Density 
   15.Lx 16.Ly 17.Lz 18.Xlo 19.Xhi 20.Ylo 21.Yhi 22.Zlo 23.Zhi '''

class VolFile:
    '''A class for volume files'''
    volfname = ''
        
    def __init__(self, fname):
         self.volfname = fname

         if '.vol' in fname:   self.get_dim_from_vol(fname)
         elif '.out' in fname: self.get_dim_from_out(fname)
         else:                 print("Error, file ext not recognized!")


    def get_dim_from_vol(self, fname):
        '''Method to get dimensions from volume file'''
        f = open(fname, "r")                                                         
        time, dims = [], []                                                            
        for line in f:                                                                 
            tmp = line.split()                                                         
            time.append([int(tmp[0]), int(tmp[1])])                                    
            dims.append([float(tmp[2]),float(tmp[3]),float(tmp[4]),                    
                         float(tmp[5]),float(tmp[6]),float(tmp[7])])                   
        f.close()                                                                      
        self.time = np.array(time)
        self.dims = np.array(dims)

    def get_dim_from_out(self, filename):
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

        self.time = np.array(time)
        self.dims = np.array(dims)

    def print_box_dim(self, vol_out):
        '''For all times in run file, print x, y, z min and max'''
        f = open(vol_out, "w")
        for i in range(len(self.time)):
            f.write("{0} {1} {2} {3} {4} {5} {6} {7}\n".format(
                    *(self.time[i].tolist()+ self.dims[i].tolist())))
        f.close()

    def get_x_rng(self):
       '''Return range of x for first snap'''
       return (self.dims[0][1] - self.dims[0][0])


def main():
    filename=sys.argv[1]
    vC = VolFile(filename)
    vC.print_box_dim(filename[:-3]+"vol")

if __name__=="__main__":
    main()


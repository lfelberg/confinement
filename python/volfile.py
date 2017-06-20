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
         else:                 self.volfname = ""


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
              #if ((float(tmp[0])<5000000) and (float(tmp[0])>=4000000) and (float(tmp[0])%500==0)) or \
               if ((float(tmp[0])>=2500000) and (float(tmp[0])%500==0)) or \
                    int(tmp[0])==0:
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

    def print_quarters(self, volf_pref):
        '''For times in run file print out quarters'''
        print(len(self.time))
        for q in range(4):
            f = open(volf_pref+"_q"+str(q+1)+".vol", "w")
            # always include the first snapshot
            f.write("{0} {1} {2} {3} {4} {5} {6} {7}\n".format(
                    *(self.time[0].tolist()+ self.dims[0].tolist())))
            start = q*2500+1; stop = (q+1)*2500+1
            for i in range(start,stop):
                f.write("{0} {1} {2} {3} {4} {5} {6} {7}\n".format(
                        *(self.time[i].tolist()+ self.dims[i].tolist())))
            f.close()


    def get_x_max(self):
        '''Return range of x for first snap'''
        xdif = self.dims[:,1] - self.dims[:,0]
        return max(xdif)

    def get_x_rng(self):
        '''Return range of x for first snap'''
        return (self.dims[0][1] - self.dims[0][0])

    def get_rng(self):
        '''Return range of all snaps'''
        return np.array([self.dims[:,1]-self.dims[:,0],
                         self.dims[:,3]-self.dims[:,2],
                         self.dims[:,5]-self.dims[:,4]]).T

    def get_x_rng_i(self, i):
        '''Return range of x for ith snap'''
        return (self.dims[i][1] - self.dims[i][0])

    def get_y_rng_i(self, i):
        '''Return range of y for ith snap'''
        return (self.dims[i][3] - self.dims[i][2])

    def get_z_rng_i(self, i):
        '''Return range of z for ith snap'''
        return (self.dims[i][5] - self.dims[i][4])

    def get_rng_i(self, i):
        '''Return array of range for xy and z for step i'''
        return np.array([self.get_x_rng_i(i), self.get_y_rng_i(i),
                         self.get_z_rng_i(i)])

    def get_vol_i(self, i):
        '''Return volume for ith snap'''
        return self.get_x_rng_i(i)*self.get_y_rng_i(i)*self.get_z_rng_i(i)

    def get_vol(self):
        '''Return average volume'''
        return self.get_x_len()*self.get_y_len()*self.get_z_len()

    def get_x_len(self):
        '''Return average length of x dim'''
        xdif = self.dims[:,1] - self.dims[:,0]
        return np.mean(xdif)

    def get_y_len(self):
        '''Return average length of y dim'''
        ydif = self.dims[:,3] - self.dims[:,2]
        return np.mean(ydif)

    def get_z_len(self):
        '''Return average length of y dim'''
        zdif = self.dims[:,5] - self.dims[:,4]
        return np.mean(zdif)


def main():
    filename=sys.argv[1]
    vC = VolFile(filename)
    vC.print_quarters("run"+filename[3:-6])
    vC.print_box_dim("run"+filename[3:-4]+".vol")

if __name__=="__main__":
    main()


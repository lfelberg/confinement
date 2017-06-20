import sys
import re
import numpy as np

'''0. Step 1.Time 2.Temp 3.Press 4.PotEng 5.KinEng 6.TotEng 7.E_vdwl 
   8.E_coul 9.E_pair 10.E_bond 11.E_angle 12.E_dihed 13.Volume 14.Density 
   15.Lx 16.Ly 17.Lz 18.Xlo 19.Xhi 20.Ylo 21.Yhi 22.Zlo 23.Zhi 24.Pxx 
   25.Pyy 26.Pzz 27.Pxy 28.Pxz 29.Pyz 30.c_2 31.c_4
'''

class Energy:
    '''A class for energy data'''
    enfname = ''
        
    def __init__(self, fname):
         self.enfname = fname
         if '.vol' in fname:    self.get_dim_from_en(fname)
         elif '.out' in fname: self.get_dim_from_out(fname)
         else:                 self.enfname = ""

    def get_dim_from_en(self, fname):
        '''Method to get energy from en file'''
        f = open(fname, "r"); time, dims = [], []                                                            
        for line in f:                                                                 
            tmp = line.split()                                                         
            if int(tmp[0]) == 0: continue
            time.append(int(tmp[0]))
           #dims.append([float(tmp[1]),float(tmp[2])])
            dims.append([float(tmp[3])-float(tmp[2]), float(tmp[5])-float(tmp[4]),
                         float(tmp[7])-float(tmp[6])])
        f.close()
        self.time = np.array(time); self.dims  = np.array(dims)

    def get_dim_from_out(self, filename):
        '''Method to get the energies of system
           from lammps .out file'''
        f=open(filename, "r"); time, dims = [], []
        for line in f:
           if (len(line) > 200) and bool(re.search(r'\d', line)) == True:
               tmp = line.split()
               if ((float(tmp[0])>=3500000) and (float(tmp[0])%1000==0)) :
                   time.append(int(tmp[0]))
                   dims.append([float(tmp[15]),float(tmp[16]),float(tmp[17])])
        f.close()

        self.time = np.array(time); self.dims = np.array(dims)

    def print_energies(self, en_out):
        '''For all times in run file, print c2 and c4 energies'''
        ANG_TO_METER = 1e-10
        A2_TO_M2 = ANG_TO_METER * ANG_TO_METER
        ym = np.mean(self.dims[:,1]);#print(np.shape(self.dims[:,1]), ym)
        zm = np.mean(self.dims[:,2]);#print(np.shape(self.dims[:,2]), zm)
        xx = self.dims[:,0]
        dx2 = (np.mean(np.power(xx,2)) - np.power(np.mean(xx),2))*A2_TO_M2
        xm = np.mean(xx)*ANG_TO_METER
        a2 = ym*zm*A2_TO_M2  # Area in m^2
        kT = 4.11e-21 # kT at 298 K in J
        print("{0:>4s},{1:.7f}".format(self.enfname.split("_")[0][3:],
                                     (a2*dx2*1e11)/(kT*xm)))
        f = open(en_out, "w"); f.write("Atim,lx\n")
        for i in range(len(self.time)):
            f.write("{0},{1}\n".format(self.time[i],*self.dims[i]))
        f.close()

def main():
    filename=sys.argv[1]; enC = Energy(filename)
    enC.print_energies(filename[:-3]+"lx")

if __name__=="__main__":
    main()


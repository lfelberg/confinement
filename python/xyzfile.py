import sys
import numpy as np

from volfile import VolFile

WOXY = 1; WHYD = 2; GRAPHENE = 3

class XYZFile:
    '''A class for xyz files'''
    xyzfname = ''

    def __init__(self, fname, VolFile = VolFile(""), xyz = [], ty = [],
                 nsol = 0, sol_ty = ""):
        self.xyzfname = fname
        self.nsol = nsol
        self.sol_ty = sol_ty
        if xyz == []:
            if VolFile.volfname == '': self.get_coords_types(fname,[],[])
            else:
                self.get_coords_types(fname,VolFile.time,VolFile.dims)
                self.half_x = VolFile.get_x_rng_i(0)/2.0
        else:
            self.print_coords(fname, xyz, ty)

    def get_coords_types(self, filename, times, dims):
        '''Method to open xyz file and save coords and types'''
        time, atoms, atom, types = [], [], [], []
        grab_snap, t_ct, t_grab = 0, 0, 0
        f=open(filename, "r")
        f.readline()
        line = f.readline()
        if times == [] or times[0][0] == int(line.split()[-1]):
            time.append(int(line.split()[-1]))
            t_ct += 1; grab_snap = 1

        self.edges = dims[:,1::2] - dims[:,0::2]
        '''First, we read in the coordinates from the trajectory file and
           save them into atom=[[coords],[coords]...,[coords]]'''
        for line in f:
            if times != [] and t_ct >= len(times): break

            if ("ITEM" not in line and len(line.split())==1 and t_grab==1 
                and times[t_ct][0]==int(line)):
                time.append(int(line)); t_grab = 0
                t_ct += 1; grab_snap = 1
            if "ITEM: TIMESTEP" in line: t_grab = 1

            #Gets past # and atoms lines
            if line[0]=="A" or "ITEM: ATOMS" in line: # for vel or xyz file
                if grab_snap == 1:
                    if len(atoms) > 0: atom+=[atoms]
                    atoms = []
                if "vel" not in self.xyzfname:
                    if (times==[] or times[t_ct][0] == int(line.split()[-1])):
                        time.append(int(line.split()[-1]))
                        t_ct += 1; grab_snap = 1
                    else: grab_snap = 0

            # Only grab snapshots that have volume data
            elif len(line.split()) > 3 and "ITEM" not in line and grab_snap == 1:
                #crds is a list of the coordinates in string formats
                crds=line.split()
                if t_ct == 1: types.append(int(crds[0]))
                if dims == []:
                    tmp = []
                    for i in range(1,len(crds)): tmp.append(float(crds[i]))
                    atoms.append(tmp)
                else: # Move coords so that box is from 0 -> xmin
                    atoms.append([float(crds[-3])-dims[t_ct][0],
                                  float(crds[-2])-dims[t_ct][2],
                                  float(crds[-1])-dims[t_ct][4]])
        #For last timestep, since it doesn't have a line about Atoms after it.
        if len(atoms) == len(atom[-1]): atom += [atoms]
        self.time  = np.array(time); self.atom = np.array(atom)
        self.types = np.array(types)
        print("xyz", filename, len(times), self.atom.shape, self.types.shape)
        f.close()

    def print_coords(self, fname, xyz = [], ty = []):
        '''If given coordinates, write them to file'''
        if xyz == []: xyz = self.atoms
        if ty == []:  ty  = self.types
        f, nat = open(fname, 'a'), 0

        for ti in range(len(xyz)):
            if nat == 0: nat = len(xyz[ti])
            if nat != len(xyz[ti]) : print("Different number of atoms")
            f.write("{0}\nAtoms\n".format(nat))
            for i in range(nat):
                strn = "{0}".format(ty[i])
                for j in range(len(xyz[ti][i])):
                    strn+=" {0}".format(xyz[ti][i][j])
                strn+="\n"
                f.write(strn)

        f.close()

    def get_nsnap(self):
        '''Get the number of timeshots'''
        return len(self.atom)

    def get_ct_i(self, i):
        '''Get count of type i'''
        return np.sum((self.types == i).astype(int))

    def get_type_i(self, i):
        '''Get indices of type i'''
        return self.types == i

    def get_graph_wall(self, wall_no=0):
        '''Get the indices of graphene wall 0 or 1'''
        grap = self.types == GRAPHENE
        if wall_no == 0:
            xl = self.atom[0,:,0] < self.half_x*0.95
            return np.all(np.array([grap, xl]), axis=0) # first wall
        else:
            xg = self.atom[0,:,0] > self.half_x*0.95
            return np.all(np.array([grap, xg]), axis=0) # second wall

    def get_wall_i_xv(self, i=0):
        '''Get the x location of graphene wall 0 or 1, from first snap'''
        wl_ids = self.get_graph_wall(i)
        return np.mean(self.atom[0,wl_ids,0])

    def get_inner_ats(self):
        '''Get indices of atoms whose x values are btw the 2 walls'''
        w0 = self.get_wall_i_xv(0)
        w1 = self.get_wall_i_xv(1)
        x_m0  = self.atom[0,:,0] > w0; x_m1 = self.atom[0,:,0] < w1

        return np.all(np.array([x_m0, x_m1]), axis=0)

    def get_outer_ats(self):
        '''Get indices of atoms whose x values are outside the 2 walls'''
        w0 = self.get_wall_i_xv(0)
        w1 = self.get_wall_i_xv(1)
        x_m0  = self.atom[0,:,0] > w0; x_m1 = self.atom[0,:,0] < w1

        wall = np.all(np.array([x_m0, x_m1]), axis=0)
        return wall == False

    def get_inner_wat(self):
        '''Get array of water atoms (ox+h) that are inside 2 walls'''
        if "_37_" in self.xyzfname: return self.get_type_i(WOXY),self.get_type_i(WHYD)

        inside = self.get_inner_ats()
        oxy = self.get_type_i(WOXY)
        hyd = self.get_type_i(WHYD)
        return np.all(np.array([oxy, inside]),axis=0), np.all(np.array([hyd, inside]),axis=0)

    def get_outer_wat(self):
        '''Get array of water atoms (ox+h) that are outside 2 walls'''
        if "_37_" in self.xyzfname: return [False]*len(self.types),[False]*len(self.types)
        out = self.get_outer_ats()
        oxy = self.get_type_i(WOXY)
        hyd = self.get_type_i(WHYD)
        return np.all(np.array([oxy, out]),axis=0), np.all(np.array([hyd, out]),axis=0)

    def get_inner_sol(self):
        '''Get array of solute atoms that are inside 2 walls'''
        if self.nsol == 0: return [False]*len(self.types)
        if "_37_" in self.xyzfname: return self.get_sol()

        inside = self.get_inner_ats()
        sl = self.get_sol()
        return np.all(np.array([sl, inside]),axis=0)

    def get_outer_sol(self):
        '''Get array of solute atoms that are outside 2 walls'''
        if self.nsol == 0: return [False]*len(self.types)
        if "_37_" in self.xyzfname: return [False]*len(self.types)

        outside = self.get_outer_ats()
        sl = self.get_sol()
        return np.all(np.array([sl, outside]),axis=0)

    def get_sol(self):
        '''Depending on the type of solute, for now either ion or benz,
           return either both types 4 and 5 (ion) or just 1 atom from each mol
           (benz)'''
        if self.sol_ty == "ion":
           ion1 = self.get_type_i(4)
           ion2 = self.get_type_i(5)
           return np.any(np.array([ion1,ion2]),axis=0)
        elif self.sol_ty == "benz":
           return self.get_benzene()

        return [False]*len(self.types)

    def get_benzene(self):
        '''Get a list len(nats) that is true just for the first C atom of each
           benzene ring'''
        benz_c = self.get_type_i(4); c_ct = 0
        for i in range(len(benz_c)):
            if benz_c[i] == True: 
                c_ct += 1
                if c_ct % 6 != 0: benz_c[i] = False
        return benz_c

    def get_spacing_for_interlayer(self):
        '''Depending on the spacing of the graphene walls, return array of
           positions for computing layered properties on, e.g. g_{2D}(R) or
           3B angles '''
        nm = np.linspace(0,1000.,0.6)

        # depending on the system, there will be varying # interlayers
        if "_6_" in self.xyzfname or "_7_" in self.xyzfname or "_8_" in self.xyzfname:
             nm = np.array([0.0,10.0])
        elif ("_8." in self.xyzfname or "_9" in self.xyzfname or 
              "_10_" in self.xyzfname or "_11_" in self.xyzfname
               or "_11." in self.xyzfname 
               or "_12_" in self.xyzfname):
             nm = np.array([0,0.55,1.04,3.0,3.30,3.90,5.1,5.6,6.7,10.0])
        elif ("_12." in self.xyzfname or "_13_" in self.xyzfname or "_14_" in self.xyzfname):
             nm = np.array([0,0.35,1.24,3.24,3.90,5.9,6.7,10.0])
        elif ("_16_" in self.xyzfname):#nm = 7
             nm = np.array([0,0.55,1.24,3.24,3.70,5.4,6.1,8.1,8.7,14.0])
        elif ("_20_" in self.xyzfname):
             nm = np.array([0,0.65,1.04,3.14,3.6,5.3,5.9,7.05,7.45,9.3,9.8,12.05,12.6,18.])
        elif ("_37_" in self.xyzfname):
             nm = np.linspace(0,38,38)

        return nm

def main():
    '''This is mostly for testing'''
    filename=sys.argv[1]; volf=sys.argv[2]
    vC = VolFile(volf)
    zf = XYZFile(filename, vC)

if __name__=="__main__":
    main()

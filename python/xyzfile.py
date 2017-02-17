import numpy as np

from volfile import VolFile

WOXY = 1; WHYD = 2; GRAPHENE = 3

class XYZFile:
     '''A class for xyz files'''
     xyzfname = ''
         
     def __init__(self, fname, VolFile):
         self.xyzfname = fname
         if VolFile.volfname == '': self.get_coords_types(fname, [], [])
         else: 
             self.get_coords_types(fname, VolFile.time, VolFile.dims)
             self.half_x = VolFile.get_x_max()/2.0

     def get_coords_types(self, filename, times, dims):
         '''Method to open xyz file and save coords and types'''
         time, atoms, atom, types = [], [], [], []                                      
         grab_snap, t_ct = 0, 0                                                         
         f=open(filename, "r")                                                          
         f.readline()                                                                   
         line = f.readline()                                                            
         if times == [] or times[0][0] == int(line.split()[-1]):                        
             time.append(int(line.split()[-1]))                                         
             t_ct += 1; grab_snap = 1                                                   
                                                                                        
         '''First, we read in the coordinates from the trajectory file and              
            save them into atom=[[coords],[coords]...,[coords]]'''                      
         for line in f:                                                                 
             if times != [] and t_ct >= len(times): break                               
                                                                                        
             #Gets past # and atoms lines                                               
             if line[0]=="A":                                                           
                 if grab_snap == 1:                                                     
                     atom+=[atoms]                                                      
                     atoms = []                                                         
                                                                                        
                 if times == [] or times[t_ct][0] == int(line.split()[-1]):             
                     time.append(int(line.split()[-1]))                                 
                     t_ct += 1; grab_snap = 1                                           
                 else: grab_snap = 0                                                    
             # Only grab snapshots that have volume data                                
             elif len(line.split()) > 1 and grab_snap == 1:                             
                 #crds is a list of the coordinates in string formats                   
                 crds=line.split()                                                      
                 if t_ct == 1: types.append(int(crds[0]))                               
                 if dims == []:                                                         
                     tmp = []
                     for i in range(1,len(crds)): tmp.append(float(crds[i]))
                     atoms.append(tmp)
                 else: # Move coords so that box is from 0 -> xmin                      
                     atoms.append([float(crds[1])-dims[t_ct][0],                        
                                   float(crds[2])-dims[t_ct][2],                        
                                   float(crds[3])-dims[t_ct][4]])                       
         #For last timestep, since it doesn't have a line about Atoms after it.         
         if atoms != []: atom += [atoms]
         self.time  = np.array(time); self.atom = np.array(atom)                                                         
         self.types = np.array(types)                                                        
         print("xyz", filename, len(times), self.atom.shape, self.types.shape)                                      
         f.close()                      
   

     def get_type_i(self, i):
         '''Get indices of type i'''
         return self.types == i
     
     def get_graph_wall(self, wall_no=0):
         '''Get the indices of graphene wall 0 or 1'''
         grap = self.types == GRAPHENE
     
         if wall_no == 0:
             xl = self.atom[0,:,0] < self.half_x
             return np.all(np.array([grap, xl]), axis=0) # first wall
         else: 
             xg = self.atom[0,:,0] > self.half_x
             return np.all(np.array([grap, xgr]), axis=0) # second wall
     
     def get_wall_i_xv(self, i=0):
         '''Get the x location of graphene wall 0 or 1, from first snap'''
         wall_idxs = self.get_graph_wall(i)
         return np.mean(self.atom[0,wall_idxs,0])
     
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




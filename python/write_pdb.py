import sys
import numpy as np

from get_xyz import xyz_reader
from get_box_dim import get_box_dim_from_vol

typ_dic = { 
            1: ['O WAT', -1.0484,], 
            2: ['H WAT',  0.5242,], 
            3: ['C GRA',  0.0,], 
            4: ['C BEN', -0.115,], 
            5: ['H BEN',  0.115,],
          }

def write_to_pdb(fname, dims, types, coords):
    '''Given coordinates, PBC box dimensions, and type list,
       write out a pdb trajectory file'''
    f = open(fname, 'w')
    for tim in range(len(coords)):
        res_ct = 0
        for i in range(len(types)):
            if types[i] == 1: res_ct += 1
            elif types[i] == 3 and types[i-1] == 2: res_ct += 1
            elif types[i] == 4 and types[i-1] != types[i]: res_ct += 1

            f.write("ATOM{0:7d}    {1:s} {2:5d}    {3:8.3f}{4:8.3f}{5:8.3f} {ll:5.2f}  0.00\n".format(
                     i, typ_dic[types[i]][0], res_ct, #coords[tim,i,0], coords[tim,i,1], coords[tim,i,2],\
                     *coords[tim, i, :], 
                     ll=typ_dic[types[i]][1]))
        f.write("END\n")
    f.close()

def main():
    xyzname = sys.argv[1]
    sep = sys.argv[2]
    itr = sys.argv[3]
    time, dims = get_box_dim_from_vol("run"+str(sep)+"_"+str(itr)+".vol")
    tim_xyz, coords, types = xyz_reader(xyzname, time, dims)
    coords = np.array(coords); types = np.array(types); dims = np.array(dims)

    write_to_pdb(xyzname[:-3]+"pdb", dims, types, coords)

if __name__=="__main__":
    main()

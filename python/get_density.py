import sys
import numpy as np


GRAPHENE = 3

type_wt = {
            1: 16.00,
            2:  1.01,
            3: 12.00,
            4: 12.00,
            5:  1.01,
          }


def histog_dist():
    ''' Get a histogram of distances and transform into distance from plate'''


def main():                                                                        
    xyzname = sys.argv[1]                                                          
    sep = sys.argv[2]                                                              
    itr = sys.argv[3]                                                              
    time, dims = get_box_dim_from_vol("run"+str(sep)+"_"+str(itr)+".vol")
    tim_xyz, coords_ar, types_ar = xyz_reader(xyzname, time, dims)
    dims_ar = np.array(dims)       

    histog_dist()


                                                                                   
if __name__=="__main__":                                                           
    main() 

import sys
import numpy as np

fname = sys.argv[1]

f = open(fname,"r").readlines()
fo = open("angle_g_bin_ct.csv", "w")
fo.write("bin,ct\n")

RMAX = 9.0; binsiz_ra = 0.1
bns_ra = np.arange(0.0,RMAX,binsiz_ra)

for l in f:
    if "bins" in l:
        tmp = l.split()
        for i in range(1, len(tmp)): 
            fo.write("{0:.4f},{1}\n".format(bns_ra[i-1],tmp[i]))

fo.close()

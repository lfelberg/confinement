import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cmx 
from pylab import *

from xyzfile import XYZFile
from volfile import VolFile

colorL = [[0,0,0], [0,0,1], [1, 0,0], [.93, .53, .18]]

def plot_scatter(plt_nm, csvL, loc):
    '''Using data from a histogram, plot several'''
   #mn = np.mean(np.mean(coords[:,:,2]))
   #rng = np.std(coords[:,:,2])
   #xrn = [np.amin(coords[:,:,0]), np.amax(coords[:,:,0])]
   #yrn = [np.amin(coords[:,:,1]), np.amax(coords[:,:,1])]
   #print(xrn, yrn)
    f = plt.figure(1, figsize = (3.0, 3.0))
    ax, ct = f.add_subplot(111), 0

    for i in range(len(csvL)):
        for j in range(len(csvL[i])):
            for k in range(len(csvL[i][k])):
                shft = max(csvL[i][j][k].dat[0])/2.0
                ax.scatter(csvL[i][j][k].dat[0]-shft, 
                           csvL[i][j][k].dat[loc], color=colorL[ct])
                ct += 1
    ax.set_xlim(xrn); ax.set_ylim(yrn)
    plt.savefig(plt_nm+str(j)+'.png', format='png',                       
                    bbox_inches = 'tight', dpi=300) 
    plt.close()

def main():                                                                        
    '''For a collection of data, get info from csv and then plot'''
    csvname = sys.argv[1]; nsep = int(sys.argv[2]); nlen = int(sys.argv[3]);
    niter = int(sys.argv[4]); st_other = 5
    sep, ln, itr = [], [], []
    print("This is nsep {0}, nlen {1}, niter {2}".format(nsep, nlen, niter))

    for i in range(st_other,nsep+st_other+1):
       sep.append(int(sys.argv[i]))

    for i in range(nsep+st_other+1,nsep+st_other+nlen+1):
       ln.append(int(sys.argv[i]))
    
    for i in range(nsep+st_other+nlen+1,nsep+st_other+niter+nlen+1):
       itr.append(int(sys.argv[i]))
    print(sep,ln, itr)
    
    ext = sys.argv[st_other+nsep+nlen+niter+1]
    loc = sys.argv[st_other+nsep+nlen+niter+2]
    print("Ext {0}, loc {1}".format(ext, loc))
    csvL = []
    for i in range(nsep):
        csL = []
        for j in range(nlen):
            cs = []
            for k in range(niter):
                cs.append(CSVFile(fname+str(sep[i])+"_"+str(ln[j])+"_"
                                  +str(itr[k])+ext))
            csL.append(cs)
        csvL.append(csL)
    plot_scatter(csvname[:-3], csvL, loc)

if __name__=="__main__":
    main()

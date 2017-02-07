import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cmx 
from pylab import *

from xyzfile import XYZFile
from volfile import VolFile

color = [[0,0,0], [0,0,1], [1, 0,0], [.93, .53, .18]]

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
            shft = max(csvL[i][j].dat[0])/2.0
            ax.scatter(csvL[i][j].dat[0]-shft, csvL[i][j].dat[loc], color[ct])
            ct += 1
    ax.set_xlim(xrn); ax.set_ylim(yrn)
    plt.savefig(plt_nm+str(j)+'.png', format='png',                       
                    bbox_inches = 'tight', dpi=300) 
    plt.close()

def main():                                                                        
    '''For a collection of data, get info from csv and then plot'''
    csvname = sys.argv[1]
    nsep = int(sys.argv[2])
    niter = int(sys.argv[3])
    sep, itr = [], [] 

    for i in range(4,nsep+1):
       sep.append(int(sys.argv[i]))

    for i in range(nsep+1,nsep+niter+1):
       itr.append(int(sys.argv[i]))
    
    ext = sys.argv[nsep+niter+1]
    loc = sys.argv[nsep+niter+2]
    csvL = []
    for i in range(nsep):
        cs = []
        for j in range(niter):
            cs.append(CSVFile(fname+str(sep[i])+"_"+str(itr[j])+ext))
        csvL.append(cs)
    plot_scatter(csvname[:-3], csvL, loc)
                                                                                   
if __name__=="__main__":                                                           
    main()

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cmx 
from pylab import *

from xyzfile import XYZFile

def plot_graph_dev(plt_nm, coords):
    '''For each time snapshot, plot a 2D heatmap of deviation of x dim'''
    mn = np.mean(coords[0,:,2])
    rng = np.std(coords[0,:,2])
    for j in range(len(coords)):
        print(np.mean(coords[j,:,2]), np.std(coords[j,:,2]))
        f = plt.figure(1, figsize = (3.0, 3.0))
        ax = f.add_subplot(111)
        xran = [min(coords[j,:,0]), max(coords[j,:,0])]
        yran = [min(coords[j,:,1]), max(coords[j,:,1])]
        cNorm = matplotlib.colors.Normalize(vmin=mn-5.*rng, vmax=mn+5.*rng)                                  
        cm = plt.get_cmap('seismic_r')                                                 
        scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)   
        ax.scatter(coords[j,:,0], coords[j,:,1], s = 45, lw = 0,
                   c=scalarMap.to_rgba(coords[j,:,2]))
        ax.set_xlim(xran); ax.set_ylim(yran)
        scalarMap.set_array(coords[j,:,2])
        cbaxes = f.add_axes([0.96, 0.11, 0.025, 0.75])
        cb = plt.colorbar(scalarMap, cax = cbaxes)                                     
       #cbaxes.set_xlabel(units, fontname='Arial',                                     
       #                  fontsize='small', labelpad=5.) 
        plt.savefig(plt_nm+str(j)+'.png', format='png',                       
                        bbox_inches = 'tight', dpi=300) 
        plt.close()

def main():                                                                        
    xyzname = sys.argv[1]; #sep = sys.argv[2]; itr = sys.argv[3]
    xyzF = XYZFile(xyzname)
    plot_graph_dev(xyzname[:-3], xyzF.atom)
                                                                                   
                                                                                   
if __name__=="__main__":                                                           
    main()

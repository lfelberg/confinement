import sys
import numpy as np
import matplotlib, scipy
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema
from matplotlib.ticker import *
MaxNLocator.default_params['nbins']=4

from csvfile import CSVFile

colorL = [[0,0,0], [0,0,1], [1, 0,0], [.93, .53, .18]]

def plot_scatter(csv):
    '''Using data from a histogram, plot several'''
    f = plt.figure(1, figsize = (1.5, 1.5)); ax = f.add_subplot(111)
    plt.rcParams.update({'font.size': 8}); flx,df = [],[];fl,ff = [],[]
    ax.xaxis.set_major_locator(MaxNLocator())
    ax.yaxis.set_major_locator(MaxNLocator())

    dens  = csv.dat[csv.find_keyword('dens')].flatten(); 
    rigid = csv.dat[csv.find_keyword('dgg_r')].flatten(); 
    flexb = csv.dat[csv.find_keyword('dgg_f')]
    for j in range(flexb.shape[1]): 
        if flexb[1][j] != 0.:
            ff.append(dens[j]); fl.append(flexb[0][j])
            ff.append(dens[j]); fl.append(flexb[1][j])
        else:  df.append(dens[j]); flx.append(flexb[0][j])
    plt.plot(ff, fl, 'mD', label = "multi")
    plt.plot(df, flx, 'k.', label = "flexible")
    plt.plot(dens, rigid, 'b.',label = "rigid")
    ax.set_xlim([0.,0.35]);#ax.set_ylim([0,0.02])

    ax.legend(loc = 2, ncol = 1, columnspacing = 0.4,fontsize=7)
    ax.set_xlabel(r"$\rho_{2D}$",fontsize= 10)
    ax.set_ylabel(r"$\langle d_{gg} \rangle \,\, (\AA)$",fontsize= 10)
    plt.savefig(csv.csvfname[:-4]+'_r.png',bbox_inches = 'tight',) #transparent=True)
    plt.close()

def main():                                                                        
    '''For a collection of data, get info from csv and then plot,
       usage: plot_vs_x.py csvStart nsep nlen niter sep1 sep2... len1 len2... 
                           iter1 iter2... ext datLoc'''
    csvname = sys.argv[1]
    plot_scatter(CSVFile(csvname))

if __name__=="__main__":
    main()

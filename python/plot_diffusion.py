import sys
import numpy as np
import matplotlib, scipy
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema
from matplotlib.ticker import *
MaxNLocator.default_params['nbins']=4

from csvfile import CSVFile
from gauss_fits import skew, double, normal

colorL = [[0,0,0], [0,0,1], [1, 0,0], [.93, .53, .18]]

def plot_scatter(csv):
    '''Using data from a histogram, plot several'''
    nbins = 150
    f = plt.figure(1, figsize = (1.5, 1.5)); ax = f.add_subplot(111)
    plt.rcParams.update({'font.size': 8})
    ax.xaxis.set_major_locator(MaxNLocator())
    ax.yaxis.set_major_locator(MaxNLocator())

    dens  = csv.dat[csv.find_keyword('dens')].flatten(); 
    rigid = csv.dat[csv.find_keyword('dr_m2_s')].flatten(); 
    flexb = csv.dat[csv.find_keyword('df_m2_s')].flatten(); 
   #plt.plot(dens, flexb, 'k.--', label = "flexible")
    plt.plot(dens, rigid, 'b.--',label = "rigid")
    ax.set_xlim([0.,0.35]); ax.set_ylim([0,6.02])
   #ax.set_yscale("log", nonposy='clip')

   #ax.legend(loc = 0, ncol = 1, fontsize=7, columnspacing = 0.4,)
    ax.set_xlabel(r"$\rho_{2D}$",fontsize= 10)
    ax.set_ylabel(r"$\mathcal{D}_{||} \times 10^9\, (m^2/s)$",fontsize= 10)
    plt.savefig(csv.csvfname[:-4]+'_r.png',bbox_inches = 'tight',transparent=True)
    plt.close()

def main():                                                                        
    '''For a collection of data, get info from csv and then plot,
       usage: plot_vs_x.py csvStart nsep nlen niter sep1 sep2... len1 len2... 
                           iter1 iter2... ext datLoc'''
    csvname = sys.argv[1]
    plot_scatter(CSVFile(csvname))

if __name__=="__main__":
    main()

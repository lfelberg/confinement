import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import *
MaxNLocator.default_params['nbins']=3

from csvfile import CSVFile
from gauss_fits import skew, double, normal

colorL = [[0,0,0], [0,0,1], [1,0,0], [.93, .53, .18]]

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
    plt.plot(dens,flexb,'k.--',dashes = (2,1),label = "flexible")
    plt.plot(dens,rigid,'b.--',dashes = (2,1),label = "rigid")
    # plot calculated tip4p EW D_{||} @ 298 K
    x = np.linspace(0,1,70); y = np.ones(70)*2.4725
    plt.plot(x,y, 'r--',dashes=(1.5,0.9),linewidth=0.7,label="bulk,\n 298K")

    ax.set_xlim([0.,0.40]); ax.set_ylim([0,4])
    ax.legend(ncol = 3, fontsize=6, columnspacing = 0.2,
              bbox_to_anchor = (1.05, 1.25), borderaxespad= 0.2)
    ax.set_xlabel(r"$\rho_{2D}$",fontsize= 10)
    ax.set_ylabel(r"$\mathcal{D}_{||} \times 10^9\, (m^2/s)$",fontsize= 10)
    plt.savefig(csv.csvfname[:-4]+'.png',bbox_inches = 'tight')#,transparent=True)
    plt.close()

def main():                                                                        
    '''For a collection of data, get info from csv and then plot,
       usage: plot_vs_x.py csvStart nsep nlen niter sep1 sep2... len1 len2... 
                           iter1 iter2... ext datLoc'''
    csvname = sys.argv[1]
    plot_scatter(CSVFile(csvname))

if __name__=="__main__":
    main()

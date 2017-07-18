import sys
import numpy as np
import matplotlib, scipy
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema
from matplotlib.ticker import *
MaxNLocator.default_params['nbins']=5

from csvfile import CSVFile
from gauss_fits import skew, double, normal

colorL = [[0,0,0], [0,0,1], [1, 0,0], [.93, .53, .18]]

def plot_scatter(csv, sep, ln, itr):
    '''Using data from a histogram, plot several'''
    nbins = 50
    f = plt.figure(1, figsize = (1.5, 1.5)); ax = f.add_subplot(111)
    ax.xaxis.set_major_locator(MaxNLocator())
    ax.yaxis.set_major_locator(MaxNLocator())

    angles = csv.dat[csv.find_keyword('the1')].flatten(); rng = [0, 180]
    y,x = np.histogram(angles, bins=nbins, range=rng); dx = x[1]-x[0]
    x = x[:-1] + dx/2.;y = y/sum(y.astype(float))/dx
    plt.plot(x, y, color = "k")
    ax.set_xlim(rng); ax.set_ylim([0,0.02])

    ax.set_xlabel(r"$\theta$",fontsize= 12)
    ax.set_ylabel(r"PDF",fontsize= 12)
    plt.savefig(csv.csvfname[:-4]+'.png',bbox_inches = 'tight',)

    plt.plot((90.,90.), (0,10), 'y-')
    plt.plot((109.,109.), (0,10), 'r-')
    plt.plot((160.,160.), (0,10), 'b-')
    params = [90.0, 120., 1.0, 1.0, 1.0, 1.0]
    print("Angle maxes: ", x[argrelextrema(y, np.greater)])
   #plt.savefig(csv.csvfname[:-4]+'_fit.png',bbox_inches='tight'); plt.close()
    plt.close()

    ft_fl = open(csv.csvfname[:-11]+'_hist.csv', 'w')
    ft_fl.write("bin,hist\n")
    for xx in range(len(x)): 
        ft_fl.write("{0:.3f},{1:.4f}\n".format(x[xx],y[xx]))
    ft_fl.close()


def main():                                                                        
    '''For a collection of data, get info from csv and then plot,
       usage: plot_vs_x.py csvStart nsep nlen niter sep1 sep2... len1 len2... 
                           iter1 iter2... ext datLoc'''
    csvname = sys.argv[1]; sep = float(sys.argv[2]); ln = int(sys.argv[3]);
    itr = sys.argv[4]
    
    plot_scatter(CSVFile(csvname), sep, ln, itr)

if __name__=="__main__":
    main()

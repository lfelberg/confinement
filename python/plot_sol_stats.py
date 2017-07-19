import sys
import numpy as np
import matplotlib, scipy
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema
from matplotlib.ticker import *
MaxNLocator.default_params['nbins']=5

from csvfile import CSVFile

def plot_scatter(csv, sep, ln, itr):
    '''Using data from a histogram, plot several'''
    leg = str(sep)+r'$\AA$ Sep, '+'L='+str(ln)+r'$\AA$'
    neigh = csv.dat[csv.find_keyword("neigh")]
    ang  = csv.dat[csv.find_keyword("ang")].flatten()
    dist = csv.dat[csv.find_keyword("dist")].flatten()

   #print(min(neigh), max(neigh), min(ang), max(ang))

    # finding solute-water oxygen within cutoff
    f = plt.figure(1, figsize = (1.0, 1.0))
    ax, ct, leg = f.add_subplot(111), 0, []
    ax.yaxis.set_ticklabels([])
    ax.xaxis.set_major_locator(MaxNLocator())
    ax.yaxis.set_major_locator(MaxNLocator())
    matplotlib.rcParams['font.size'] = 5;
    nbins = int(neigh.max())+1; ra = (-0.5, nbins-0.5)
    vac, surr = 0., 0.
    for i in range(len(neigh)):
        y,x = np.histogram(neigh[i], bins=nbins, range=ra); dx = x[1]-x[0]
        neig_m = np.mean(x[argrelextrema(y, np.greater)])
        print("sol {0}: neigh {1:.2f}".format(i, neig_m))
        if neig_m > 18: surr += 1.
        else:           vac  += 1.
    print("{0},{1},{2},{3:.6f}".format(csv.csvfname,surr,vac,vac/(surr+vac)))
        
    y,x = np.histogram(neigh.flatten(), bins=nbins, range=ra); dx = x[1]-x[0]
    x = x[:-1] + dx/2.;y = y/sum(y.astype(float))/dx
    plt.plot(x, y, color = "k")
    ax.set_xlabel(r"$O_W < 10 \AA$",fontsize=7); ax.set_xlim([0,25])
    ax.set_ylabel("Probability",fontsize=7);#ax.set_ylim([0.,2.0])
    plt.savefig(csv.csvfname[:-3]+'neigh.png',bbox_inches = 'tight',)
    plt.close()

    nbins = 25
    f = plt.figure(1, figsize = (1.0, 1.0))
    ax, ct, leg = f.add_subplot(111), 0, []
    ax.yaxis.set_ticklabels([])
    ax.xaxis.set_major_locator(MaxNLocator())
    ax.yaxis.set_major_locator(MaxNLocator())
    y,x = np.histogram(ang, bins=nbins); dx = x[1]-x[0]
    x = x[:-1] + dx/2.;y = y/sum(y.astype(float))/dx
    print("Ang position maxes: ", x[argrelextrema(y, np.greater)])
    plt.plot(x, y, color = "k")
    ax.set_xlabel(r"$\theta_B$",fontsize=7); ax.set_xlim([0,180])
    ax.set_ylabel("Probability",fontsize=7);#ax.set_ylim([0.,2.0])
    plt.savefig(csv.csvfname[:-3]+'ang.png',bbox_inches = 'tight',)
    plt.close()

    nbins = 30
    y,x = np.histogram(dist, bins=nbins); dx = x[1]-x[0]
    x = x[:-1] + dx/2.;y = y/sum(y.astype(float))/dx
    print("dgg position maxes: ", x[argrelextrema(y, np.greater)])
    f = plt.figure(1, figsize = (1.0, 1.0))
    ax, ct, leg = f.add_subplot(111), 0, []
    ax.yaxis.set_ticklabels([])
    ax.xaxis.set_major_locator(MaxNLocator())
    ax.yaxis.set_major_locator(MaxNLocator())
    plt.plot(x, y, color = "k")
    ax.set_xlabel(r"$d_{{gg}}^B$",fontsize=7); ax.set_xlim([3,18])
    ax.set_ylabel("Probability",fontsize=7);#ax.set_ylim([0.,2.0])
    plt.savefig(csv.csvfname[:-3]+'dist.png',bbox_inches = 'tight',)
    plt.close()

def main():                                                                        
    '''For a collection of data, get info from csv and then plot,
       usage: plot_vs_x.py csvStart nsep nlen niter sep1 sep2... len1 len2... 
                           iter1 iter2... ext datLoc'''
    csvname = sys.argv[1]; sep = float(sys.argv[2]); ln = int(sys.argv[3]);
    itr = sys.argv[4]
    plot_scatter(CSVFile(csvname), sep, ln, itr)

if __name__=="__main__":
    main()

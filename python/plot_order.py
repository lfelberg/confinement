import sys
import numpy as np
import matplotlib, scipy
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema
from matplotlib.ticker import *
MaxNLocator.default_params['nbins']=5

from csvfile import CSVFile

colorL = [[0,0,0], [0,0,1], [1, 0,0], [.93, .53, .18]]

def plot_scatter(csv, sep, ln, itr):
    '''Using data from a histogram, plot several'''
    nbins = 200
    f = plt.figure(1, figsize = (1.0, 1.0))
    ax, ct, leg = f.add_subplot(111), 0, []
    ax.xaxis.set_major_locator(MaxNLocator())
    ax.yaxis.set_major_locator(MaxNLocator())
    
    leg = str(sep)+r'$\AA$ Sep, '+'L='+str(ln)+r'$\AA$'
    q4b = csv.dat[csv.key.index("q4_bar")]
    q6b = csv.dat[csv.key.index("q6_bar")]
    print(min(q4b), max(q4b), min(q6b), max(q6b))
    plt.hist2d(q4b, q6b, bins=(nbins,nbins), range=([0,1.0], [0,1.0]),cmap=plt.get_cmap('plasma'))
   #ax.set_xlim([0.2,0.8]); ax.set_ylim(yr)

    ax.set_xlabel("$\\bar{q}_{4}$",fontsize= 8)
    ax.set_ylabel("$\\bar{q}_{6}$",fontsize= 8)
    plt.savefig(csv.csvfname[:-3]+'.png',bbox_inches = 'tight',)
    plt.close()

    # finding graphene separation dist(s)
    nbins = 30
    y,x = np.histogram(q4b, bins=nbins); dx = x[1]-x[0]
    x = x[:-1] + dx/2.;y = y/sum(y.astype(float))/dx
    print("X position maxes: ", x[argrelextrema(y, np.greater)])
    f = plt.figure(1, figsize = (1.0, 1.0))
    ax, ct, leg = f.add_subplot(111), 0, []
    ax.xaxis.set_major_locator(MaxNLocator())
    ax.yaxis.set_major_locator(MaxNLocator())
    matplotlib.rcParams['font.size'] = 5;
    plt.plot(x, y, color = "k")
    ax.set_xlabel("$\\bar{q}_{4}$",fontsize=7); ax.set_xlim([0,1.0])
    ax.set_ylabel("Probability",fontsize=7);#ax.set_ylim([0.,2.0])
    plt.savefig(csv.csvfname[:-3]+'q4b.png',bbox_inches = 'tight',)
    plt.close()

    y,x = np.histogram(q6b, bins=nbins); dx = x[1]-x[0]
    x = x[:-1] + dx/2.;y = y/sum(y.astype(float))/dx
    print("X position maxes: ", x[argrelextrema(y, np.greater)])
    f = plt.figure(1, figsize = (1.0, 1.0))
    ax, ct, leg = f.add_subplot(111), 0, []
    ax.xaxis.set_major_locator(MaxNLocator())
    ax.yaxis.set_major_locator(MaxNLocator())
    plt.plot(x, y, color = "k")
    ax.set_xlabel("$\\bar{q}_{6}$",fontsize=7); ax.set_xlim([0,1.0])
    ax.set_ylabel("Probability",fontsize=7);#ax.set_ylim([0.,2.0])
    plt.savefig(csv.csvfname[:-3]+'q6b.png',bbox_inches = 'tight',)
    plt.close()

def main():                                                                        
    '''For a collection of data, get info from csv and then plot,
       usage: plot_vs_x.py csvStart nsep nlen niter sep1 sep2... len1 len2... 
                           iter1 iter2... ext datLoc'''
    csvname = sys.argv[1]; sep = int(sys.argv[2]); ln = int(sys.argv[3]);
    itr = sys.argv[4]
    
    plot_scatter(CSVFile(csvname), sep, ln, itr)

if __name__=="__main__":
    main()

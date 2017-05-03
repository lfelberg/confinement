import sys
import numpy as np
import matplotlib, scipy
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema
from matplotlib.ticker import *
MaxNLocator.default_params['nbins']=5

from csvfile import CSVFile
from gauss_fits import skew, double, triple

colorL = [[0,0,0], [0,0,1], [1, 0,0], [.93, .53, .18]]

def plot_scatter(csv, sep, ln, itr):
    '''Using data from a histogram, plot several'''
    nbins = 200
    f = plt.figure(1, figsize = (1.0, 1.0))
    ax, ct, leg = f.add_subplot(111), 0, []
    ax.xaxis.set_major_locator(MaxNLocator())
    ax.yaxis.set_major_locator(MaxNLocator())
    
    leg = str(sep)+r'$\AA$ Sep, '+'L='+str(ln)+r'$\AA$'
    Y = csv.dat[csv.key.index("dgg")]; yr = [5,17]
    yy = Y < sep + float(sep)*0.3; Y = Y[yy]
    X = csv.dat[csv.key.index("dg1"),yy]/Y.astype(float)
    print(min(Y), max(Y), min(X), max(X))
    plt.hist2d(X, Y, bins=(nbins,nbins), range=([0,1], yr),cmap=plt.get_cmap('plasma'))
    ax.set_xlim([0.2,0.8]); ax.set_ylim(yr)

    ax.set_xlabel("$x/d_{gg}$",fontsize= 8)
    ax.set_ylabel("$d_{gg} \, (\AA)$",fontsize= 8)
    plt.savefig(csv.csvfname[:-3]+'.png',bbox_inches = 'tight',)
    plt.close()

    # finding graphene separation dist(s)
    bil = [9,12] #list of initial seps that form 2 diff layers in flexible
    nbins = 50
    y,x = np.histogram(Y, bins=nbins); dx = x[1]-x[0]
    x = x[:-1] + dx/2.;y = y/sum(y.astype(float))/dx
    f = plt.figure(1, figsize = (1.0, 1.0))
    ax, ct, leg = f.add_subplot(111), 0, []
    ax.xaxis.set_major_locator(MaxNLocator())
    ax.yaxis.set_major_locator(MaxNLocator())
    matplotlib.rcParams['font.size'] = 5;
    plt.plot(x, y, color = "k")
    ax.set_xlabel("$d_{gg} \, (\AA$)",fontsize=7);ax.set_xlim([6.,17.])
    ax.set_ylabel("Probability $(1/\AA)$",fontsize=7);ax.set_ylim([0.,2.0])
    plt.savefig(csv.csvfname[:-3]+'g_sep_fit.png',bbox_inches = 'tight',)
    plt.close()

   #ft_fl = open(csv.csvfname[:-3]+'g_sep_fit.dat', 'w')
   #ft_fl.write("bin,dat,fit\n")
   #for xx in range(len(x)): 
   #    ft_fl.write("{0:.3f},{1:.4f},{2:.4f}\n".format(x[xx],y[xx],fit[xx]))
   #ft_fl.close()

    bil = [9,12] #list of initial seps that form 2 diff layers in flexible
    y,x = np.histogram(X, bins=nbins); dx = x[1]-x[0]
    x = x[:-1] + dx/2.;y = y/sum(y.astype(float))/dx
    print("X position maxes: ", x[argrelextrema(y, np.greater)])
    f = plt.figure(1, figsize = (1.0, 1.0))
    ax, ct, leg = f.add_subplot(111), 0, []
    ax.xaxis.set_major_locator(MaxNLocator())
    ax.yaxis.set_major_locator(MaxNLocator())
    plt.plot(x, y, color = "k")
    ax.set_xlabel("$x/d_{gg}$",fontsize=7); ax.set_xlim([0.2,0.8])
    ax.set_ylabel("Probability",fontsize=7); ax.set_ylim([0.,15])
    plt.savefig(csv.csvfname[:-3]+'o_sep_fit.png',bbox_inches = 'tight',)
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

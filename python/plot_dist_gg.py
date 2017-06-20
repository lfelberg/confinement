import sys
import numpy as np
import matplotlib, scipy
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema
from matplotlib.ticker import *
MaxNLocator.default_params['nbins']=5

from csvfile import CSVFile
from gauss_fits import skew, double, triple

def plot_scatter(csv, sep, ln, itr):
    '''Using data from a histogram, plot several'''
    nbins = 200
    f = plt.figure(1, figsize = (1.0, 1.0))
    ax, ct, leg = f.add_subplot(111), 0, []
    ax.xaxis.set_major_locator(MaxNLocator())
    ax.yaxis.set_major_locator(MaxNLocator())
    
    leg = str(sep)+r'$\AA$ Sep, '+'L='+str(ln)+r'$\AA$'
    Y = csv.dat[csv.key.index("dgg")]; yr = [5,17] #REST #[10,20] 
    yy = np.all(np.array([Y > 5.8, Y < sep + float(sep)*0.3]),axis =0); Y=Y[yy]
    X = csv.dat[csv.key.index("dg1"),yy]/Y.astype(float)
    xx = X > 1.0; X[xx] = X[xx] - 1
    xx = X > 1.0; X[xx] = X[xx] - 1
    print(min(Y), max(Y), min(X), max(X))
    plt.hist2d(X, Y, bins=(nbins,nbins), range=([0,1], yr),cmap=plt.get_cmap('plasma'))
    ax.set_xlim([0.2,0.8]); ax.set_ylim(yr)

    ax.set_xlabel("$x/d_{gg}$",fontsize= 8)
    ax.set_ylabel("$d_{gg} \, (\AA)$",fontsize= 8)
    plt.savefig(csv.csvfname[:-3]+'.png',bbox_inches = 'tight',)
    plt.close()

    # finding graphene separation dist(s)
    nbins = 50
    y,x = np.histogram(Y, bins=nbins); dx = x[1]-x[0]
    x = x[:-1] + dx/2.;y = y/sum(y.astype(float))/dx
    print("Dgg maxes: ", x[argrelextrema(y, np.greater)])
    f = plt.figure(1, figsize = (1.0, 1.0))
    ax, ct, leg = f.add_subplot(111), 0, []
    ax.xaxis.set_major_locator(MaxNLocator())
    ax.yaxis.set_major_locator(MaxNLocator())
    matplotlib.rcParams['font.size'] = 5;
    plt.plot(x, y, color = "k")
    plt.plot(np.ones(70)*10.5, np.linspace(0,2,70), color = "r")
    ax.set_xlabel("$d_{gg} \, (\AA$)",fontsize=7);ax.set_xlim([6.,20.])
    ax.set_ylabel("Probability $(1/\AA)$",fontsize=7);ax.set_ylim([0.,2.0])
    plt.savefig(csv.csvfname[:-3]+'g_sep_fit.png',bbox_inches = 'tight',)
    plt.close()

   #ml = Y < 10.5
   #bl = Y > 10.5
   #Xml = X[ml]#/Y[ml].astype(float)
   #Xbl = X[bl]#/Y[bl].astype(float)
   #y_ml,x_ml = np.histogram(Xml, bins=nbins); dx = x_ml[1]-x_ml[0]
   #x_ml = x_ml[:-1] + dx/2.;y_ml = y_ml/sum(y_ml.astype(float))/dx

   #f = open("run121_46_1.diso_sep_fit.dat", "w")
   #f.write("x_dgg,x_val\n")
   #for i in range(len(x_ml)): f.write("{0:.5f},{1:.5f}\n".format(x_ml[i],y_ml[i]))
   #f.close()

   #y_bl,x_bl = np.histogram(Xbl, bins=nbins); dx = x_bl[1]-x_bl[0]
   #x_bl = x_bl[:-1] + dx/2.;y_bl = y_bl/sum(y_bl.astype(float))/dx

   #f = open("run122_46_1.diso_sep_fit.dat", "w")
   #f.write("x_dgg,x_val\n")
   #for i in range(len(x_bl)): f.write("{0:.5f},{1:.5f}\n".format(x_bl[i],y_bl[i]))
   #f.close()

    y,x = np.histogram(X, bins=nbins); dx = x[1]-x[0]
    x = x[:-1] + dx/2.;y = y/sum(y.astype(float))/dx
    print("X position maxes: ", x[argrelextrema(y, np.greater)])
    f = plt.figure(1, figsize = (1.0, 1.0))
    ax, ct, leg = f.add_subplot(111), 0, []
    ax.xaxis.set_major_locator(MaxNLocator())
    ax.yaxis.set_major_locator(MaxNLocator())
    plt.plot(x, y, color = "k")
    ax.set_xlabel("$x/d_{gg}$",fontsize=7); ax.set_xlim([0.2,0.8])
    ax.set_ylabel("Probability",fontsize=7);#ax.set_ylim([0.,15])
    plt.savefig(csv.csvfname[:-3]+'o_sep_fit.png',bbox_inches = 'tight',)
    plt.close()
    f = open(csv.csvfname[:-3]+"o_sep_fit.dat", "w")
    f.write("x_dgg,x_val\n")
    for i in range(len(x)): f.write("{0:.5f},{1:.5f}\n".format(x[i],y[i]))
    f.close()

    ## This is for 6A flex system which has some G-G contacts
   #nbins = 30
   #y,x = np.histogram(Y[X==0.0], bins=nbins); dx = x[1]-x[0]
   #x = x[:-1] + dx/2.;y = y/sum(y.astype(float))/dx
   #print("X position maxes: ", x[argrelextrema(y, np.greater)])
   #f = plt.figure(1, figsize = (1.0, 1.0))
   #ax, ct, leg = f.add_subplot(111), 0, []
   #ax.xaxis.set_major_locator(MaxNLocator())
   #ax.yaxis.set_major_locator(MaxNLocator())
   #plt.plot(x, y, color = "k")
   #ax.set_xlabel("$d_{gg}$",fontsize=7); #ax.set_xlim([0.2,0.8])
   #ax.set_ylabel("Probability",fontsize=7); #ax.set_ylim([0.,15])
   #plt.savefig(csv.csvfname[:-3]+'o_00.png',bbox_inches = 'tight',)

def main():                                                                        
    '''For a collection of data, get info from csv and then plot,
       usage: plot_vs_x.py csvStart nsep nlen niter sep1 sep2... len1 len2... 
                           iter1 iter2... ext datLoc'''
    csvname = sys.argv[1]; sep = int(sys.argv[2]); ln = int(sys.argv[3]);
    itr = sys.argv[4]
    
    plot_scatter(CSVFile(csvname), sep, ln, itr)

if __name__=="__main__":
    main()

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib, scipy
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

    Y = csv.dat[csv.key.index("dgg")]; 
    print(min(Y), max(Y))
   #y_okay = Y > sep-3.5; Y = Y[y_okay]
   #X = csv.dat[csv.key.index("dg1"),y_okay]/Y.astype(float)
    X = csv.dat[csv.key.index("dg1")]/Y.astype(float)
    plt.hist2d(X, Y, bins=(nbins,nbins), range=([0,1], [5,14]),cmap=plt.get_cmap('plasma'))
    ax.set_xlim([0.2,0.8]); ax.set_ylim([5,14])

    ax.set_xlabel("$x/d_{gg}$",fontsize= 8)
    ax.set_ylabel("$d_{gg} \, (\AA$)",fontsize= 8)
    plt.savefig(csv.csvfname[:-3]+'.png',bbox_inches = 'tight',)
    plt.close()

    # finding graphene separation dist(s)
    bil = [9,12] #list of initial seps that form 2 diff layers in flexible
    y,x = np.histogram(Y, bins=nbins); dx = x[1]-x[0]
    x = x[:-1] + dx/2.;y = y/sum(y.astype(float))/dx
    if sep in bil and itr != 'r':
        params = [sep-2, sep-1, 1.0, 1.0, 1.0, 1.0]
        fp,_ = scipy.optimize.curve_fit(double, x, y, p0=params)
        fit = double(x, *fp);
        tit = "$d_{{gg1}}$: {0:.4f} $\AA$, $d_{{gg2}}$: {1:.4f} $\AA$".format(fp[0],fp[1])
    else:
        params = [np.mean(Y), 1.0, -1.0]
        fp,_ = scipy.optimize.curve_fit(skew, x, y, p0=params)
        fit = skew(x, *fp);
        tit = "$d_{{gg}}$: {0:.4f} $\AA$".format(fp[0])
    print(tit)
    f = plt.figure(1, figsize = (1.0, 1.0))
    ax, ct, leg = f.add_subplot(111), 0, []
    ax.xaxis.set_major_locator(MaxNLocator())
    ax.yaxis.set_major_locator(MaxNLocator())
    matplotlib.rcParams['font.size'] = 5; nbins = 200
    ax.bar(x, y, width=dx,color = "y",edgecolor = "none"); ax.plot(x,fit)
    ax.set_xlabel("$d_{gg} \, (\AA$)",fontsize=7);ax.set_xlim([6.,16.])
    ax.set_ylabel("Probability ($1/\AA$)",fontsize=7);ax.set_ylim([0.,2.5])
    ax.set_title(tit)
    plt.savefig(csv.csvfname[:-3]+'g_sep_fit.png',bbox_inches = 'tight',)
    plt.close()

    ft_fl = open(csv.csvfname[:-3]+'g_sep_fit.dat', 'w')
    ft_fl.write("bin,dat,fit\n")
    for xx in range(len(x)): 
        ft_fl.write("{0:.3f},{1:.4f},{2:.4f}\n".format(x[xx],y[xx],fit[xx]))
    ft_fl.close()

    bil = [9,12] #list of initial seps that form 2 diff layers in flexible
    y,x = np.histogram(X, bins=nbins); dx = x[1]-x[0]
    x = x[:-1] + dx/2.;y = y/sum(y.astype(float))/dx
    if sep > 11 or (sep == 9 and itr != 'r'):
        params = [0.37, 0.5, 0.60, 0.01, 0.01, 0.01, 0.15, 0.45, 0.45]
        fp,_=scipy.optimize.curve_fit(triple,x,y,p0=params); fit=triple(x, *fp);
        print(fp)
        tit = "$x/d_{{gg1}}$: {0:.2f}$\AA$, $x/d_{{gg2}}$: {1:.2f}$\AA$, $x/d_{{gg3}}$: {2:.2f}$\AA$".format(fp[0],fp[1],fp[2])
    elif sep > 8:
        params = [0.30, 0.69, 0.01, 0.01, 0.4, 0.4]
        fp,_ = scipy.optimize.curve_fit(double, x, y, p0=params)
        fit = double(x, *fp);
        print(fp)
        tit = "$x/d_{{gg1}}$: {0:.2f}$\AA$, $x/d_{{gg2}}$: {1:.2f}$\AA$".format(fp[0],fp[1])
    else:
        params = [np.mean(X), 1.0, -1.0]
        fp,_ = scipy.optimize.curve_fit(skew, x, y, p0=params)
        fit = skew(x, *fp);
        tit = "$x/d_{{gg}}$: {0:.2f}$\AA$".format(fp[0])
    print(tit)
    f = plt.figure(1, figsize = (1.0, 1.0))
    ax, ct, leg = f.add_subplot(111), 0, []
    ax.xaxis.set_major_locator(MaxNLocator())
    ax.yaxis.set_major_locator(MaxNLocator())
    matplotlib.rcParams['font.size'] = 5; nbins = 200
    ax.bar(x, y, width=dx,color = "y",edgecolor = "none"); ax.plot(x,fit)
    ax.set_xlabel("$x/d_{gg}$",fontsize=7); ax.set_xlim([0.2,0.8])
    ax.set_ylabel("Probability",fontsize=7); ax.set_ylim([0.,15])
    ax.set_title(tit,fontsize=5)
    plt.savefig(csv.csvfname[:-3]+'o_sep_fit.png',bbox_inches = 'tight',)
    plt.close()

    ft_fl = open(csv.csvfname[:-3]+'o_sep_fit.dat', 'w')
    ft_fl.write("bin,dat,fit\n")
    for xx in range(len(x)): 
        ft_fl.write("{0:.3f},{1:.4f},{2:.4f}\n".format(x[xx],y[xx],fit[xx]))
    ft_fl.close()

def main():                                                                        
    '''For a collection of data, get info from csv and then plot,
       usage: plot_vs_x.py csvStart nsep nlen niter sep1 sep2... len1 len2... 
                           iter1 iter2... ext datLoc'''
    csvname = sys.argv[1]; sep = int(sys.argv[2]); ln = int(sys.argv[3]);
    itr = sys.argv[4]
    
    plot_scatter(CSVFile(csvname), sep, ln, itr)

if __name__=="__main__":
    main()

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib, scipy
from matplotlib.ticker import *
MaxNLocator.default_params['nbins']=5

from csvfile import CSVFile
from gauss_fits import skew, double, normal

colorL = [[0,0,0], [0,0,1], [1, 0,0], [.93, .53, .18]]

def plot_scatter(csv, sep, ln, itr):
    '''Using data from a histogram, plot several'''
    nbins = 150
    f = plt.figure(1, figsize = (1.5, 1.5)); ax = f.add_subplot(111)
    ax.xaxis.set_major_locator(MaxNLocator())
    ax.yaxis.set_major_locator(MaxNLocator())

    angles = csv.dat[csv.find_keyword('the1')].flatten(); rng = [0, 180]
    y,x = np.histogram(angles, bins=nbins, range=rng); dx = x[1]-x[0]
    x = x[:-1] + dx/2.;y = y/sum(y.astype(float))/dx
    ax.bar(x, y, width=dx, color = "y",edgecolor = "none");
    ax.set_xlim(rng); ax.set_ylim([0,0.03])

    ax.set_xlabel(r"3 body angle ($^{\circ}$)",fontsize= 8)
    ax.set_ylabel(r"Probability (1/$^{\circ}$)",fontsize= 8)
    plt.savefig(csv.csvfname[:-4]+'.png',bbox_inches = 'tight',)

   #params = [np.pi/2.,np.pi*(2./3.), 1.0, 1.0, 1.0, 1.0]
   #fp,_ = scipy.optimize.curve_fit(double, x, y, p0=params)
   #fit = double(x, *fp);  ax.plot(x,fit)
   #tit = "$A_{{3B_1}}$: {0:.2f} rad, $A_{{3B_2}}$: {1:.2f} rad".format(
   #       fp[0],fp[1])
    params = [np.mean(angles), 21.5,-0.01]
    plt.plot((90.,90.), (0,10), 'k-')
    plt.plot((109.,109.), (0,10), 'k-')
    fp,_ = scipy.optimize.curve_fit(skew, x, y, p0=params)
    print(np.mean(angles))
    fit = skew(x, *params);  ax.plot(x,fit)
    tit = r'$A_{{3B}}$: {0:.2f}$^{{\circ}}$'.format(fp[0])
    print(tit, fp);  ax.set_title(tit)
    plt.savefig(csv.csvfname[:-4]+'_fit.png',bbox_inches='tight'); plt.close()

    # Plotting and fitting NN distribution
    f = plt.figure(1, figsize = (1.5, 1.5)); ax = f.add_subplot(111)
    ax.xaxis.set_major_locator(MaxNLocator())
    ax.yaxis.set_major_locator(MaxNLocator())
    dists = csv.dat[csv.find_keyword('dis')].flatten()
    print(dists.shape)
    y,x = np.histogram(dists, bins=nbins, range=[1.5,4.]); dx = x[1]-x[0]
    x = x[:-1] + dx/2.;y = y/sum(y.astype(float))/dx
    ax.bar(x, y, width=dx, color = "y",edgecolor = "none");
    ax.set_xlim([1.5,4.0]); ax.set_ylim([0, 6])

    ax.set_xlabel("Near neigh dist ($\AA$)",fontsize= 8)
    ax.set_ylabel("Probability ($1/\AA$)",fontsize= 8)
    plt.savefig(csv.csvfname[:-10]+'dists.png',bbox_inches = 'tight',)

    params = [np.mean(dists), 1.0, -1.0]
    fp,_ = scipy.optimize.curve_fit(skew, x, y, p0=params)
    fit = skew(x, *fp); tit = "$d_{{NN}}$: {0:.2f}$\AA$".format(fp[0])
    print(tit, fp); ax.plot(x,fit);  ax.set_title(tit)
    plt.savefig(csv.csvfname[:-10]+'dists_fit.png',bbox_inches = 'tight',)

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

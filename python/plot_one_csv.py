import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib, scipy

from csvfile import CSVFile
from gauss_fits import skew, double

colorL = [[0,0,0], [0,0,1], [1, 0,0], [.93, .53, .18]]

def plot_scatter(csv, sep, ln, itr):
    '''Using data from a histogram, plot several'''
    matplotlib.rcParams['font.size'] = 5; nbins = 200
    f = plt.figure(1, figsize = (1.5, 1.5))
    ax, ct, leg = f.add_subplot(111), 0, []
    leg = str(sep)+r'$\AA$ Sep, '+'L='+str(ln)+r'$\AA$'

    Y = csv.dat[csv.key.index("dgg")]; 
    y_okay = Y < sep+4.0; Y = Y[y_okay]
    X = csv.dat[csv.key.index("dg1"),y_okay]/Y.astype(float)
    plt.hist2d(X, Y, bins=(nbins,nbins), range=([0,1], [5,14]),cmap=plt.get_cmap('plasma'))
    ax.set_xlim([0.2,0.8]); ax.set_ylim([5,14])

    ax.set_title(leg, fontsize=10)
    ax.set_xlabel("$x/d_{gg}$",fontsize=10)
    ax.set_ylabel("$d_{gg} \, (\AA$)",fontsize=10)
    plt.savefig(csv.csvfname[:-3]+'.png', format='png',
                    bbox_inches = 'tight', dpi=300) 
    plt.close()

    # finding graphene separation dist(s)
    y,x = np.histogram(Y, bins=nbins); dx = x[1]-x[0]
    x = x[:-1] + dx/2.;y = y/sum(y.astype(float))/dx
    if sep > 8:
        params = [sep-1.,sep-0.5, 1.0, 1.0, 1.0, 1.0]
        fp,_ = scipy.optimize.curve_fit(double, x, y, p0=params)
        fit = double(x, *fp)
        tit = "$d_{{gg1}}$: {0:.1f}, $d_{{gg2}}$: {1:.1f}".format(fp[0],fp[1])
        print(tit)
    else:
        params = [sep, 1.0, -1.0]
        fp,_ = scipy.optimize.curve_fit(skew, x, y, p0=params)
        fit = skew(x, *fp)
        print(fp)
        tit = "$d_{{gg}}$: {0:.4f}".format(fp[0])
    f = plt.figure(1, figsize = (1.5, 1.5))
    ax, ct, leg = f.add_subplot(111), 0, []
    ax.bar(x, y, width=dx, color = "y");  ax.plot(x,fit)
    ax.set_xlabel("$d_{gg} \, (\AA$)",fontsize=10)
    ax.set_ylabel("Frequency",fontsize=10)
    ax.set_title(tit, fontsize=10)
    plt.savefig(csv.csvfname[:-3]+'g_sep_fit.png', format='png',
                    bbox_inches = 'tight', dpi=300) 

def main():                                                                        
    '''For a collection of data, get info from csv and then plot,
       usage: plot_vs_x.py csvStart nsep nlen niter sep1 sep2... len1 len2... 
                           iter1 iter2... ext datLoc'''
    csvname = sys.argv[1]; sep = int(sys.argv[2]); ln = int(sys.argv[3]);
    itr = int(sys.argv[4]);
    
    plot_scatter(CSVFile(csvname), sep, ln, itr)

if __name__=="__main__":
    main()

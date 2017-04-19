import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

from csvfile import CSVFile

colorL = [[0,0,0], [0,0,1], [1, 0,0], [.93, .53, .18]]

def plot_scatter(csv, sep, ln, itr):
    '''Using data from a histogram, plot several'''
    matplotlib.rcParams['font.size'] = 5; nbins = 70
    f = plt.figure(1, figsize = (2.0, 2.0))
    ax, ct, leg = f.add_subplot(111), 0, []

    # distances = X, angle = Y
    Y = csv.dat[5::3].flatten()
    X = csv.dat[3:-1:3].flatten()/Y
    xbins = np.linspace(0, 1, nbins)
    ybins = np.linspace(0, 6, nbins)
    heatmap, xedges, yedges = np.histogram2d(X, Y, bins=(xbins,ybins))
    print(heatmap.shape)
    heatmap = heatmap.astype(float)
    extent = [0, 1.0, 0., 6.]
    ax.imshow(heatmap.T,aspect='auto', extent=extent)
    leg = str(sep)+r'$\AA$ Sep, '+'L='+str(ln)+r'$\AA$'
    ax.set_xlim([0,1])
    titl = leg
    ax.set_title(titl)
   #ax.legend(leg, loc = 9, ncol = 1,
   #    columnspacing = 0.4,
   #    fontsize =  7 ,
   #    handletextpad = 0.2,
   #    handlelength = 1.3,
   #    borderaxespad = -0.9,
   #   #bbox_to_anchor = bbox
   #    )
   #plt.show()
    plt.savefig(csv.csvfname[:-3]+'.png', format='png',
                    bbox_inches = 'tight', dpi=200) 
    plt.close()

def main():                                                                        
    '''For a collection of data, get info from csv and then plot,
       usage: plot_vs_x.py csvStart nsep nlen niter sep1 sep2... len1 len2... 
                           iter1 iter2... ext datLoc'''
    csvname = sys.argv[1]; sep = int(sys.argv[2]); ln = int(sys.argv[3]);
    itr = int(sys.argv[4]);
    
    plot_scatter(CSVFile(csvname), sep, ln, itr)

if __name__=="__main__":
    main()

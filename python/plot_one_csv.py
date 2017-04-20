import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

from csvfile import CSVFile

colorL = [[0,0,0], [0,0,1], [1, 0,0], [.93, .53, .18]]

def plot_scatter(csv, sep, ln, itr):
    '''Using data from a histogram, plot several'''
    matplotlib.rcParams['font.size'] =  8; nbins = 200
    f = plt.figure(1, figsize = (1.5, 1.5))
    ax, ct, leg = f.add_subplot(111), 0, []

    Y = csv.dat[csv.key.index("dgg")]; 
    X = csv.dat[csv.key.index("dg1")]/Y.astype(float)
   #print(max(Y), max(X), min(Y), min(X), np.mean(X), np.mean(Y))
   #heatmap, xedges, yedges = np.histogram2d(X, Y, bins=(nbins,nbins), range=([0,1], [4,14]))
   #heatmap = heatmap.astype(float)
   #extent = [xedges[0], xedges[-1],yedges[0], yedges[-1]]
   #print(extent, heatmap.shape)#heatmap)
   #ax.imshow(heatmap, aspect='auto', extent=extent, cmap=plt.get_cmap('plasma'))
    plt.hist2d(X, Y, bins=(nbins,nbins), range=([0,1], [5,14]),cmap=plt.get_cmap('plasma'))
   #plt.scatter(X,Y)
   #plt.hist(Y,bins=nbins)
    leg = str(sep)+r'$\AA$ Sep, '+'L='+str(ln)+r'$\AA$'
    ax.set_xlim([0.2,0.8])
    ax.set_ylim([5,14])
    ax.set_title(leg)
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

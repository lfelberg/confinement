import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import *
MaxNLocator.default_params['nbins']=5

from csvfile import CSVFile
from dics import colorL, dens
A2_TO_M2 = 1e-20
PS_TO_S  = 1e12
A2_PS_TO_M2_S = A2_TO_M2*PS_TO_S 

def plot_scatter(csv, sep, nsol):
    '''Using data from a histogram, plot several'''
    f = plt.figure(1, figsize = (1.5, 1.5))
    ax, ct, leg = f.add_subplot(111), 0, []; mav = 0
    ax.xaxis.set_major_locator(MaxNLocator()); ax.yaxis.set_major_locator(MaxNLocator())

    st = 400
    en = 800
    xd = csv.dat[0,st:en]; yd = csv.dat[1,st:en]
    print(csv.key)

    # fitting and plotting 2D water MSD
    m,b = np.polyfit(xd,yd,1)
    print("Water D: {0}".format(m/4.*A2_PS_TO_M2_S*(1e9)))
    ax.plot(csv.dat[0], csv.dat[1])
   #ax.set_xlim([0,800])
   #ax.set_ylim([0,800])
    ax.set_xlabel("Time (ps)",fontsize=10); 
    ax.set_ylabel("$MSD_{||} \,\, (\AA^2)$",fontsize=10)
    plt.savefig(csv.csvfname[:-4]+'_water.png', bbox_inches = 'tight',) 
    plt.close()

    # Fitting and plotting 2D solute MSD
    f = plt.figure(1, figsize = (1.5, 1.5))
    ax, leg = f.add_subplot(111), []; nsam = float(len(csv.dat)-3)
    ax.xaxis.set_major_locator(MaxNLocator()); ax.yaxis.set_major_locator(MaxNLocator())
    for i in range(3,len(csv.dat)):
        xd = csv.dat[0,st:en]; yd = csv.dat[i,st:en]
        m,b = np.polyfit(xd,yd,1)
        print("{0} D: {1}".format(csv.key[i],m/4.*A2_PS_TO_M2_S*(1e9)))
        mav += m/4.*A2_PS_TO_M2_S*(1e9)
        ax.plot(csv.dat[0], csv.dat[i])
    print("{0},{1:.7f}".format(dens[float(sep)][0],mav/nsam))

    ax.set_xlim([0,800])
    ax.set_ylim([0,800])
    ax.set_xlabel("Time (ps)",fontsize=10); 
    ax.set_ylabel("$MSD_{||} \,\, (\AA^2)$",fontsize=10)
    plt.savefig(csv.csvfname[:-4]+'_solute.png', bbox_inches = 'tight',) 
    plt.close()

def main():
    '''For a collection of data, get info from csv and then plot,
       usage: plot_msd_many.py csvfn sep nsol'''
    csvname = sys.argv[1]; sep = int(sys.argv[2]); ln = int(sys.argv[3]);
    
    cs = CSVFile(csvname)
   #print("Ext {0}, loc {1} = {2}".format(ext, loc, cs[-1].key), sep, ln, itr)
    plot_scatter(cs, sep, ln)

if __name__=="__main__":
    main()

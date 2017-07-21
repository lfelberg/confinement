import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import *
MaxNLocator.default_params['nbins']=4

from csvfile import CSVFile
from dics import dens, nsol, colorL
from gauss_fits import skew, double, normal

def plot_scatter(csv):
    '''Using data from a histogram, plot several'''
    f = plt.figure(1, figsize = (1.5, 1.5)); ax = f.add_subplot(111)
   #f = plt.figure(1, figsize = (3.3, 3.3)); ax = f.add_subplot(111)
    plt.rcParams.update({'font.size': 8}); ncl = 3; ymin = 0
    ax.xaxis.set_major_locator(MaxNLocator())
    ax.yaxis.set_major_locator(MaxNLocator())

    qs = csv.dat[csv.find_keyword('q')]; qm = np.mean(qs, axis = 0)
    seps = csv.dat[csv.find_keyword("sep")][0]; 
    lsep = list(set(list(seps.flatten()))); lsep.sort()
    nsol=csv.dat[csv.find_keyword("nsol")][0]

    # compressibility
    if "compressibility" in csv.csvfname:
        ymx = 65; bbox = (1.01, 1.30); texty=ymx*0.89; plt_l = 1;bulk_v = 48.1
        ax.set_ylabel(r"$\kappa_T \times 10^{-11}\,\,\, (Pa^{-1})$",fontsize= 8)

    # graphene graphene separation
    elif "d_gg" in csv.csvfname:
        bbox, ncl = (0.59,0.96), 1; plt_l = 0
        qm = csv.dat[csv.find_keyword('q')][0]
        ax.set_xlim([0.,16]);ax.set_ylim([6,15])
        ax.set_ylabel(r"$\langle d_{gg} \rangle \,\, (\AA)$",fontsize= 10)

    # diffusion coefficients
    else:
        ymx = 5; bbox = (1.00,0.65); texty = 6.4; plt_l = 1; bulk_v = 2.4725
        ax.set_ylabel(r"$\mathcal{D}_{||} \times 10^9\, (m^2/s)$",fontsize= 10)
        ax.set_yscale("log", nonposy='clip'); ncl = 1; ymin = 1e-6

   ## PLotting the delineation between layers 1-2 (rho=0.1465), 2-3 (0.265),3-4
    if plt_l == 1:
        x = np.linspace(0,25,70); y = np.ones(70)*bulk_v
        plt.plot(x,y, 'r--',dashes=(1.5,0.9),linewidth=0.7,label="bulk,\n298K")
        ax.set_xlim([0,16]); ax.set_ylim([ymin,ymx])

    for i in range(len(lsep)):
        rho = dens[lsep[i]][0]; col = colorL[i*2];
        plt.plot(nsol[seps == lsep[i]], qm[seps == lsep[i]], color = col,
                 label="{0:.3f}".format(rho))
                #label=r"$\rho_{{2D}}^W={0:.3f}$".format(rho))
    ax.set_xlim([0,16])
   #ax.legend(ncol = ncl, fontsize=5, columnspacing = 0.2,handletextpad= 0.25,
   #          bbox_to_anchor = bbox, handlelength = 1.1,numpoints = 1,
   #          borderaxespad= 0.2)
    ax.set_xlabel(r"$\rho_{2D}^B$",fontsize= 10)
    plt.savefig(csv.csvfname[:-4]+'.png',bbox_inches = 'tight')
    plt.close()

def main():                                                                        
    '''For a collection of data, get info from csv and then plot,
       '''
    csvname = sys.argv[1]
    plot_scatter(CSVFile(csvname))

if __name__=="__main__":
    main()

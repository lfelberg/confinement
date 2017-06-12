import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import *
MaxNLocator.default_params['nbins']=2

from csvfile import CSVFile
from gauss_fits import skew, double, normal

def plot_scatter(csv):
    '''Using data from a histogram, plot several'''
    f = plt.figure(1, figsize = (1.5, 1.5)); ax = f.add_subplot(111)
    plt.rcParams.update({'font.size': 8}); ncl = 3
    ax.xaxis.set_major_locator(MaxNLocator())
    ax.yaxis.set_major_locator(MaxNLocator())

    dens  = csv.dat[csv.find_keyword('dens')].flatten(); 
    rigid = csv.dat[-1].flatten(); flexb = csv.dat[1].flatten(); 
    print("rigid col: {0}\n flex col: {1}".format(csv.key[-1],csv.key[1]))

    # compressibility
    if "compressibility" in csv.csvfname:
        xm = 65; bbox = (1.01, 1.25); texty=xm*0.89; plt_l = 1;bulk_v = 48.1
        ax.set_ylabel(r"$\kappa_T \times 10^{-11}\,\,\, (Pa^{-1})$",fontsize= 8)

    # graphene graphene separation
    elif "d_gg" in csv.csvfname:
        plt_l = 0; ff, df, flx, fl = [], [], [], []; bbox, ncl = (0.52,0.99), 1
        flexb = csv.dat[csv.find_keyword('dgg_f')]
        for j in range(flexb.shape[1]): 
            if flexb[1][j] != 0.:
                ff.append(dens[j]); fl.append(flexb[0][j])
                ff.append(dens[j]); fl.append(flexb[1][j])
            else:  df.append(dens[j]); flx.append(flexb[0][j])
        plt.plot(ff, fl, 'mD', label = "multi")
        plt.plot(df, flx, 'k.', label = "flexible")
        plt.plot(dens, rigid, 'b.',label = "rigid")
        ax.set_xlim([0.,0.40]);ax.set_ylim([3,15])
        ax.set_ylabel(r"$\langle d_{gg} \rangle \,\, (\AA)$",fontsize= 10)

    # diffusion coefficients
    else:
        xm = 30; bbox = (1.02,1.25); texty = 6.4; plt_l = 1; bulk_v = 2.4725
        ax.set_ylabel(r"$\mathcal{D}_{||} \times 10^9\, (m^2/s)$",fontsize= 10)
        ax.set_yscale("log", nonposy='clip')

    # PLotting the delineation between layers 1-2 (rho=0.1465), 2-3 (0.265),3-4
    if plt_l == 1:
        plt.plot(dens,flexb,'k.--',dashes = (2,1),label = "flexible")
        plt.plot(dens,rigid,'b.--',dashes = (2,1),label = "rigid")
        x = np.linspace(0,1,70); y = np.ones(70)*bulk_v
        plt.plot(x,y, 'r--',dashes=(1.5,0.9),linewidth=0.7,label="bulk,\n298K")
        y=np.linspace(0,xm,70); ax.set_xlim([0.,0.55]); ax.set_ylim([0,xm])
        lay = [0.1465, 0.245, 0.3525]; txl = [0.06,0.17,0.27, 0.405]
        clr = [0.439, 0.502, 0.565]
        for ll in range(len(lay)):
            plt.plot(np.ones(70)*lay[ll],y,'--',color=clr,
                     dashes=(2.5,0.9),linewidth=0.7)
            plt.text(txl[ll],texty,str(ll)+'L',fontsize=10, color = clr)
        plt.text(txl[ll],texty,str(ll)+'L',fontsize=10, color = clr)

    ax.legend(ncol = ncl, fontsize=6, columnspacing = 0.2,handletextpad= 0.25,
              bbox_to_anchor = bbox, handlelength = 2.1,
              borderaxespad= 0.2)
    ax.set_xlabel(r"$\rho_{2D}$",fontsize= 10)
    plt.savefig(csv.csvfname[:-4]+'.png',bbox_inches = 'tight')
    plt.close()

def main():                                                                        
    '''For a collection of data, get info from csv and then plot,
       '''
    csvname = sys.argv[1]
    plot_scatter(CSVFile(csvname))

if __name__=="__main__":
    main()

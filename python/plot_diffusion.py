import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import *
MaxNLocator.default_params['nbins']=2

from csvfile import CSVFile
from gauss_fits import skew, double, normal

colorL = [[0,0,0], [0,0,1], [1,0,0], [.93, .53, .18]]

def plot_scatter(csv):
    '''Using data from a histogram, plot several'''
    nbins = 150
    f = plt.figure(1, figsize = (1.5, 1.5)); ax = f.add_subplot(111)
    plt.rcParams.update({'font.size': 8})
    ax.xaxis.set_major_locator(MaxNLocator())
    ax.yaxis.set_major_locator(MaxNLocator())

    dens  = csv.dat[csv.find_keyword('dens')].flatten(); 
    rigid = csv.dat[2].flatten();
    flexb = csv.dat[1].flatten(); 
    print("rigid col: {0}\n flex col: {1}".format(csv.key[2],csv.key[1]))
    plt.plot(dens,flexb,'k.--',dashes = (2,1),label = "flexible")
    plt.plot(dens,rigid,'b.--',dashes = (2,1),label = "rigid")

    # compressibility
    if "compressibility" in csv.csvfname:
        xm = 900; bbox = (1.01, 1.15); texty=xm*0.89
        ax.set_ylabel(r"$(\langle l_x^2 \rangle - \langle l_x \rangle^2)^{-1}\, (m^{-2})$",fontsize= 8)

    # diffusion coefficients
    else:
        xm = 30; bbox = (1.02,1.30); texty = 6.4
        x = np.linspace(0,1,70); y = np.ones(70)*2.4725
        plt.plot(x,y, 'r--',dashes=(1.5,0.9),linewidth=0.7,label="bulk,\n298K")
        ax.set_ylabel(r"$\mathcal{D}_{||} \times 10^9\, (m^2/s)$",fontsize= 10)
        ax.set_yscale("log", nonposy='clip')

    # PLotting the delineation between layers 1-2 (rho=0.1465), 2-3 (0.265),3-4
    y=np.linspace(0,xm,70);
    x1=np.ones(70)*0.1465; x2=np.ones(70)*0.26;x3=np.ones(70)*0.3525
    plt.plot(x1,y,'--',color=[0.439, 0.502, 0.565],dashes=(2.5,0.9),linewidth=0.7)
    plt.plot(x2,y,'--',color=[0.439, 0.502, 0.565],dashes=(2.5,0.9),linewidth=0.7)
    plt.plot(x3,y,'--',color=[0.439, 0.502, 0.565],dashes=(2.5,0.9),linewidth=0.7)
    plt.text(0.06, texty,'1L',fontsize=10, color = [0.439, 0.502, 0.565])
    plt.text(0.170,texty, '2L', fontsize=10, color = [0.439, 0.502, 0.565])
    plt.text(0.270,texty, '3L', fontsize=10, color = [0.439, 0.502, 0.565])
    plt.text(0.405,texty, '4L', fontsize=10, color = [0.439, 0.502, 0.565])

    ax.set_xlim([0.,0.55]); ax.set_ylim([0,xm])
    ax.legend(ncol = 3, fontsize=6, columnspacing = 0.2,handletextpad= 0.15,
              bbox_to_anchor = bbox, 
              borderaxespad= 0.2)
    ax.set_xlabel(r"$\rho_{2D}$",fontsize= 10)
   #ax.set_ylabel(r"$\rho_{3D} \, (g/cm^3)$",fontsize= 10)
    plt.savefig(csv.csvfname[:-4]+'.png',bbox_inches = 'tight')#,transparent=True)
    plt.close()

def main():                                                                        
    '''For a collection of data, get info from csv and then plot,
       usage: plot_vs_x.py csvStart nsep nlen niter sep1 sep2... len1 len2... 
                           iter1 iter2... ext datLoc'''
    csvname = sys.argv[1]
    plot_scatter(CSVFile(csvname))

if __name__=="__main__":
    main()

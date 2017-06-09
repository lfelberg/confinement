import sys
import numpy as np
from scipy.integrate import simps
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import *
MaxNLocator.default_params['nbins']=2

from csvfile import CSVFile
from util import get_gr_lst

def fourier(grcl):
    '''Fourier transform the g(r)'''
    dt_idx,_ = get_gr_lst(grcl.csvfname); n = len(grcl.dat[0])
    dx = grcl.dat[0,1]-grcl.dat[0,0]; r = grcl.dat[0][np.newaxis]
    gr = np.mean(grcl.dat[dt_idx],axis = 0)

    q = np.arange(0.15,5,0.01)[:,np.newaxis]; 
    pref=(np.sin(q*r)/q)*(gr[np.newaxis]-1.)
    sk = 1 + 2*np.pi*simps(pref,r,axis = -1)

    ft = (1./float(n))*np.fft.fft(1-gr)
    nu = np.fft.fftshift(np.fft.fftfreq(n, dx))
    Fk = np.fft.fftshift(ft);
    f,ax = plt.subplots(3,1, figsize = (1.5, 1.5), sharex=True)
   #ax, ct, leg = f.add_subpsot(211), 0, []
   #ax2, ct, leg = f.add_subplot(212), 0, []
    
    matplotlib.rcParams.update({'font.size': 8}) 
    ax[0].plot(nu,1+Fk.real)
    ax[1].plot(nu,1+Fk.imag)
    ax[2].plot(q,sk);#ax[2].set_ylim([-5,20])
    for i in range(len(ax)): 
        ax[i].set_xlim([0,4])
       #ax[i].set_ylim([0.0,1.1])
        ax[i].xaxis.set_major_locator(MaxNLocator())
        ax[i].yaxis.set_major_locator(MaxNLocator())
    plt.savefig("ft.png", format='png', bbox_inches = 'tight', dpi=300)
   #print(ft)
    plt.close()

def main():
    ''' Given a list of pairs for g(r) calcs, do this as a function as 
        distance from wall, but only 2D and only for your wall side '''
    grname=sys.argv[1];#sep=sys.argv[2]; ln=sys.argv[3]; itr=sys.argv[4]

   #nm = str(sep)+"_"+str(ln)+"_"+str(itr)
    grcl = CSVFile(grname)
    fourier(grcl)


if __name__=="__main__":
    main()

import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import *
MaxNLocator.default_params['nbins']=2

from csvfile import CSVFile
from dics import colorL, dens

def plot_scatter(plt_nm, csvL, sep, ln):
    '''Using data from a histogram, plot several'''
    f = plt.figure(1, figsize = (1.5, 1.5))
    ax, ct, leg, mn = f.add_subplot(111), 0, [], ""
   #matplotlib.rcParams.update({'font.size': 8})
    ax.xaxis.set_major_locator(MaxNLocator())
    ax.yaxis.set_major_locator(MaxNLocator())
    for i in range(len(csvL)):
        for j in range(len(csvL[i])):
           fn = csvL[i][j].csvfname.split("_")
           sep = int(fn[0][3:])
           if "bulk" != dens[sep][0]: lg = "{0:.3f}".format(dens[sep][0])
           else:                          lg = "bulk"
           ax.plot(csvL[i][j].dat[0],csvL[i][j].dat[1],
                   color = dens[sep][1], label=lg)
           mn += fn[0][3:] + "_";  ct += 1
   #ax.legend(ncol = 1, columnspacing = -0.4,
   #    fontsize =  6 , handletextpad = -0.1,
   #    handlelength = 1.2, borderaxespad = -0.9,
   #    bbox_to_anchor = (0.32,0.91),
   #    )
    ax.set_xlim([0,180]); ax.set_ylim([0,0.02])
    ax.set_xlabel(r"3 body angle ($^{\circ}$)",fontsize= 8)
    ax.set_ylabel(r"Probability (1/$^{\circ}$)",fontsize= 8)
    ax.yaxis.labelpad = -0.6; ax.xaxis.labelpad = -0.6
    fname = "ang_3B_"+mn[:-1]+".png"
    plt.savefig(fname, format='png', bbox_inches = 'tight', dpi=300) 
    plt.close()

def main():                                                                        
    '''For a collection of data, get info from csv and then plot,
       usage: plot_vs_x.py csvStart nsep nlen iter sep1 sep2... len1 len2... 
                           ext datLoc'''
    csvname = sys.argv[1]; nsep = int(sys.argv[2]); nlen = int(sys.argv[3]);
    itr=int(sys.argv[4]); spS=5; spE=spS+nsep; lnE=spE+nlen; sep,ln = [], []
    
    for i in range(spS, spE): sep.append(int(sys.argv[i]))
    for i in range(spE, lnE): ln.append(int(sys.argv[i]))
    ext = sys.argv[lnE]
    csvL = []
    for i in range(nsep):
        csL = []
        for j in range(nlen):
            csL.append(CSVFile(csvname+str(sep[i])+"_"+str(ln[j])+"_"+str(itr)+ext))
        csvL.append(csL)
    plot_scatter(csvname, csvL, sep, ln)

if __name__=="__main__":
    main()

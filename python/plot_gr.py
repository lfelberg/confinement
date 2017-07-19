import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import *
MaxNLocator.default_params['nbins']=4

from csvfile import CSVFile
from plot_util import get_gr_lst
from dics import colorL, dens

def plot_scatter(plt_nm, csvL, sep, ln):
    '''Using data from a histogram, plot several'''
    f = plt.figure(1, figsize = (1.5, 1.5))
    ax, ct, leg, mn = f.add_subplot(111), 0, [], ""
    matplotlib.rcParams.update({'font.size': 8})
    ax.xaxis.set_major_locator(MaxNLocator())
    ax.yaxis.set_major_locator(MaxNLocator())
    for i in range(len(csvL)):
        for j in range(len(csvL[i])):
           nfct = 1.0
           cls,nfct=get_gr_lst(csvL[i][j].csvfname); 
           dat = csvL[i][j].dat
           fn = csvL[i][j].csvfname.split("_")
           mn += fn[3] + "_"; dn = float(fn[3])
           ct=0;#print(cls)

          #cls = []; print(len(dat))
          #for k in range(1,len(dat)):
          #    if max(dat[k]) <25.5 and sum(dat[k]) > 0.1: cls.append(k)
          #ct=0; print(cls)
          #for k in cls:
          #    print(csvL[i][j].key[k])
          #    ax.plot(dat[0], dat[k], label=str(k), color = colorL[ct+3])
          #    ct += 1

           if type(dens[dn][0]) == float: lg = "{0:.3f}".format(dens[dn][0])
           else:               lg = dens[dn][0]
           
           if dens[dn][2] == "--": dsh = (2,1)
           else:                   dsh = (None,None)
           ax.plot(dat[0], np.mean(dat[cls],axis=0)/nfct, dens[dn][2],
                   color = dens[dn][1], dashes = dsh, label=lg)
           ct += 1
    if dn > 13 and len(sep) > 1: 
        ax.legend(ncol = 5, columnspacing = 0.4,
            fontsize =  5 , handletextpad = 0.2,
            handlelength = 1.3, borderaxespad = -0.9,
            bbox_to_anchor = (0.90,0.90),
            )
    ax.set_xlim([0,12.0]); ax.set_ylim([0,4.0])
    if "_1_2.csv" in csvL[i][j].csvfname:
        ax.set_ylabel(r"$g_{{OH}}(r)$",fontsize=12)
    else: ax.set_ylabel(r"$g_{{OO}}(r)$",fontsize=12)
    ax.set_xlabel("$r \, (\AA)$",fontsize=12)
    ax.yaxis.labelpad = -0.6; ax.xaxis.labelpad = -0.6
    fname = plt_nm+mn+csvL[i][j].csvfname[-7:-4]+".png"
    plt.savefig(fname, format='png', bbox_inches = 'tight', dpi=300) 
    plt.close()

def main():                                                                        
    '''For a collection of data, get info from csv and then plot,
       usage: plot_gr.py csvStart nsep nlen iter sep1 sep2... len1 len2... 
                           ext datLoc'''
    csvname = sys.argv[1]; nsep = int(sys.argv[2]); nlen = int(sys.argv[3]);
    itr=sys.argv[4]; spS=5; spE=spS+nsep; lnE=spE+nlen; sep,ln = [], []
    
    for i in range(spS, spE): sep.append(sys.argv[i])
    for i in range(spE, lnE): ln.append(int(sys.argv[i]))
    ext = sys.argv[lnE]
    csvL = []
    for i in range(nsep):
        csL = []
        for j in range(nlen):
            csL.append(CSVFile(csvname+sep[i]+"_"+str(ln[j])+"_"+itr+ext))
        csvL.append(csL)
    plot_scatter(csvname, csvL, sep, ln)

if __name__=="__main__":
    main()

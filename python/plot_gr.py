import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import *
MaxNLocator.default_params['nbins']=4

from csvfile import CSVFile
from dics import colorL, dens

def get_plt_lst(csv_nm):
    ''' Given a filename for g(r) csv, find sep size and use that to get
        columns to plot'''
    sep_siz = int(csv_nm.split("_")[3])
    dim = csv_nm.split("_")[2]
    if dim == "3D": cls = [1]
    # flexible
    elif sep_siz == 9:  cls = [2,4,15,28,30,32]
    elif sep_siz == 10: cls = [4, 6, 7, 17, 18,]
    elif sep_siz == 12: cls = [9,21,22,41,]
    # rigid
   #elif sep_siz == 9:  cls = [5,15,18,27]
   #elif sep_siz == 10: cls = [3,4,5,14,]
   #elif sep_siz == 12: cls = [3,18,19,31,32,7,8,17,20,]

    elif sep_siz == 11: cls = [4, 7, 18, 36]
    elif sep_siz == 13: cls = [5, 8, 12, 13, 22, 25, 34, 45]
    elif sep_siz == 14: cls = [2,19,29,38,4,7,8,11,12,23,24,] #27]
    elif sep_siz == 16: cls = [3,6,7,11,14,19,31,32,40,46]

    else: cls = [1] # dont have this sep size saved
    return cls

def plot_scatter(plt_nm, csvL, sep, ln):
    '''Using data from a histogram, plot several'''
    f = plt.figure(1, figsize = (1.5, 1.5))
    ax, ct, leg, mn = f.add_subplot(111), 0, [], ""
    matplotlib.rcParams.update({'font.size': 8})
    ax.xaxis.set_major_locator(MaxNLocator())
    ax.yaxis.set_major_locator(MaxNLocator())
    cos = [colorL[0], colorL[1], colorL[-3]]
    for i in range(len(csvL)):
        for j in range(len(csvL[i])):
           cls = get_plt_lst(csvL[i][j].csvfname); dat = csvL[i][j].dat
           fn = csvL[i][j].csvfname.split("_")
           mn += fn[3] + "_"
          #cls = []
          #for k in range(len(dat)):
          #    if max(dat[k]) <10.5 and sum(dat[k]) > 0.1: cls.append(k)
          #ct=0; print(cls)
          #for k in cls:
          #    ax.plot(dat[0], dat[k]+0 , label="x="+csvL[i][j].key[k][1:-2], color = colorL[ct])
          #    ct += 1
           ax.plot(dat[0], np.mean(dat[cls],axis=0), 
                   color = cos[ct], label="{0:.3f}".format(dens[int(fn[3])]))
          #ax.plot(dat[0], np.mean(dat[cls],axis=0)+ct*4,  color = cos[ct])
          #plt.text(6.0, ct*4+2.05, r'$\rho_{{2D}}={0:.3f}$'.format(dens[int(fn[3])]),
          #         color = cos[ct])
           ct += 1
   #ax.legend(ncol = 1, columnspacing = 0.4,
   #    fontsize =  5 , handletextpad = 0.2,
   #    handlelength = 1.3, borderaxespad = -0.9,
   #    bbox_to_anchor = (0.9,0.9),
   #    )
    ax.set_xlim([0,12.0]); ax.set_ylim([0,4.0])
   #ax.set_xlim([0,12.0]); ax.set_ylim([0,13.0])
    if "_1_2.csv" in csvL[i][j].csvfname:
        ax.set_ylabel(r"$g_{{O-H}}(R)$ (a.u.)",fontsize=12)
    else: ax.set_ylabel(r"$g_{{\mathrm{{O-O}}}}$(R) (a.u.)",fontsize=12)
    ax.set_xlabel("R ($\AA$)",fontsize=12)
    ax.yaxis.labelpad = -0.6; ax.xaxis.labelpad = -0.6
    fname = plt_nm+mn+csvL[i][j].csvfname[-7:-4]+".png"
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

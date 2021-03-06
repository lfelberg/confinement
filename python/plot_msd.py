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

def plot_scatter(plt_nm, csvL, loc, sep, ln, itr):
    '''Using data from a histogram, plot several'''
    f = plt.figure(1, figsize = (1.5, 1.5))
    ax, ct, leg, stn = f.add_subplot(111), 0, [], ""
    ax.xaxis.set_major_locator(MaxNLocator()); ax.yaxis.set_major_locator(MaxNLocator())

    for i in range(len(csvL)):
        for j in range(len(csvL[i])):
            mav = 0.0 
            for k in range(len(csvL[i][j])):
                st = 800
                en = 1600
                xd = csvL[i][j][k].dat[0,st:en]; yd = csvL[i][j][k].dat[loc,st:en]
                nsam = float(len(csvL[i][j]))

                if sep[i] == 37:
                     xms = csvL[i][j][k].dat[3,st:en]
                     yms = csvL[i][j][k].dat[5,st:en]
                     zms = csvL[i][j][k].dat[9,st:en]
                     m,b = np.polyfit(xd,yms+xms,1)
                     print(m/4.*A2_PS_TO_M2_S*(1e9))
                     m,b = np.polyfit(xd,xms+zms,1)
                     mav += m/4.*A2_PS_TO_M2_S*(1e9)
                     print(m/4.*A2_PS_TO_M2_S*(1e9))
                     m,b = np.polyfit(xd,yms+zms,1)
                     mav += m/4.*A2_PS_TO_M2_S*(1e9)
                     print(m/4.*A2_PS_TO_M2_S*(1e9))
                     mav += m/4.*A2_PS_TO_M2_S*(1e9)
                     nsam = len(csvL[i][j])*3.

                m,b = np.polyfit(xd,yd,1)
                stn += ",{0:.7f}".format(m/4.*A2_PS_TO_M2_S*(1e9))
                mav += m/4.*A2_PS_TO_M2_S*(1e9)
                ax.plot(csvL[i][j][k].dat[0], 
                        csvL[i][j][k].dat[loc], color=colorL[ct])
               #ax.plot(xd, xd*m + b, "r")
               #leg += [r'$\rho_{{2D}}=${0:.2f}'.format(dens[sep[i]])]
                ct += 1
            p_all = "{0},{1},{2}{3}".format(sep[i],ln[j],itr[0].split("_")[0],stn)
            print(p_all)  #mav/nsam
   #ax.set_xlim([0,800])
   #ax.set_ylim([0,800])
    bbox = [1.1, 0.95]
   #ax.legend(leg, loc = 2, ncol = 1, columnspacing = 0.4,
   #    fontsize =  7 , handletextpad = 0.2, handlelength = 1.3,
   #    borderaxespad =  0.2, )
    ax.set_xlabel("Time (ps)",fontsize=10); 
    ax.set_ylabel("$MSD_{||} \,\, (\AA^2)$",fontsize=10)
    plt.savefig(plt_nm+csvL[i][j][k].key[loc]+'.png', bbox_inches = 'tight',) 

def main():
    '''For a collection of data, get info from csv and then plot,
       usage: plot_vs_x.py csvStart nsep nlen niter sep1 sep2... len1 len2... 
                           iter1 iter2... ext datLoc'''
    csvname = sys.argv[1]; nsep = int(sys.argv[2]); nlen = int(sys.argv[3]);
    niter = int(sys.argv[4]); spS = 5;
    spE = spS + nsep; lnE = spE + nlen; itE = lnE + niter
    sep, ln, itr = [], [], []
    
    for i in range(spS, spE): sep.append(sys.argv[i])
    for i in range(spE, lnE): ln.append(int(sys.argv[i]))
    for i in range(lnE, itE): itr.append(sys.argv[i])
    
    ext = sys.argv[itE]; loc = int(sys.argv[itE+1])
    csvL = []; 
    for i in range(nsep):
        csL = []
        for j in range(nlen):
            cs = []
            for k in range(niter):
                csvn = csvname + str(sep[i])+"_"+str(ln[j])+"_"+str(itr[k])
                cs.append(CSVFile(csvn+ext))
            csL.append(cs)
        csvL.append(csL)
   #print("Ext {0}, loc {1} = {2}".format(ext, loc, cs[-1].key[loc]), sep, ln, itr)
    plot_scatter(csvname+"_", csvL, loc, sep, ln, itr)

if __name__=="__main__":
    main()

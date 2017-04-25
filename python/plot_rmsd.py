import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import *
MaxNLocator.default_params['nbins']=5

from csvfile import CSVFile

colorL = [[0,0,0], [0,0,1], [1, 0,0], [.93, .53, .18]]

def plot_scatter(plt_nm, csvL, loc, sep, ln, itr):
    '''Using data from a histogram, plot several'''
    f = plt.figure(1, figsize = (1.5, 1.5))
    ax, ct, leg = f.add_subplot(111), 0, []
    matplotlib.rcParams.update({'font.size': 4})
    ax.xaxis.set_major_locator(MaxNLocator())
    ax.yaxis.set_major_locator(MaxNLocator())

    for i in range(len(csvL)):
        for j in range(len(csvL[i])):
            for k in range(len(csvL[i][j])):
                print("This is plot val mean, std: {0} {1}".format(
                      np.mean(csvL[i][j][k].dat[loc]),
                      np.std(csvL[i][j][k].dat[loc])))
                m,b = np.polyfit(csvL[i][j][k].dat[0]*2./1000.,csvL[i][j][k].dat[loc],1)
                print("Linear fit slope {0}, intercept {1}".format(m/6.,b))
                ax.plot(range(len(csvL[i][j][k].dat[0])), 
                        csvL[i][j][k].dat[loc], color=colorL[ct])
                leg += [str(sep[i])+r'$\AA$ Sep, '+'L='+str(ln[j])+r'$\AA$']
                ct += 1
   #ax.set_xlim([-9,9])
    ax.set_xlim([0,600])
    ax.set_ylim([0,6000])
    ax.legend(leg, loc = 9, ncol = 1,
        columnspacing = 0.4,
        fontsize =  7 ,
        handletextpad = 0.2,
        handlelength = 1.3,
        borderaxespad = -0.9,
       #bbox_to_anchor = bbox
        )
    ax.set_xlabel("Time (ps)",fontsize=12)
    ax.set_ylabel("MSD ($\AA$)",fontsize=12)
    ax.yaxis.labelpad = -0.6; ax.xaxis.labelpad = -0.6
    plt.savefig(plt_nm+csvL[i][j][k].key[loc]+'.png', format='png',
                    bbox_inches = 'tight', dpi=300) 
    plt.close()

def main():
    '''For a collection of data, get info from csv and then plot,
       usage: plot_vs_x.py csvStart nsep nlen niter sep1 sep2... len1 len2... 
                           iter1 iter2... ext datLoc'''
    csvname = sys.argv[1]; nsep = int(sys.argv[2]); nlen = int(sys.argv[3]);
    niter = int(sys.argv[4]); spS = 5;
    spE = spS + nsep; lnE = spE + nlen; itE = lnE + niter
    sep, ln, itr = [], [], []
    
    for i in range(spS, spE): sep.append(int(sys.argv[i]))
    for i in range(spE, lnE): ln.append(int(sys.argv[i]))
    for i in range(lnE, itE): itr.append(int(sys.argv[i]))
    
    ext = sys.argv[itE]; loc = int(sys.argv[itE+1])
    csvL = []; csvname += (str(sep[0])+"_"+str(ln[0])+"_"+str(itr[0]))
    for i in range(nsep):
        csL = []
        for j in range(nlen):
            cs = []
            for k in range(niter):
                cs.append(CSVFile(csvname+ext))
            csL.append(cs)
        csvL.append(csL)
    print("Ext {0}, loc {1} = {2}".format(ext, loc, cs[-1].key[loc]))
    plot_scatter(csvname+"_", csvL, loc, sep, ln, itr)

if __name__=="__main__":
    main()

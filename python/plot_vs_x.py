import sys
import numpy as np
import matplotlib.pyplot as plt

from csvfile import CSVFile

colorL = [[0,0,0], [0,0,1], [1, 0,0], [.93, .53, .18]]

def plot_scatter(plt_nm, csvL, loc):
    '''Using data from a histogram, plot several'''
    f = plt.figure(1, figsize = (3.0, 3.0))
    ax, ct = f.add_subplot(111), 0

    for i in range(len(csvL)):
        for j in range(len(csvL[i])):
            for k in range(len(csvL[i][k])):
                shft = max(csvL[i][j][k].dat[0])/2.0
                ax.scatter(csvL[i][j][k].dat[0]-shft, 
                           csvL[i][j][k].dat[loc], color=colorL[ct])
                ct += 1
    ax.set_xlim(xrn); ax.set_ylim(yrn)
    plt.savefig(plt_nm+str(j)+'.png', format='png',                       
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
    print("This is nsep {0}, nlen {1}, niter {2}".format(nsep, nlen, niter))
    
    for i in range(spS, spE): sep.append(int(sys.argv[i]))
    for i in range(spE, lnE): ln.append(int(sys.argv[i]))
    for i in range(lnE, itE): itr.append(int(sys.argv[i]))
    print(sep,ln, itr)
    
    ext = sys.argv[itE]; loc = sys.argv[itE+1]
    print("Ext {0}, loc {1}".format(ext, loc))
    csvL = []
    for i in range(nsep):
        csL = []
        for j in range(nlen):
            cs = []
            for k in range(niter):
                cs.append(CSVFile(csvname+str(sep[i])+"_"+str(ln[j])+"_"
                                  +str(itr[k])+ext))
            csL.append(cs)
        csvL.append(csL)
    plot_scatter(csvname[:-3], csvL, loc)

if __name__=="__main__":
    main()

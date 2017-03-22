import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

from csvfile import CSVFile

colorL = [[0,0,0], [0,0,1], [1, 0,0], [.93, .53, .18]]

def plot_scatter(plt_nm, csvL, loc, sep, ln, itr):
    '''Using data from a histogram, plot several'''
    matplotlib.rcParams['font.size'] = 5
    f = plt.figure(1, figsize = (3.0, 3.0))
    f1= plt.figure(2, figsize = (3.0, 3.0))
    ax, ax1, ct, leg = f.add_subplot(111),f1.add_subplot(111), 0, []

    for i in range(len(csvL)):
        for j in range(len(csvL[i])):
            for k in range(len(csvL[i][j])):
                # distances = X, angle = Y
                X = csvL[i][j][k].dat[8]; Y = csvL[i][j][k].dat[loc]*np.sin(csvL[i][j][k].dat[loc])
                nbins = 20
                xbins = np.linspace(min(X), max(X), 7)
                print(min(X), max(X), xbins)
                ybins = np.linspace(min(Y), max(Y), nbins)
                heatmap, xedges, yedges = np.histogram2d(X, Y, bins=(xbins,ybins))
               #print(heatmap.shape)
                heatmap = heatmap.astype(float)
                for b in range(len(xedges)-1): ax1.plot(yedges[:-1], heatmap[b])
                extent = [min(X), max(X), min(Y), max(Y)]
               #ax.imshow(heatmap.T,aspect='auto', extent=extent)
                leg += [str(sep[i])+r'$\AA$ Sep, '+'L='+str(ln[j])+r'$\AA$']
                ct += 1
   #ax.set_xlim([-9,9])
    titl = leg[0]+" "+csvL[0][0][0].key[8]+" vs "+csvL[0][0][0].key[loc]
    ax.set_title(titl)
   #ax.legend(leg, loc = 9, ncol = 1,
   #    columnspacing = 0.4,
   #    fontsize =  7 ,
   #    handletextpad = 0.2,
   #    handlelength = 1.3,
   #    borderaxespad = -0.9,
   #   #bbox_to_anchor = bbox
   #    )
   #plt.show()
    plt.savefig(csvL[0][0][0].key[loc]+'.png', format='png',
                    bbox_inches = 'tight', dpi=200) 
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
    plot_scatter(ext[1:], csvL, loc, sep, ln, itr)

if __name__=="__main__":
    main()

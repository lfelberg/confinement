import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import *
MaxNLocator.default_params['nbins']=5

from csvfile import CSVFile
A2_TO_M2 = 1e-20
PS_TO_S  = 1e12
A2_PS_TO_M2_S = A2_TO_M2*PS_TO_S 

colorL = [
          [0,0,0], 
          [1, 0, 0], 
          [0.859, 0.078, 0.234], 
          [.93, .53, .18],
          [0.801, 0.586, 0.047],
          [0., 0.391, 0.], 
          [0.289, 0.461, 1.],
          [0.289, 0., 0.508], 
          [0.777, 0.441, 0.441],
          [0,0,1], 
          [0.777, 0.379, 0.078], 
          [0., 0.9297, 0.9297]
         ]

dens = {
         6   :   0.065,
         7   :   0.097,
         8   :   0.130,
         9   :   0.163,
         10  :   0.196,
         11  :   0.227,
         12  :   0.260,
         13  :   0.293,
         14  :   0.325,
         16  :   0.390,
         20  :   0.521,
         37  :   0.99,   # for bulk
}

def plot_scatter(plt_nm, csvL, loc, sep, ln, itr):
    '''Using data from a histogram, plot several'''
    f = plt.figure(1, figsize = (1.5, 1.5))
    ax, ct, leg = f.add_subplot(111), 0, []
    ax.xaxis.set_major_locator(MaxNLocator()); ax.yaxis.set_major_locator(MaxNLocator())

    for i in range(len(csvL)):
        for j in range(len(csvL[i])):
            for k in range(len(csvL[i][j])):
                m,b = np.polyfit(csvL[i][j][k].dat[0,150:300],
                                 csvL[i][j][k].dat[loc,150:300],1)
               #m,b = np.polyfit(csvL[i][j][k].dat[0,20:],
               #                 csvL[i][j][k].dat[loc,20:],1)
                print("{0:.3f},{1:.5f}".format(dens[sep[i]],m/4.*A2_PS_TO_M2_S*(1e9)))
                ax.plot(range(len(csvL[i][j][k].dat[0])), 
                        csvL[i][j][k].dat[loc], color=colorL[ct])
                leg += [r'$\rho_{{2D}}=${0:.2f}'.format(dens[sep[i]])]
                ct += 1
    ax.set_xlim([0,1000])
    ax.set_ylim([0,6000])
    bbox = [1.1, 0.95]
    ax.legend(leg, loc = 2, ncol = 1, columnspacing = 0.4,
        fontsize =  5 , handletextpad = 0.2, handlelength = 1.3,
        borderaxespad = -0.9, bbox_to_anchor = bbox   )
    ax.set_xlabel("Time (ps)",fontsize=10); 
    ax.set_ylabel("$MSD_{||} \,\, (\AA^2)$",fontsize=10)
    plt.savefig(plt_nm+csvL[i][j][k].key[loc]+'.png', bbox_inches = 'tight',) 

    fn = 'diffusion_coeff_2D.png'
    im = plt.imread(fn, format='png')
    xl = ax.get_xlim(); yl = ax.get_ylim()
   #newax = f.add_axes([0.135, 0.40, 0.48, 0.97], anchor='SW',) # for flex
    newax = f.add_axes([0.295, 0.23, 0.60, 0.97], anchor='SW',) # rigid
    newax.imshow(im, extent=[0, 1,0., 1.]) #, zorder=19)
    newax.axis('off')
    plt.savefig(plt_nm+csvL[i][j][k].key[loc]+'_insert.png',bbox_inches='tight',)
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
    print("Ext {0}, loc {1} = {2}".format(ext, loc, cs[-1].key[loc]), sep, ln, itr)
    plot_scatter(csvname+"_", csvL, loc, sep, ln, itr)

if __name__=="__main__":
    main()

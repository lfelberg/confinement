import sys
import numpy as np
import matplotlib.pyplot as plt

from csvfile import CSVFile

colorL = [[0,0,0], [0,0,1], [1, 0,0], [.93, .53, .18],
          [0.859, 0.078, 0.234], [0.801, 0.586, 0.047],                            
            [0., 0.391, 0.], [0.289, 0.461, 1.],                                   
            [0.289, 0., 0.508], [0.777, 0.441, 0.441],                             
            [0.777, 0.379, 0.078], [0., 0.9297, 0.9297]]

def plot_scatter(plt_nm, csvL, sep, ln):
    '''Using data from a histogram, plot several'''
    f = plt.figure(1, figsize = (3.0, 3.0))
    ax, ct, leg = f.add_subplot(111), 0, []

    for i in range(len(csvL)):
        for j in range(len(csvL[i])):
           ax.plot(csvL[i][j].dat[0], 
                   csvL[i][j].dat[1], color=colorL[ct])
           leg += [r"{0}$\AA$ sep".format(sep[i])]
           ct += 1
           #for k in range(len(csvL[i][j].dat[:-1])):
           #    ax.plot(csvL[i][j].dat[-1][2:], 
           #            csvL[i][j].dat[k][2:], color=colorL[ct])
           #    leg += [r"{0}$\AA$ sep, x={1:.1f}$\AA$".format(sep[i], 
           #                             float(csvL[i][j].key[k]))]
           #    ct += 1
   #ax.plot(csvL[i][j].dat[0], 
   #        np.repeat(864, len(csvL[i][j].dat[0])), color=colorL[ct])
    ax.legend(leg, loc = 9, ncol = 2,
        columnspacing = 0.4,
        fontsize =  4 ,
        handletextpad = 0.2,
        handlelength = 1.3,
        borderaxespad = -0.9,
        )
   #ax.set_ylim([0,0.5])
    plt.savefig(plt_nm+'.png', format='png',                       
                    bbox_inches = 'tight', dpi=300) 
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

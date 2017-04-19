import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from csvfile import CSVFile

def fourier(grcl):
    '''Fourier transform the g(r)'''
    ft = np.fft.fft(grcl.dat[2])
    f = plt.figure(1, figsize = (1.5, 1.5))
    ax, ct, leg = f.add_subplot(111), 0, []
    matplotlib.rcParams.update({'font.size': 8}) 
    ax.plot(np.arange(len(ft)),ft.real)
    plt.show()

def print_gr(x, grs, fname):
    '''Print distances to carbon wall in xyz like format'''
    print(grs.shape)
    f = open(fname, 'w'); 
    f.write("Bin,"); st = ""
    dimbins = HIS[:-1]+dr/2.
    for i in range(len(x)): st += "x{0:.3f},".format(x[i])
    f.write("{0}\n".format(st[:-1]))
    for i in range(len(dimbins)-1):
        st = ""
        for j in range(len(grs[i])): st += "{0:.5f},".format(grs[i][j])
        f.write("{0:.4f},{1}\n".format(dimbins[i], st[:-1]))
    f.close()

def main():
    ''' Given a list of pairs for g(r) calcs, do this as a function as 
        distance from wall, but only 2D and only for your wall side '''
    grname=sys.argv[1]; sep=sys.argv[2]; ln=sys.argv[3]; itr=sys.argv[4]

    nm = str(sep)+"_"+str(ln)+"_"+str(itr)
    grcl = CSVFile(grname)
    fourier(grcl)


if __name__=="__main__":
    main()

import sys
import re
import math
import numpy as np
import itertools
import matplotlib.pyplot as plt 

from scipy.integrate import simps

from csvfile import CSVFile
from volfile import VolFile
from util    import d_pbc

RMAX = 15.
NANGLES = 5
KB = 1.38064852e-23   # joules/kelvin
MOL = 6.022140857e23  # avogadro's number
J_TO_CAL = 1.0/4.1868 # 1 calorie/4.1868 Joules 
KB_CAL_MOL = KB * J_TO_CAL * MOL

def g_of_r(bns, dist, density = 1.0):
    ''' Given an array of distances, histogram to compute g(r) '''
    hist, bins = np.histogram(dist, bins = bns)
    upp, low = bns[1:], bns[:-1]
    ndens = float(len(dist)) / density
    gr = hist/(4.0/3.0*np.pi*(np.power(upp, 3.) - np.power(low,3.)))/ndens
    return gr, ndens

def g_of_r_multi(bns, dist):
    ''' Given an array of distances, histogram to compute g(r) 
        in multi dimensions, currently used only for orientation, so 
        we normalize in a different way, also going to assume each dim has same
        bins!
    '''
    ar = 1.0; dim = dist.shape[1] # number of dimensions for histogram
    bn_t = tuple([len(bns)-1 for x in range(dim)]) # nu bins from bns input
    rg_t = tuple([tuple([min(bns),max(bns)]) for x in range(dim)])
    hist, bins = np.histogramdd(dist, bins = bn_t, range = rg_t)
    for i in range(dim): ar *= (bins[i][1]-bins[i][0])
    ar *= np.sum(hist)
    return hist/ar

def trans_entropy(dists, vl):
    '''Compute the translational entropy from a list of pairwise distances
    '''
    ntimes, npairs = dists.shape[0], dists.shape[1]
    binsiz_r, grs, ndens = 0.03, [], np.zeros(ntimes)
    bns_r = np.arange( 2.5, RMAX, binsiz_r)
    radi_r  = binsiz_r/2.+bns_r[:-1] # center of each histogram bin
    nb_r = len(bns_r)-1

   #f = plt.figure(1, figsize = (3.0, 3.0))
   #ax = f.add_subplot(111)

    # calculating g(R) for each time snapshot
    grs = np.zeros((ntimes, nb_r))
    for t in range(ntimes):
        grs[t], ndens[t] = g_of_r(bns_r, dists[t], vl.get_vol_i(t))
       #ax.plot(radi_r, grs[t])
   #plt.show()

    gr = open("trans_gr.csv", "w")
    st = "bin,"
    for t in range(ntimes): st += ("time"+str(t)+",")
    gr.write(st[:-1]+"\n")
    for b in range(nb_r):
        st = "{0:.4f},".format(radi_r[b])
        for t in range(ntimes): st += "{0:.6f},".format(grs[t][b])
        gr.write(st[:-1]+"\n")
    gr.close()

    grs_avg = np.mean(grs, axis = 0)
    s_t_integrand = grs_avg*np.log(grs_avg) - grs_avg + 1.0
    ent_t = simps(s_t_integrand, radi_r)
    return -0.5 * ent_t * KB_CAL_MOL * np.mean(ndens)

def orien_order_entropy(order, dists, angles, vl):
    '''Compute the orientational entropy from a list of pairwise distances
       and angles to the order approximation
    '''
    ntimes, npairs, angle_combos = dists.shape[0], dists.shape[1], []
    binsiz_ra, binsiz_a, grs, nden = 0.10, 0.174533, [], np.zeros(ntimes)
    bns_ra = np.arange(2.5,RMAX,binsiz_ra); 
    bns_a = np.arange(0,np.pi+binsiz_a,binsiz_a)
    radi_a  = binsiz_a/2.+bns_a[:-1] # center of each histogram bin
    radi_ra = binsiz_ra/2.+bns_ra[:-1] # center of each histogram bin
    nb_a, nb_ra = len(radi_a), len(radi_ra) # Number of bins for each hist

    for subset in itertools.combinations(range(NANGLES), order):
        angle_combos.append(list(subset))
    ncombos = len(angle_combos)
    print(angle_combos, ncombos, len(angles))
    # calculating g(R) for each time snapshot
    dim_ang_gr = [ntimes, nb_ra, ncombos] + [nb_a for o in range(order)]
    grs = np.zeros(tuple(dim_ang_gr))
    grs_r = np.zeros((ntimes, nb_ra))
    sh_shell = np.zeros((nb_ra, ncombos))
    subset = np.zeros((order,npairs))
    
    for t in range(ntimes):
        r_bins = np.digitize(dists[t], bns_ra)
        grs_r[t], nden[t] = g_of_r(bns_ra, dists[t], vl.get_vol_i(t))
        for bn in range(nb_ra): #for each r bin, compute the g(angle)
             this_r = r_bins == bn
             if sum(this_r.astype(int)) == 0: sh_shell[bn] = np.zeros(ncombos)
             else:
                  # array for data to be histogrammed, shape = [nsamp, #ang]
                  an_dat = np.zeros((sum(this_r.astype(int)),order))
                  for an in range(ncombos):
                      an_rng = [t*NANGLES+i for i in angle_combos[an]]
                      for o in range(order): an_dat[:,o] = angles[an_rng[o], this_r]
                      grs[t][bn][an] = g_of_r_multi(bns_a, an_dat)
                  gr_mean = np.mean(grs[:,bn], axis = 0)
                  gr_mean[gr_mean==0] = 1.0
                  integ = gr_mean*np.log(gr_mean)
                  for i in range(order):  # integration over all dims of hist
                      integ = simps(integ, radi_a)/np.pi
                  sh_shell[bn] = -integ

    for bn in range(nb_ra): #for each r bin and angle combo
        for an in range(ncombos):
            ans, st = "", ""
            for i in angle_combos[an]: 
                ans += str(i)+"_"
                st += "bin"+str(i)+","
            f = open("angle_g_o{0}_bn{1:.4f}.csv".
                      format(ans[:-1],radi_ra[bn]), "w")
            for t in range(ntimes): st += ("time"+str(t)+",")
            f.write(st[:-1]+"\n")

            if order == 1:
                inds = np.mgrid[0:nb_a:1][:,np.newaxis]
                xy = radi_a[:,np.newaxis]
            if order == 2:
                inds = np.mgrid[0:nb_a:1,0:nb_a:1].reshape(2,-1).T
                xy = np.mgrid[0:radi_a[-1]+binsiz_a:binsiz_a,
                              0:radi_a[-1]+binsiz_a:binsiz_a].reshape(2,-1).T
            if order == 3:
                inds = np.mgrid[0:nb_a:1,0:nb_a:1,0:nb_a:1].reshape(3,-1).T
                xy = np.mgrid[0:radi_a[-1]+binsiz_a:binsiz_a,
                              0:radi_a[-1]+binsiz_a:binsiz_a,
                              0:radi_a[-1]+binsiz_a:binsiz_a].reshape(3,-1).T

            for ab in range(inds.shape[0]):
                st = ""
                for di in range(inds.shape[1]): 
                    st+="{0:.4f},".format(xy[ab][di])
                for t in range(ntimes): 
                    dat = grs[t][bn][an]
                    st += "{0:.6f},".format(dat[tuple(inds[ab])])
                f.write(st[:-1]+"\n")
            f.close()
            
    grs_avg = np.mean(grs_r, axis = 0)
    s_o_integrand = grs_avg[:,np.newaxis]*sh_shell
    radi_ra_ord = np.repeat(radi_ra[:,np.newaxis],ncombos,axis=1)
    ent_o = simps(s_o_integrand,radi_ra_ord,axis = 0)
    print(-0.5 * ent_o * KB_CAL_MOL * np.mean(nden), -0.5*sum(ent_o)*KB_CAL_MOL*np.mean(nden))
    return -0.5 * ent_o * KB_CAL_MOL * np.mean(nden)

def main():
    ''' Given a csv file with 5 angles, oxygen pair distance and dist from wall
        calculate translational and rotational entropy
    '''
    angname=sys.argv[1]; sep=sys.argv[2]; ln=sys.argv[3]; itr=sys.argv[4]
    ent_type = sys.argv[5]
    nm = str(sep)+"_"+str(ln)+"_"+str(itr)
    angC = CSVFile(angname)
    volC = VolFile("run"+nm+".vol")
    
    dis_loc = angC.find_keyword("dis")
    other_loc = angC.find_not_keyword("dis")

    if ent_type == "trans" or ent_type == "both":
        ent_t  = trans_entropy(angC.dat[dis_loc], volC)
        print("This is the translational entropy: {0:.4f}".format(ent_t))
   #ent_or = orien_entropy(angC.dat[dis_loc], angC.dat[other_loc], volC)
    if ent_type == "orien" or ent_type == "both":
        nord = int(sys.argv[6])
        if nord >= 1:
            ent_or_1 = orien_order_entropy(1,angC.dat[dis_loc],angC.dat[other_loc],volC)
        if nord >= 2:
            ent_or_2 = orien_order_entropy(2,angC.dat[dis_loc],angC.dat[other_loc],volC)
        if nord >= 3:
            ent_or_3 = orien_order_entropy(3,angC.dat[dis_loc],angC.dat[other_loc],volC)
        if nord >= 4:
            ent_or_4 = orien_order_entropy(4,angC.dat[dis_loc],angC.dat[other_loc],volC)

if __name__=="__main__":
    main()

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

RMAX = 6. 

# Angle list in order from csv
NANGLES = 5
ANG_NM = [ 'chi1', 'chi2', 'phi', 'the1', 'the2']
ANG_NORM = [ np.pi, np.pi, np.pi, 2.0, 2.0] # THIS IS FOR A RNG OF 0-pi for all

## Constants for entropy eval @ final, k_B and rho
KB = 1.38064852e-23   # joules/kelvin
MOL = 6.022140857e23  # avogadro's number
J_TO_CAL = 1.0/4.1868 # 1 calorie/4.1868 Joules 
KB_CAL_MOL = KB * J_TO_CAL * MOL
NUM_DENSITY_H2O = 1728./52534.8456042 # num of wat mols per A^3 (TIP4P EW 298K)

def g_of_r(bns, dist, density = 1.0):
    ''' Given an array of distances, histogram to compute g(r) '''
    hist, bins = np.histogram(dist, bins = bns)
    upp, low = bns[1:], bns[:-1]
    ndens = float(len(dist)) / density
    gr = hist/(4.0/3.0*np.pi*(np.power(upp, 3.) - np.power(low,3.)))/ndens
    return gr

def trans_entropy(dists, vl):
    '''Compute the translational entropy from a list of pairwise distances
    '''
    ntimes, npairs = dists.shape[0], dists.shape[1]
    binsiz_r, grs, ndens = 0.03, [], np.zeros(ntimes)
    bns_r = np.arange(0.0, RMAX, binsiz_r)
    radi_r  = binsiz_r/2.+bns_r[:-1] # center of each histogram bin
    nb_r = len(bns_r)-1

    # calculating g(R) for each time snapshot
    grs = np.zeros((ntimes, nb_r))
    for t in range(ntimes): grs[t] = g_of_r(bns_r, dists[t], vl[t])

    print_gr(ntimes,radi_r, grs) # printing for if needed later
    gr_av = np.mean(grs, axis = 0)
    nzer = gr_av != 0.0
    s_t_integrand = np.zeros(len(gr_av))
    s_t_integrand[nzer] = gr_av[nzer]*np.log(gr_av[nzer])-gr_av[nzer]+1.0
    s_t_integrand[gr_av == 0.0] = 1.0
    ent_t = simps(s_t_integrand*np.power(radi_r,2.0)*4.*np.pi, radi_r)
    return -0.5 * ent_t * KB_CAL_MOL * NUM_DENSITY_H2O

def trans_gr(gr_dat):
    '''Given a g(r), compute the translational entropy'''
    gr_av = np.mean(gr_dat[1:], axis = 0)
    nzer = gr_av != 0.0
    s_t_integrand = np.zeros(len(gr_av))
    s_t_integrand[nzer] = gr_av[nzer]*np.log(gr_av[nzer])-gr_av[nzer]+1.0
    s_t_integrand[gr_av == 0.0] = 1.0
    # integrate with r^2 * 4 * pi for volume
    ent_t = simps(s_t_integrand*np.power(gr_dat[0],2.0)*4.*np.pi, gr_dat[0])
    print(ent_t, KB_CAL_MOL, NUM_DENSITY_H2O)
    return -0.5 * ent_t * KB_CAL_MOL * NUM_DENSITY_H2O

def print_gr(ntimes, radi_r, grs):
    '''Print the gr for each timestep for use later if needed '''
    gr = open("trans_gr.csv", "w");  st = "bin,"
    for t in range(ntimes): st += ("time"+str(t)+",")
    gr.write(st[:-1]+"\n")
    for b in range(nb_r):
        st = "{0:.4f},".format(radi_r[b])
        for t in range(ntimes): st += "{0:.6f},".format(grs[t][b])
        gr.write(st[:-1]+"\n")
    gr.close()

def g_of_r_multi(bns_r, bns_a, dat, nfact, ang_no):
    ''' Given an array of distances, histogram to compute g(r) 
        in multi dimensions, currently used only for orientation, so 
        we normalize in a different way, also going to assume each dim has same
        bins!
    '''
    dim = dat.shape[1] - 1 # number of angles for histogram
    bn_t = tuple([len(bns_a)-1 for x in range(dim)]) # nu bins from bns input
    bn_t = tuple([len(bns_r)-1]) + bn_t
    rg_t = tuple([tuple([min(bns_a),max(bns_a)]) for x in range(dim)])
    rg_t = tuple([tuple([min(bns_r),max(bns_r)])]) + rg_t

    for i in range(len(ang_no)): # transforming theta vars, first col is dist
        if ang_no[i] == 4: dat[:,i+1] = np.pi - dat[:,i+1]
    hist, bins = np.histogramdd(dat, bins = bn_t, range = rg_t)

    sfact = 1.0
    for i in range(len(ang_no)): # transforming theta vars, first col is dist
        if ang_no[i] > 2: sfact = np.sin(bins[i+1][:-1]+(bins[i+1][1]-bins[i+1][0])/2.)

    # normalizing the histogram by angle nfact*binwid*nsamp
    for i in range(dim): nfact /= (bins[i+1][1]-bins[i+1][0])
    nfact /= np.sum(hist.astype(float), axis = 1)[:, np.newaxis]
    return hist.astype(float) * nfact / sfact

def orien_order_entropy(order, keys, dists, angles, vl):
    '''Compute the orientational entropy from a list of pairwise distances
       and angles to the order approximation
    '''
    ntimes, npairs, angle_combos = dists.shape[0],dists.shape[1],[]
   #binsiz_ra, binsiz_a, grs, nden = 0.10, 0.174533, [], np.zeros(ntimes)
    binsiz_ra, binsiz_a, grs, nden = 1.00, 0.174533, [], np.zeros(ntimes)
    bns_ra = np.arange(1.8,RMAX,binsiz_ra); 
    bns_a = np.arange(0,np.pi+binsiz_a,binsiz_a)
    radi_a  = binsiz_a/2.+bns_a[:-1] # center of each histogram bin
    radi_ra = binsiz_ra/2.+bns_ra[:-1] # center of each histogram bin
    nb_a, nb_ra = len(radi_a), len(radi_ra) # Number of bins for each hist

    # Making a list of all possible combinations of angle groups
    for subset in itertools.combinations(range(NANGLES), order):
        angle_combos.append(list(subset))
    ncombos = len(angle_combos); nfacts = np.ones(ncombos)
    ifacts = np.ones((ncombos, nb_a, order))
    for i in range(ncombos): 
        for j in range(order): 
            nfacts[i] *= ANG_NORM[angle_combos[i][j]]
            if angle_combos[i][j] > 2: # for theta angles, need to mult by sin(theta)
                ifacts[i,:,j] *= np.sin(radi_a)

    f, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2, sharex='col', sharey='row')

    # calculating g(R) for each time snapshot
    dim_ang_gr = [ntimes, ncombos, nb_ra] + [nb_a for o in range(order)]
    grs = np.zeros(tuple(dim_ang_gr)); grs_r = np.zeros((ntimes, nb_ra))
    sh_shell = np.zeros((nb_ra, ncombos)); subset = np.zeros((order,npairs))

    # For each timestep, calc the histogram for all combos of angles
    for t in range(ntimes):
        r_bins = np.digitize(dists[t], bns_ra)
        grs_r[t] = g_of_r(bns_ra, dists[t], vl[t])
        for bn in range(nb_ra): #for each r bin, compute the g(angle)
            this_r = r_bins == bn
            if sum(this_r.astype(int)) == 0: sh_shell[bn] = np.zeros(ncombos)
        an_dat = np.zeros((npairs, order+1))
        an_dat[:,0] = dists[t]
        for an in range(ncombos):
            an_rng = [i*ntimes+t for i in angle_combos[an]]
            # adding correct angles for time t to hist input
            for a in range(order): an_dat[:,a+1] = angles[an_rng[a],:].T
            grs[t][an] = g_of_r_multi(bns_ra, bns_a, an_dat, nfacts[an], angle_combos[an])

    print(ifacts, ifacts.shape)

    for bn in range(nb_ra): #for each r bin, integrate the g(angle)
        gr_mean = np.mean(grs[:,:,bn], axis = 0)
        ax1.plot(radi_a, gr_mean[0], label = str(radi_ra[bn]))
        ax2.plot(radi_a, gr_mean[1], label = str(radi_ra[bn]))
        ax3.plot(radi_a, gr_mean[2], label = str(radi_ra[bn]))
        ax4.plot(radi_a, gr_mean[3], label = str(radi_ra[bn]))
        ax5.plot(radi_a, gr_mean[4], label = str(radi_ra[bn]))
        ax1.legend(); ax2.legend();ax3.legend();ax4.legend();ax5.legend();
        gr_mean[gr_mean==0] = 1.0
        
        integ = gr_mean*np.log(gr_mean)
        print(integ.shape)
        for i in range(order):  # integration over all dims of hist
            integ = simps(integ*ifacts[:,:,i], radi_a)/nfacts
        sh_shell[bn] = integ

    print(gr_mean.shape)
    plt.show()
   #print_gang(nb_ra, nb_a, ncombos, ntimes, angle_combos, radi_ra, radi_a,
   #           binsiz_a, grs, order)
            
    grs_avg = np.mean(grs_r, axis = 0)
    s_o_integrand = grs_avg[:,np.newaxis]*sh_shell
    radi_ra_ord = np.repeat(radi_ra[:,np.newaxis],ncombos,axis=1)
    s_o_integrand *= (np.power(radi_ra_ord,2.0)*4.*np.pi)
    ent_o = simps(s_o_integrand,radi_ra_ord,axis = 0)
    print(ent_o, -0.5 * ent_o * KB_CAL_MOL * NUM_DENSITY_H2O,
          -0.5*sum(ent_o)*KB_CAL_MOL*NUM_DENSITY_H2O)
    return -0.5 * ent_o * KB_CAL_MOL * NUM_DENSITY_H2O

def print_gang(nb_ra, nb_a, ncombos, ntimes, angle_combos, radi_ra, radi_a,
               binsiz_a, grs, order):
    '''Plot the histogram for each group for each timestep for multi
       dimensional g(\omega) '''
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
    

def main():
    ''' Given a csv file with 5 angles, oxygen pair distance and dist from wall
        calculate translational and rotational entropy
    '''
    angname=sys.argv[1]; sep=sys.argv[2]; ln=sys.argv[3]; itr=sys.argv[4]
    ent_type = sys.argv[5]
    nm = str(sep)+"_"+str(ln)+"_"+str(itr)
    angC = CSVFile(angname)
    dis_loc = angC.find_keyword("dis")
    vol_loc = angC.find_keyword("vol")

    if ent_type == "trans" or ent_type == "both":
        if "gr" in angname:    ent_t = trans_gr(angC.dat) 
        else:  ent_t = trans_entropy(angC.dat[dis_loc], angC.dat[vol_loc,0])
        print("This is the translational entropy: {0:.7f}".format(ent_t))
    if ent_type == "orien" or ent_type == "both":
        nord = int(sys.argv[6])
        other_loc = angC.find_not_keyword("dis")
        oth_key = [angC.key[i] for i in other_loc]
        if nord >= 1:
            ent_or_1 = orien_order_entropy(1,oth_key,angC.dat[dis_loc],
                                           angC.dat[other_loc],angC.dat[vol_loc,0])
        if nord >= 2:
            ent_or_2 = orien_order_entropy(2,angC.dat[dis_loc],
                                           angC.dat[other_loc],angC.dat[vol_loc,0])
        if nord >= 3:
            ent_or_3 = orien_order_entropy(3,angC.dat[dis_loc],
                                           angC.dat[other_loc],angC.dat[vol_loc,0])
        if nord >= 4:
            ent_or_4 = orien_order_entropy(4,angC.dat[dis_loc],
                                           angC.dat[other_loc],angC.dat[vol_loc,0])

if __name__=="__main__":
    main()

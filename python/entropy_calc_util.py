import sys
import re
import math
import numpy as np
import itertools
import matplotlib
import matplotlib.pyplot as plt 
from scipy.integrate import simps

from csvfile import CSVFile
from volfile import VolFile
from util    import d_pbc

RMAX = 12. 

# Angle list in order from csv
NANGLES = 5
ANG_NM = [r'$\chi_1$',r'$\chi_2$',r'$\phi$',r'$\theta_1$',r'$\theta_2$']
ANG_NORM = [ np.pi, np.pi, np.pi, 2.0, 2.0] # THIS IS FOR A RNG OF 0-pi for all

## Constants for entropy eval @ final, k_B and rho
KB = 1.38064852e-23   # joules/kelvin
MOL = 6.022140857e23  # avogadro's number
J_TO_CAL = 1.0/4.1868 # 1 calorie/4.1868 Joules 
KB_CAL_MOL = KB * J_TO_CAL * MOL
NUM_DENSITY_H2O = 1728./52534.8456042 # num of wat mols per A^3 (TIP4P EW 298K)


class STrans:
    def __init__(self, fname, dim):
        self.dim = dim   


def g_of_r(bns, dist, dim = 3, density = 1.0):
    ''' Given an array of distances, histogram to compute g(r) '''
    hist, bins = np.histogram(dist, bins = bns)
    upp, low = bns[1:], bns[:-1]
    ndens = float(len(dist)) / density
    if dim == 3: nfact = 4.0/3.0*np.pi*(np.power(upp, 3.) - np.power(low,3.)) 
    else:        nfact = np.pi*(np.power(upp, 2.) - np.power(low,2.))
    gr = hist/(nfact*ndens)
    return gr

def trans_entropy(dists, vl, dim = 3):
    '''Compute the translational entropy from a list of pairwise distances
    '''
    print("This is my dat ", dists.shape, vl)
    ntimes, npairs = dists.shape[0], dists.shape[1]
    binsiz_r, grs, ndens = 0.03, [], np.zeros(ntimes)
    bns_r = np.arange(0.0, RMAX, binsiz_r)
    radi_r  = binsiz_r/2.+bns_r[:-1] # center of each histogram bin
    nb_r = len(bns_r)-1

    # calculating g(R) for each time snapshot
    grs = np.zeros((ntimes, nb_r))
    for t in range(ntimes): grs[t] = g_of_r(bns_r, dists[t], dim, vl[t])

    print_gr(ntimes, nb_r, radi_r, grs) # printing for if needed later
    return integ_rg(radi_r, np.mean(grs,axis=0), dim)

def trans_gr(gr_dat, dim = 3):
    '''Given a g(r), compute the translational entropy'''
    radi_r = gr_dat[0]; gr_av = np.mean(gr_dat[1:], axis = 0)
    return integ_rg(gr_dat[0],np.mean(gr_dat[1:],axis=0), dim)

def integ_rg(radi_r, gr_av, dim):
    ''' Given distances and gr, integrate to calc entropy'''
    plot_gr(radi_r, gr_av)
    nzer = gr_av != 0.0
    s_t_integrand = np.zeros(len(gr_av))
    s_t_integrand[nzer] = gr_av[nzer]*np.log(gr_av[nzer])-gr_av[nzer]+1.0
    s_t_integrand[gr_av == 0.0] = 1.0
    # integrate with 4pi*r^2 for volume, 2*pi*r for area
    if dim == 3: r_int = np.power(radi_r,2.0)*4.*np.pi
    else:        r_int = radi_r*2.*np.pi
    ent_t = simps(s_t_integrand*r_int, radi_r)
    print("Etrans before conv {0}, conv {1}".format(ent_t, 
          -0.5 * KB_CAL_MOL * NUM_DENSITY_H2O))
    return -0.5 * ent_t * KB_CAL_MOL * NUM_DENSITY_H2O

def plot_gr(radi_r, gr_av):
    '''Plot the computed gr'''
    matplotlib.rcParams.update({'font.size': 8})
    f = plt.figure(1, figsize = (1.5, 1.5))
    ax = f.add_subplot(111)
    ax.plot(radi_r, gr_av)
    ax.set_xlim(0, RMAX); #ax.set_ylim(0, 3);
    ax.yaxis.labelpad = -0.6; ax.xaxis.labelpad = -0.6
    ax.set_ylabel("g(R)",fontsize=10)
    ax.set_xlabel("Distance ($\AA$)",fontsize=10)
    plt.savefig('g_r.png',format='png',
                bbox_inches='tight',dpi=300)
    plt.close()

def print_gr(ntimes, nb_r, radi_r, grs):
    '''Print the gr for each timestep for use later if needed '''
    gr = open("trans_gr.csv", "w");  st = "bin,"
    for t in range(ntimes): st += ("time"+str(t)+",")
    gr.write(st[:-1]+"\n")
    for b in range(nb_r):
        st = "{0:.4f},".format(radi_r[b])
        for t in range(ntimes): st += "{0:.6f},".format(grs[t][b])
        gr.write(st[:-1]+"\n")
    gr.close()

def integ_fact(shape, angle_nos, angle_rng):
    ''' Given a list of angles and an array shape, will create an array of
        shape = shape, to scale data by. Assumes that d1 = distance,
        and the other dimensions are angles 
    '''
    scale_factor = np.ones(shape)
    for i in range(len(angle_nos)): # Xfrm theta vars, first col is dist
        if angle_nos[i] > 2: 
            sft = np.sin(angle_rng)
            for j in range(len(angle_nos)):
                if i > j:   sft = np.expand_dims(sft, axis = 0)
                elif i < j: sft = np.expand_dims(sft, axis = -1)
            sft = np.expand_dims(sft, axis = 0)
            scale_factor *= sft
    return scale_factor

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
    
    sfact = integ_fact(hist.shape, ang_no,
                       bins[i+1][:-1]+(bins[i+1][1]-bins[i+1][0])/2.)
    # normalizing the histogram by angle nfact*binwid*nsamp
    r_bins_tot = hist.astype(float)
    for i in range(dim): 
        nfact /= (bins[i+1][1]-bins[i+1][0])
        r_bins_tot = np.sum(r_bins_tot, axis = -1)
    for i in range(len(ang_no)): r_bins_tot=np.expand_dims(r_bins_tot,axis=-1)
    return hist.astype(float) * nfact / sfact, r_bins_tot

def ang_cmbs(order):
    '''Given the desired order of entropy calc, return array of all combos'''
    angle_combos = []
    for subset in itertools.combinations(range(NANGLES), order):
        angle_combos.append(list(subset))
    return angle_combos, len(angle_combos)

def orien_order_entropy(order, keys, dists, angles, vl, dim):
    '''Compute the orientational entropy from a list of pairwise distances
       and angles to the order approximation
    '''
    ntimes, npairs = dists.shape[0],dists.shape[1]
    binsiz_ra, binsiz_a, grs, nden = 0.10, 0.174533, [], np.zeros(ntimes)
   #binsiz_ra, binsiz_a, grs, nden = 1.10, 0.7853981634, [], np.zeros(ntimes)
    bns_ra = np.arange(0.0,RMAX,binsiz_ra); 
    bns_a = np.arange(0,np.pi+binsiz_a,binsiz_a)
    radi_a  = binsiz_a/2.+bns_a[:-1] # center of each histogram bin
    radi_ra = binsiz_ra/2.+bns_ra[:-1] # center of each histogram bin
    nb_a, nb_ra = len(radi_a), len(radi_ra) # Number of bins for each hist

    angle_combos, ncombos = ang_cmbs(order); nfacts = np.ones(ncombos)
    hist_shape = tuple([nb_ra]) + tuple([nb_a for x in range(order)])
    ifacts = np.ones(tuple([ncombos])+hist_shape); 
    ntot = np.zeros(tuple([ncombos,nb_ra])+tuple([1 for x in range(order)]))
    for i in range(ncombos): 
        ifacts[i] = integ_fact(hist_shape, angle_combos[i], radi_a)
        for j in range(order): 
            nfacts[i] *= ANG_NORM[angle_combos[i][j]]

    # calculating g(R) for each time snapshot
    dim_ang_gr = [ntimes+1, ncombos, nb_ra] + [nb_a for o in range(order)]
    grs = np.zeros(tuple(dim_ang_gr)); grs_r = np.zeros((ntimes, nb_ra))
    sh_shell = np.zeros((nb_ra, ncombos))

    # For each timestep, calc the histogram for all combos of angles
    for t in range(ntimes):
        r_bins = np.digitize(dists[t], bns_ra)
        grs_r[t] = g_of_r(bns_ra, dists[t], dim, vl[t])
        for bn in range(nb_ra): #for each r bin, compute the g(angle)
            this_r = r_bins == bn
            if sum(this_r.astype(int)) == 0: sh_shell[bn] = np.zeros(ncombos)
        an_dat = np.zeros((npairs, order+1))
        an_dat[:,0] = dists[t]
        for an in range(ncombos):
            an_rng = [i*ntimes+t for i in angle_combos[an]]
            # adding correct angles for time t to hist input
            for a in range(order): an_dat[:,a+1] = angles[an_rng[a],:].T
            gr_one, nhis = g_of_r_multi(bns_ra, bns_a, an_dat, nfacts[an], 
                                        angle_combos[an])
            hi_nz = np.copy(nhis); hi_nz[hi_nz == 0.] = 1.0
            grs[t][an] = gr_one/hi_nz; grs[-1][an] += gr_one; ntot[an] += nhis

    ntot[ntot == 0.0] = 1.0 # finding empty bins, setting to 1. for divide
    gr_mean = grs[-1]/ntot 

    if order < 3: plot_g_ang(order,ncombos,angle_combos,nb_ra,
                             radi_a,radi_ra,gr_mean) # plotting 1/2 ord angs
    print_gang(nb_ra, nb_a, ncombos, ntimes, angle_combos, radi_ra, radi_a,
               binsiz_a, grs, ntot, order)
    return integ_angle(order, ncombos, nb_ra, nfacts, ifacts, sh_shell, 
                       radi_ra, grs_r, radi_a, ntot, gr_mean)

def integ_angle(order, ncombos, nb_ra, nfacts, ifacts, sh_shell, radi_ra, grs_r, radi_a, ntot, gr_mean):
    ''' Integrate g(angle)*g(R) to calc S_orient'''
    gr_mean[gr_mean==0] = 1.0; st = "bins: "; tt = ntot[0].flatten()
    for bn in range(nb_ra): #for each r bin, integrate the g(angle)
        st += "{0} ".format(tt[bn])
        integ = gr_mean[:,bn]*np.log(gr_mean[:,bn])*ifacts[:,bn]
        for i in range(order):  # integration over all dims of hist
            integ = simps(integ, radi_a)
        sh_shell[bn] = integ/nfacts

    grs_avg = np.mean(grs_r, axis = 0)
    s_o_integrand = grs_avg[:,np.newaxis]*sh_shell
    radi_ra_ord = np.repeat(radi_ra[:,np.newaxis],ncombos,axis=1)
    s_o_integrand *= (np.power(radi_ra_ord,2.0)*4.*np.pi)
    ent_o = simps(s_o_integrand,radi_ra_ord,axis = 0)
    print(st); print(ent_o, -0.5 * ent_o * KB_CAL_MOL * NUM_DENSITY_H2O,
          -0.5*sum(ent_o)*KB_CAL_MOL*NUM_DENSITY_H2O)
    return -0.5 * ent_o * KB_CAL_MOL * NUM_DENSITY_H2O

def orien_gang(order, dirs, gr_dat):
    ''' Given already computed g(r) and list of angle files, compute 
        the g(angle) of specified order'''
    ncombos, angle_combos = ang_cmbs(order) 
    

def plot_g_ang(order,ncombos,angle_combos,nb_ra,radi_a,radi_ra,gr_mean):
    '''Method for plotting the g(angle) distributions for 1&2nd order'''
    matplotlib.rcParams.update({'font.size': 4})
    if order == 1: 
        f,axes=plt.subplots(5,1,sharex='col',sharey='row',figsize=(1.3,3.5))
    else:
        f,axes=plt.subplots(5,2,sharex='col',sharey='row',figsize=(2,5)) 
        print(axes.shape, gr_mean.shape)
        plt.setp(axes.flat, aspect=1.0, adjustable='box-forced')
    axes = axes.flatten()
    bin_rng = [ 27, 30, 43, 65]

    for bn in range(nb_ra): #for each r bin, integrate the g(angle)
        if (bn in bin_rng and order == 1):
            for j in range(ncombos):
                axes[j].plot(radi_a,gr_mean[j][bn],label=str(radi_ra[bn]))
        if (bn == 27 and order == 2):
            X, Y = np.meshgrid(radi_a, radi_a); c = []
            for j in range(ncombos):
                c.append(axes[j].contour(X, Y, gr_mean[j][bn].T))
                axes[j].text(.5,.8,ANG_NM[angle_combos[j][0]]+","+
                             ANG_NM[angle_combos[j][1]],
                             horizontalalignment='center',
                             transform=axes[j].transAxes)
                plt.clabel(c[j], inline=1, fontsize= 3);
    max_tc = 4
    if order == 1:
        axes[0].legend(ncol=len(bin_rng),columnspacing=-0.1,labelspacing=-1.95,
                       borderaxespad=-1.9,handlelength=1.8,fontsize=4)
        axes[-1].set_xlabel("Angle (radians)",fontsize=6)
        axes[-1].xaxis.labelpad = -1
        for j in range(ncombos):
            axes[j].xaxis.set_major_locator(plt.MaxNLocator(max_tc))
            axes[j].yaxis.set_major_locator(plt.MaxNLocator(max_tc))
            axes[j].tick_params(axis='both', which='major', pad= 1.3)
            axes[j].set_xlim(0, np.pi);axes[j].yaxis.labelpad = -1
            axes[j].set_ylim(0, 3);
            axes[j].set_ylabel("g("+ANG_NM[angle_combos[j][0]]+")",fontsize=6)
            axes[j].text(.7,.6,ANG_NM[angle_combos[j][0]], 
                         horizontalalignment='center',
                         transform=axes[j].transAxes,fontsize=10);
    plt.subplots_adjust(wspace=0.14, hspace=0.14)
    plt.savefig("gangle_"+str(order)+'.png',format='png',
                bbox_inches='tight',dpi=300)

def grid_for_print(nb_a, radi_a, binsiz_a, order):
    '''Grid for a varying dimensional g(angle)'''
    if order == 1:
        inds = np.mgrid[0:nb_a:1][:,np.newaxis]
        xy = radi_a[:,np.newaxis]
    if order == 2:
        inds = np.mgrid[0:nb_a:1,0:nb_a:1].reshape(2,-1).T
        xy = np.mgrid[radi_a[0]:radi_a[-1]+binsiz_a:binsiz_a,
                      radi_a[0]:radi_a[-1]+binsiz_a:binsiz_a].reshape(2,-1).T
    if order == 3:
        inds=np.mgrid[0:nb_a:1,0:nb_a:1,0:nb_a:1].reshape(3,-1).T
        xy = np.mgrid[radi_a[0]:radi_a[-1]+binsiz_a:binsiz_a,
                      radi_a[0]:radi_a[-1]+binsiz_a:binsiz_a,
                      radi_a[0]:radi_a[-1]+binsiz_a:binsiz_a].reshape(3,-1).T
    if order == 4:
        inds=np.mgrid[0:nb_a:1,0:nb_a:1,0:nb_a:1,0:nb_a:1].reshape(4,-1).T
        xy = np.mgrid[radi_a[0]:radi_a[-1]+binsiz_a:binsiz_a,
                      radi_a[0]:radi_a[-1]+binsiz_a:binsiz_a,
                      radi_a[0]:radi_a[-1]+binsiz_a:binsiz_a,
                      radi_a[0]:radi_a[-1]+binsiz_a:binsiz_a].reshape(4,-1).T
    if order == 5:
        inds=np.mgrid[0:nb_a:1,0:nb_a:1,0:nb_a:1,0:nb_a:1,
                      0:nb_a:1].reshape(5,-1).T
        xy = np.mgrid[radi_a[0]:radi_a[-1]+binsiz_a:binsiz_a,
                      radi_a[0]:radi_a[-1]+binsiz_a:binsiz_a,
                      radi_a[0]:radi_a[-1]+binsiz_a:binsiz_a,
                      radi_a[0]:radi_a[-1]+binsiz_a:binsiz_a,
                      radi_a[0]:radi_a[-1]+binsiz_a:binsiz_a].reshape(5,-1).T
    return inds, xy

def print_gang(nb_ra, nb_a, ncombos, ntimes, angle_combos, radi_ra, radi_a,
               binsiz_a, grs, ntot, order):
    '''Plot the histogram for each group for each timestep for multi
       dimensional g(\omega) '''
    his_ct = open("angle_g_bin_ct.csv","w"); tt = ntot[0].flatten()
    his_ct.write("bin,ct\n")
    for i in range(len(tt)): 
        his_ct.write("{0:.3f},{1}\n".format(radi_ra[i],tt[i]))
    his_ct.close()

    inds,xy = grid_for_print(nb_a, radi_a, binsiz_a, order)

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

            for ab in range(inds.shape[0]):
                st = ""
                for di in range(inds.shape[1]): 
                    st+="{0:.4f},".format(xy[ab][di])
                for t in range(ntimes+1): 
                    dat = grs[t][an][bn]
                    st += "{0:.6f},".format(dat[tuple(inds[ab])])
                f.write(st[:-1]+"\n")
            f.close()


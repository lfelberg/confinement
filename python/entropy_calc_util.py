import sys
import re
import math
import numpy as np
import itertools
import matplotlib
import matplotlib.pyplot as plt 
from scipy.integrate import simps

from csvfile import CSVFile
from util    import d_pbc

RMAX = 9.

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
NUM_DENSITY_H2O_2D = 32.05/1392.4178 # num of wat mols per A^2 (TIP4P EW 298K)

class STrans:
    '''Class with variables and methods for the calc of translational entropy
    ''' 
    def __init__(self, ntimes, npairs, bnr = 0.03, dim = 3):
        '''Initialize vars for translational entropy'''
        self.dim = dim # dimension of g(r)
        self.ntimes = ntimes # number of snapshots to average across
        self.npairs = npairs # number of water pairs considered
        self.binsiz_r = bnr  # size of g(r) histogram bin
        self.init_g_bins()   # initialize bins for g(r)
        if dim == 3: self.rho = NUM_DENSITY_H2O
        else:        self.rho = NUM_DENSITY_H2O_2D

    def init_g_bins(self):
        '''Method to initialize g(val) things, like the number of bins, 
           bin size etc'''
        self.grs, self.ndens = [], np.zeros(self.ntimes)
        self.bns_r = np.arange(0.0, RMAX, self.binsiz_r)
        self.radi_r  = self.binsiz_r/2.+self.bns_r[:-1] # centr of hist bins
        self.nb_r = len(self.radi_r) # number of bins
        self.grs = np.zeros((self.ntimes, self.nb_r))

    def trans_entropy(self, dists, vl):
        '''Given a list of pair separations, compute the g(r), print out
           and integrate to get translational entropy'''
        self.many_g_of_r(dists, vl)
        return self.integ_rg(np.mean(self.grs,axis=0))

    def trans_gr(self, gr_dat, gr_dir = []):
        '''Given a g(r), compute the translational entropy'''
        if gr_dir != []: 
            nb = -1
            for i in range(len(gr_dir)): 
                grCSV = CSVFile("iter_"+str(gr_dir[i])+"/trans_gr_"
                                +str(gr_dir[i])+"_0.03.csv")
                radi_r = grCSV.dat[0]; print("This shape",grCSV.dat.shape)
                if nb == -1: 
                    nb = len(grCSV.dat[1])
                    gr_dat = np.mean(grCSV.dat[1:], axis = 0)
                else:  gr_dat += np.mean(grCSV.dat[1:], axis = 0)
            gr_av = gr_dat/float(len(gr_dir))
        else:
            radi_r = gr_dat[0]; gr_av = np.mean(gr_dat[1:], axis = 0)
        return self.integ_rg(gr_av,radi_r)

    def many_g_of_r(self, dists, vl):
        '''Given a list of pair separations, compute the g(r), print out'''
        for t in range(self.ntimes): 
            self.grs[t] = self.g_of_r(dists[t], vl[t])
        self.print_gr() # printing for if needed later

    def g_of_r(self, dist, density = 1.0):
        ''' Given an array of distances, histogram to compute g(r) '''
        hist, bins = np.histogram(dist, bins = self.bns_r)
        upp, low = self.bns_r[1:], self.bns_r[:-1]
        ndens = float(len(dist)) / density
        if self.dim == 3: nfact=4.0/3.0*np.pi*(np.power(upp,3.)-np.power(low,3.)) 
        else:        nfact = np.pi*(np.power(upp, 2.) - np.power(low,2.))
        return hist/(nfact*ndens)

    def integ_rg(self, gr_av, rad = []):
        ''' Given distances and gr, integrate to calc entropy'''
        if rad != []: self.radi_r = rad
        nzer = gr_av != 0.0;self.plot_gr(gr_av);
        
        s_t_integrand = np.zeros(len(gr_av))
        s_t_integrand[nzer] = gr_av[nzer]*np.log(gr_av[nzer])-gr_av[nzer]+1.0
        s_t_integrand[gr_av == 0.0] = 1.0
        # integrate with 4pi*r^2 for volume, 2*pi*r for area
        if self.dim == 3: r_int = np.power(self.radi_r,2.0)*4.*np.pi
        else:        r_int = self.radi_r*2.*np.pi
        ent_t = simps(s_t_integrand*r_int, self.radi_r)
        print("Etrans before conv {0}, conv {1}".format(ent_t, 
              -0.5 * KB_CAL_MOL * self.rho))
        return -0.5 * ent_t * KB_CAL_MOL * self.rho

    def plot_gr(self, gr_av):
        '''Plot the computed gr'''
        matplotlib.rcParams.update({'font.size': 8})
        f = plt.figure(1, figsize = (1.5, 1.5))
        ax = f.add_subplot(111)
        print(self.radi_r.shape, gr_av.shape)
        ax.plot(self.radi_r, gr_av)
        ax.set_xlim(0, RMAX); #ax.set_ylim(0, 3);
        ax.yaxis.labelpad = -0.6; ax.xaxis.labelpad = -0.6
        ax.set_ylabel("g(R)",fontsize=10)
        ax.set_xlabel("Distance ($\AA$)",fontsize=10)
        plt.savefig('g_r.png',format='png', bbox_inches='tight',dpi=300)
        plt.close()

    def print_gr(self):
        '''Print the gr for each timestep for use later if needed '''
        gr = open("trans_gr_{0:.2f}.csv".format(self.binsiz_r),"w"); st="bin,"
        for t in range(self.ntimes): st += ("time"+str(t)+",")
        gr.write(st[:-1]+"\n")
        for b in range(self.nb_r):
            st = "{0:.4f},".format(self.radi_r[b])
            for t in range(self.ntimes): st+="{0:.6f},".format(self.grs[t][b])
            gr.write(st[:-1]+"\n")
        gr.close()

class SOrien:
    '''Class with variables and methods for the calc of approximation
       of orientational entropy
    ''' 
    def __init__(self, ntimes, npairs, order, bnr = 0.10, dim = 3):
        '''Initialize vars for orientational entropy'''
        self.dim = dim # dimension of g(r)
        self.ntimes = ntimes # number of snapshots to average across
        self.npairs = npairs # number of water pairs considered
        self.order = order   # the order of the approximation 1-5
        self.binsiz_r = bnr  # size of g(r) histogram bin
        self.grC = STrans(ntimes, npairs, bnr, dim)
        self.init_g_bins()   # initialize bins for g(r)

    def init_g_bins(self):
        '''Method to initialize g(val) things, like the number of bins, 
           bin size etc'''
        self.grs,self.ndens,self.binsiz_a = [],np.zeros(self.ntimes),0.174533
        self.bns_a = np.arange(0,np.pi+self.binsiz_a,self.binsiz_a)
        self.radi_a  = self.binsiz_a/2.+self.bns_a[:-1] # centr of hist bins
        self.nb_a = len(self.radi_a) # number of bins

        self.ang_cmbs(); self.nfacts = np.ones(self.ncombos)
        self.ang_shape = tuple([self.nb_a for x in range(self.order)])
        self.hist_shape = tuple([self.get_rbn_ct()]) + self.ang_shape
        self.ntot = np.zeros(tuple([self.ncombos,self.get_rbn_ct()])+
                             tuple([1 for x in range(self.order)]))
        self.integ_fact()

    def orien_order_entropy(self, dists, angles, vl):
        '''Compute orientational entropy from a list of pairwise distances
           and angles to the order approximation
        '''
        self.dim_ang_gr = tuple([self.ntimes+1,self.ncombos])+self.hist_shape
        self.gs=np.zeros(self.dim_ang_gr)
        sh_shell = np.zeros((self.get_rbn_ct(), self.ncombos))

        self.grC.many_g_of_r(dists, vl)  # computing the g(r)
        # For each timestep, calc the histogram for all combos of angles
        for t in range(self.ntimes):
            an_dat=np.zeros((self.npairs,self.order+1)); an_dat[:,0]=dists[t]
            for an in range(self.ncombos):
                an_rng = [i*self.ntimes+t for i in self.angle_combos[an]]
                # adding correct angles for time t to hist input
                for a in range(self.order): an_dat[:,a+1]=angles[an_rng[a],:].T
                g_one, nhis = self.g_of_multi(an_dat,an)
                hi_nz = np.copy(nhis); hi_nz[hi_nz == 0.] = 1.0
                self.gs[t][an]=g_one/hi_nz; 
                self.gs[-1][an]+=g_one; self.ntot[an]+=nhis

        self.print_gang()
        self.ntot[self.ntot==0.0]=1. # empty bins, setting to 1. for divide
        g_mean = self.gs[-1]/self.ntot 

       #if self.order < 3: self.plot_g_ang(g_mean) # plotting 1/2 ord angs
        return self.integ_angle(sh_shell, g_mean)

    def orien_gang(self, dirs, gr_dat):
        ''' Given already computed g(r) and list of angle files, compute 
            the g(angle) of specified order'''
        # getting histogram tallies    
        nct = np.zeros((len(dirs), self.get_rbn_ct())); ndt = -1
        for i in range(len(dirs)): 
            ctCSV = CSVFile("iter_"+str(dirs[i]+"/angle_g_bin_ct.csv"))
            nct[i] = ctCSV.dat[1]
            gr = CSVFile("iter_"+str(dirs[i])+"/trans37_37_"+str(dirs[i])
                         +"_gr.0.1.csv")
            if ndt == -1: 
                ndt = len(gr.dat[1:]);
                gr_dat = np.zeros((len(dirs)*ndt,len(gr.dat[1])))
            gr_dat[i*ndt:(i+1)*ndt] = gr.dat[1:]
        nct[nct == 50.0]=0.0; nct[nct==1.0]=0.0 # THIS WILL NEED TO BE CHANGED
        nct = np.sum(nct,axis = 0)[np.newaxis]
        for i in range(self.order): nct = np.expand_dims(nct, axis = -1)
        nct[nct == 0.0] = 1.0 

        self.grC.grs = gr_dat[1:]  # computing the g(r)
       #self.grC.plot_gr(np.mean(gr_dat[1:],axis=0))
        gang = np.zeros(tuple([self.ncombos,len(dirs)])+self.hist_shape)
        for an in range(self.ncombos):
            an_nm = ""
            for j in self.angle_combos[an]: an_nm += (str(j)+"_")
            for bn in range(self.get_rbn_ct()):
                bn_nm = "bn{0:.4f}".format(self.get_r_rad()[bn])
                for i in range(len(dirs)):
                    fname = "iter_"+str(dirs[i])+"/angle_g_o"+an_nm+bn_nm+".csv"
                    dat = self.get_g_ang(fname)
                    gang[an][i][bn] = dat.reshape(self.ang_shape)

        g_an = np.sum(gang, axis = 1)/nct
       #if self.order < 3: self.plot_g_ang(g_an) # plotting 1/2 ord angs
        return self.integ_angle(np.zeros((self.get_rbn_ct(),self.ncombos)),g_an)

    def g_of_multi(self, dat, an):
        ''' Given an array of distances, histogram to compute g(r) 
            in multi dimensions, currently used only for orientation, so 
            we normalize in a different way, going to assume each dim has same
            bins!
        '''
        dim,ba,br = self.order,self.bns_a,self.get_rbn() # No. angles for hist
        bn_t = tuple([len(br)-1]) + tuple([len(ba)-1 for x in range(dim)]) #bns
        rg_t = tuple([tuple([min(ba),max(ba)]) for x in range(dim)])
        rg_t = tuple([tuple([min(br),max(br)])]) + rg_t

        for i in range(self.order): # XFRM theta vars, first col is dist
            if self.angle_combos[an][i] == 4: dat[:,i+1] = np.pi - dat[:,i+1]
        hist, bins = np.histogramdd(dat, bins = bn_t, range = rg_t)
        
        # normalizing the histogram by angle nfact*binwid*nsamp
        r_bins_tot = hist.astype(float);  scl = self.nfacts[an]
        for i in range(dim): 
            scl /= (bins[i+1][1]-bins[i+1][0])
            r_bins_tot = np.sum(r_bins_tot, axis = -1)
        for i in range(self.order): 
            r_bins_tot=np.expand_dims(r_bins_tot,axis=-1)
        return hist.astype(float) * scl / self.ifacts[an], r_bins_tot

    def get_g_ang(self, fname):
        '''From each csv file, get the data'''
        f = open(fname,'r').readlines()
        dat = np.zeros(len(f)-1)
        for l in range(1,len(f)):
            tmp = f[l].split(",")
            dat[l-1] = tmp[-1]
        return dat

    def integ_angle(self, sh_shell, g_mean):
        ''' Integrate g(angle)*g(R) to calc S_orient'''
        g_mean[g_mean==0] = 1.0; st = "bins: "; tt = self.ntot[0].flatten()
        for bn in range(self.get_rbn_ct()): #for each rbn, integ g(angle)
            st += "{0} ".format(tt[bn])
            integ = g_mean[:,bn]*np.log(g_mean[:,bn])*self.ifacts[:,bn]
            for i in range(self.order):  # integration over all dims of hist
                integ = simps(integ, self.radi_a)
            sh_shell[bn] = integ/self.nfacts

        grs_avg = np.mean(self.grC.grs, axis = 0)
        s_o_integrand = grs_avg[:,np.newaxis]*sh_shell
        rr_ord = np.repeat(self.get_r_rad()[:,np.newaxis],self.ncombos,axis=1)
        if self.dim==3: s_o_integrand*=(np.power(rr_ord,2.0)*4.*np.pi)
        else:           s_o_integrand*=(rr_ord*2.*np.pi)
        ent_o = simps(s_o_integrand,rr_ord,axis = 0)
        print(st); print(ent_o, -0.5 * ent_o * KB_CAL_MOL * NUM_DENSITY_H2O,
              -0.5*sum(ent_o)*KB_CAL_MOL*NUM_DENSITY_H2O)
        return -0.5 * ent_o * KB_CAL_MOL * NUM_DENSITY_H2O

    def ang_cmbs(self):
        '''Given desired order of entropy calc, generate array of all combos'''
        self.angle_combos = []
        for subset in itertools.combinations(range(NANGLES), self.order):
            self.angle_combos.append(list(subset))
        self.ncombos = len(self.angle_combos)

    def integ_fact(self):
        ''' Given a list of angles and an array shape, will create an array of
            shape = shape, to scale data by. Assumes that d1 = distance,
            and the other dimensions are angles 
        '''
        self.ifacts = np.ones(tuple([self.ncombos])+self.hist_shape); 
        self.nfacts = np.ones(self.ncombos)
        for c in range(self.ncombos): 
            for i in range(self.order): # Xfrm theta,1st col=dist
                self.nfacts[c] *= ANG_NORM[self.angle_combos[c][i]]
                if self.angle_combos[c][i] > 2: 
                    sft = np.sin(self.radi_a)
                    for j in range(self.order):
                        if i > j:   sft = np.expand_dims(sft, axis = 0)
                        elif i < j: sft = np.expand_dims(sft, axis = -1)
                    sft = np.expand_dims(sft, axis = 0)
                    self.ifacts[c] *= sft

    def get_rbn(self):
        '''Return bin edges of g(r)'''
        return self.grC.bns_r

    def get_rbn_ct(self):
        '''Return bin count of g(r)'''
        return self.grC.nb_r

    def get_r_rad(self):
        '''Return centers of g(r) bins'''
        return self.grC.radi_r

    def plot_g_ang(self, g_mean):
        '''Method for plotting the g(angle) distributions for 1&2nd order'''
        matplotlib.rcParams.update({'font.size': 4})
        if self.order == 1: 
            f,axes=plt.subplots(5,1,sharex='col',sharey='row',figsize=(1.3,3.5))
        else:
            f,axes=plt.subplots(5,2,sharex='col',sharey='row',figsize=(2,5)) 
            print(axes.shape, g_mean.shape)
            plt.setp(axes.flat, aspect=1.0, adjustable='box-forced')
        axes = axes.flatten()
        bin_rng = [ 27, 30, 43, 65]

        for bn in range(self.get_rbn_ct()): #for each r bin, plot g(angle)
            if (bn in bin_rng and self.order == 1):
                for j in range(self.ncombos):
                    axes[j].plot(self.radi_a,g_mean[j][bn],
                                 label=str(self.get_r_rad()[bn]))
            if (bn == 27 and self.order == 2):
                X, Y = np.meshgrid(self.radi_a, self.radi_a); c = []
                for j in range(self.ncombos):
                    nm = self.angle_combos[j]
                    c.append(axes[j].contour(X, Y, g_mean[j][bn].T))
                    axes[j].text(.5,1.04,"g("+ANG_NM[nm[0]]+", "+ANG_NM[nm[1]]+")",
                                 horizontalalignment='center',fontsize = 5,
                                 transform=axes[j].transAxes)
                    plt.clabel(c[j], inline=1, fontsize=2);
        max_tc = 4
        if self.order == 1:
            axes[0].legend(ncol=len(bin_rng),columnspacing=-0.1,labelspacing=-1.95,
                           borderaxespad=-1.9,handlelength=1.8,fontsize=4)
            axes[-1].set_xlabel("Angle (radians)",fontsize=6)
            axes[-1].xaxis.labelpad = -1
            for j in range(self.ncombos):
                axes[j].xaxis.set_major_locator(plt.MaxNLocator(max_tc))
                axes[j].yaxis.set_major_locator(plt.MaxNLocator(max_tc))
                axes[j].tick_params(axis='both', which='major', pad= 1.3)
                axes[j].set_xlim(0, np.pi);axes[j].yaxis.labelpad = -1
                axes[j].set_ylim(0, 3);
                axes[j].set_ylabel("g("+ANG_NM[self.angle_combos[j][0]]+")",
                                  fontsize=6)
                axes[j].text(.7,.6,ANG_NM[self.angle_combos[j][0]], 
                             horizontalalignment='center',
                             transform=axes[j].transAxes,fontsize=10);
        plt.subplots_adjust(wspace=0.14, hspace=0.14)
        plt.savefig("gangle_"+str(self.order)+'.png',format='png',
                    bbox_inches='tight',dpi=300)

    def grid_for_print(self):
        '''Grid for a varying dimensional g(angle)'''
        ba,ra,sa = self.nb_a,self.radi_a,self.binsiz_a # No. angles for hist
        print(ba, ra, sa)
        if self.order == 1:
            inds = np.mgrid[0:ba:1][:,np.newaxis]; xy = ra[:,np.newaxis]
        if self.order == 2:
            inds = np.mgrid[0:ba:1,0:ba:1].reshape(2,-1).T
            xy=np.mgrid[ra[0]:ra[-1]+sa:sa,ra[0]:ra[-1]+sa:sa].reshape(2,-1).T
        if self.order == 3:
            inds=np.mgrid[0:ba:1,0:ba:1,0:ba:1].reshape(3,-1).T
            xy = np.mgrid[ra[0]:ra[-1]+sa:sa,ra[0]:ra[-1]+sa:sa,
                          ra[0]:ra[-1]+sa:sa].reshape(3,-1).T
        if self.order == 4:
            inds=np.mgrid[0:ba:1,0:ba:1,0:ba:1,0:ba:1].reshape(4,-1).T
            xy = np.mgrid[ra[0]:ra[-1]+sa:sa,ra[0]:ra[-1]+sa:sa,
                          ra[0]:ra[-1]+sa:sa,ra[0]:ra[-1]+sa:sa].reshape(4,-1).T
        if self.order == 5:
            inds=np.mgrid[0:ba:1,0:ba:1,0:ba:1,0:ba:1,0:ba:1].reshape(5,-1).T
            xy = np.mgrid[ra[0]:ra[-1]+sa:sa,ra[0]:ra[-1]+sa:sa,
                          ra[0]:ra[-1]+sa:sa,ra[0]:ra[-1]+sa:sa,
                          ra[0]:ra[-1]+sa:sa].reshape(5,-1).T
        return inds, xy

    def print_gang(self):
        '''Print histogram for each group for each timestep for multi
           dimensional g(\omega) '''
        inds, xy = self.grid_for_print(); radi_ra = self.get_r_rad()
        his_ct = open("angle_g_bin_ct.csv","w"); tt = self.ntot[0].flatten()
        his_ct.write("bin,ct\n")
        for i in range(len(tt)): 
            his_ct.write("{0:.3f},{1}\n".format(radi_ra[i],tt[i]))
        his_ct.close()

        for bn in range(len(radi_ra)): #for each r bin and angle combo
            for an in range(self.ncombos):
                ans, st = "", ""
                for i in self.angle_combos[an]: 
                    ans += str(i)+"_"
                    st += "bin"+str(i)+","
                f = open("angle_g_o{0}_bn{1:.4f}.csv".
                          format(ans[:-1],radi_ra[bn]), "w")
                for t in range(self.ntimes): st += ("time"+str(t)+",")
                f.write(st+"timetot\n")

                for ab in range(inds.shape[0]):
                    st = ""
                    for di in range(inds.shape[1]): 
                        st+="{0:.4f},".format(xy[ab][di])
                    for t in range(self.ntimes+1): 
                        dat = self.gs[t][an][bn]
                        st += "{0:.6f},".format(dat[tuple(inds[ab])])
                    f.write(st[:-1]+"\n")
                f.close()


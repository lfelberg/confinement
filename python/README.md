# Python analysis scripts for confinement project

This file contains an overview of files used for the confinement project.

1.  compute_mean.py
3.  csvfile.py
3.  dics.py
3.  distance_graphene.py
3.  energy_lmp.py
4.  entropy_calc.py
0.  entropy_calc_util.py
3.  fourier.py
5.  g_of_r.py
4.  gauss_fits.py
3.  msd_water.py
2.  order_params_2D.py
3.  order_params_3D.py
3.  util.py
2.  volfile.py
3.  water_3B_angles.py
3.  water_angles.py
3.  water_angles_util.py
3.  water_deviation.py
4.  write_pdb.py
5.  xdensity.py
0.  xyzfile.py


For plotting the following scripts are in the `plotting` directory:

3.  plot_3B_angles.py
3.  plot_3B_angles_multi.py
3.  plot_dist_gg.py
4.  plot_gr.py
5.  plot_msd.py
0.  plot_order.py
3.  plot_vs_2D_density.py
5.  plot_vs_time.py

## compute_mean.py

Given a csv file, compute the mean for the specified column

**usage**: `python compute_mean.py csvfname col#`

Opens `csvfname` and computes the mean and stdev for column `col#`.


## csvfile.py

A class for reading in CSV files.

**usage**: `cC = CSVFile(csvname)`

In the object, saves the filename as `csvfname` and the data as an array, `dat`
and the keys to the data as `key`.


## dics.py

Contains two dictionaries:

1.  `colorL` is a dictionary of colors for plotting a rainbow of colors.
2.  `dens` is a dictionary of the different 2D densities in the system, stored by
their initial separation. The data is an array with entries (0) = name of the density,
(1) = the color to plot it with, (2) = the line type (filled, dashed) to plot with.


## distance_graphene.py

A program to take an XYZ file and for each snapshot, compute the graphene-graphene
separation as well as the distance of the water molecules between the sheets to a
given wall.

**usage**: `python distance_graphene.py *.xyz sep len iter`
Will open the file `*.xyz` and also expects a vol file named `run{sep}_{len}_{iter}.vol`.

**Outputs**: Prints out a file: `run{sep}_{len}_{iter}.distgg`, which is a csv file
with columns: 

1.  dg1, the distance between the given water molecule and graphene sheet 1
1.  dg2, the distance between the given water molecule and graphene sheet 2
2.  dgg, the distance between the two graphene patches that enclose the water 

Note: There is a commented section to be used with the smallest flexible system
to count regions with no water, where the graphene walls are in contact.
Also, for the rigid system, there is a lot of displacement in the x direction, 
and this program is not able to catch all instances of when the walls are wrapped
by the box, so it just continues to next if it catches some unusual results.


## energy_lmp.py

File to get info from lammps output file, specifically the `thermo` section.

**usage**: `python energy_lmp.py *.out` or `python energy_lmp.py *.ext`

**Outputs**: Will read either a lammps output file (with the custom output dump
displayed in the beginning of this python file) or a simple `*.vol` style file and
print out the energies to a csvfile and compute right now the compressibility!


## entropy_calc.py

Method to compute the entropy from oxygen pair interactions (g(R) and angle distributions)
from the paper: **Orientational correlations and entropy in liquid water** by 
Lazaridis and Karplus (1996).

From the output of a file that computes the distances and 5 angles
between all pairs of atoms in the system, compute the translational or orientational
entropy.

**usage**: `python entropy_calc.py angfname sep len iter entropy_type [orien_order]`
Entropy type is either: `trans`, `orien` or `both`. If `orien` or `both` is chosen,
you also need to specify the order of then orientational calculation (1-5).

The file `angfname` can be a csv file of all the pair distances and angles, but it
can also be a g(R) or g(angle) file, that allows the user to use these intermediates
to calculate the entropy quicker.

**Outputs**: If an angle csv file is input, the g(R or angle) is computed and printed
out. If desired, this is correlation function is also plotted. Then the translational
or orientational entropy is calculated from that correlation function as described in
the paper above. 

Notes: the translational entropy is calculated with a smaller dR bin than the 
orientational, and therefore the `both` command may not work. Also, it is important that
if using the g(\*) option that the dR used is the same in the main here as it is 
in the inputs to these methods.


## entropy_calc_util.py

This is the main body of the entropy calculations. There is a class for the translational
and the orientational entropy functions, and both have the following:

1.  Plotting method: to plot the g(\*) computed. For the orientational, you can only 
plot 1 and 2nd order curves.
2.  Printing method: print out g(\*). For orientational, one printed out for each angle combo
and each graphene separation. 
3.  Integration method: will integrate according to paper. 


## fourier.py

Compute the fourier transform of the g(R) to create an S(q) function. Currently
plots but will probably also save as csv file

**usage**: `python fourier.py *.gr sep len iter`

**Outputs**: Plot of fourier transform.


## g_of_r.py

A program to calculate the g(R) for a given pair of atom types. 

**usage**: `python g_of_r.py *xyz sep len iter #Pairs Pr0_0 Pr0_1 Pr1_0 Pr1_1 ... `

**Outputs**:
1. g_r_2D_sep_len_iter_PrN0_PrN1.csv - a csv file with histogram as a function
of x distance from graphene wall of 2D g(R) of pair.
1. g_r_3D_sep_len_iter_PrN0_PrN1.csv - a csv file with histogram as a function
of x distance from graphene wall of 3D g(R) of pair.

Note: currently plot of 3D g(R) is disabled but it should be pretty easy to 
implement.


## gauss_fits.py

File that contains a variety of gaussian functions: 
1.  Normal distribution
2.  Skewed normal
3.  Multiple gaussian curves

That can be fit to the data.


## msd_water.py

Computing the MSD of water molecules, using the displacement of water oxygens.
The program unwraps the atoms, so all atoms should progress out from original
center of mass. Also, the center of mass of each group of waters (between and outside
of walls) is subtracted to remove any collective motion effects. 

**usage**: `python msd_water.py *.xyz sep len iter`

**Outputs**: `msd_sep_len_iter.csv`, a file that contains MSD calcs and their 
standard deviations for X, Y, Z, YX and XYZ dimensions. 


## order_params_2D.py

Plotting 2D order parameters for bond order: phi_n = < sum_{Nb} exp(i*n*theta) >.
This is taken from the methodology of **Dislocation-mediated melting in two dimensions**
by David R. Nelson and B.I. Halperin (1979). There is a good discussion of the method
in **Formation of Two-Dimensional Crystals with Square Lattice
Structure from the Liquid State** by Vo Van Hoang and Nguyen Thanh Hieu (2016). 
 
**usage**: `python order_params_2D.py *.xyz sep len iter cutoff n`
Where `cutoff` is the distance between oxygen pairs to consider them neighbors.
N is the number in the exponential.

**Outputs**: Csv file of phi for each water molecule.


## order_params_3D.py

**This program may not be correct as is**.

Method to calculate the 3D order parameters q_4 bar and q_6 bar. This is for
quantifying the structure of ice in 3D. From the paper: **Bond-orientational order in liquids and glasses**
by Paul J. Steinhardt, David R. Nelson, and Marco Ronchetti (1983). A good
summary of it is also in **Accurate determination of crystal structures based on averaged local bond order
parameters** by Wolfgang Lechner and Christoph Dellago (2008). 

**usage**: `python order_params_3D.py *.xyz sep len iter`

**Outputs**: Csv file of q4 and q6 bar.


## plot_3B_angles.py 

Plotting the distribution of 3-body angle.

**usage**: `python plot_3B_angles.py *.csv sep len iter`

**Outputs**: Plot of 3B angle distribution, can also plot fits to this curve 
and distribution of distances between water nearest neighbors. Will also save
the distribution of angles as a csv file.


## plot_3B_angles_multi.py

**usage**: `python plot_3B_angles_multi.py csvStart nsep nlen iter sep1 sep2... len1 len2... ext`
For a given list of confinement systems, plot their 3-body angle distributions 
all on one plot.

**Outputs**: A plot of the angle distributions from 0 -> 180 for all systems input.


## plot_dist_gg.py 

Plot a 2D histogram of the distance of graphene-graphene separation and the 
water distances between them. Also plot a 1D histogram of dgg and xdist. X dist 
is normalized with respect to dgg, so it will range from 0 - 1. Realistically,
because of VDW, the true range ~0.25-0.75. This also sorts out any unrealistic 
values, which are not all caught in the distance_graphene program. This can
also calculate fits for d_gg and x.

**usage**: `python plot_dist_gg.py *.csv sep len iter`


## plot_gr.py

This will plot a series of curves for g(r). 

**usage**: `python plot_gr.py csvStart nsep nlen iter sep1 sep2... len1 len2... ext`

Should output a plot of the multiple datasets.


## plot_msd.py

Will plot multiple curves for MSD, and will also approximate D_2D for many
files. If you want to change from D_2D, you need to change the scaling factor.

**usage**: `python plot_msd.py csvStart nsep nlen niter sep1 sep2... len1 len2... iter1 iter2... ext datLoc`


## plot_order.py

Plot the distribution of order parameters q or phi.

**usage**: `python plot_order.py *.csv sep len iter`


## plot_vs_2D_density.py

Plot the values of something for rigid and flexible confinement versus 2D
density. This is currently either the d_gg, the compressibility or the estimated
$D_{||}$.

**usage**: `python plot_vs_2D_density.py *.csv`


## plot_vs_time.py

Will plot multiple curves for an input csv file type. Assumes that the first
column of the csv is some sort of sequential increase.

**usage**: `python plot_vs_time.py csvStart nsep nlen niter sep1 sep2... len1 len2... iter1 iter2... ext datLoc`


##  util.py

A file that contains some common methods for many of the python scripts.

## water_3B_angles.py

**usage**: `python water_3B_angles.py *.xyz sep len iter`

**Outputs**:  A csv file with the distance between an oxygen and its nearest neighbors
and the angle between them.

## water_angles.py

Program to compute pairwise angles and distance for waters in the system. This
output is then used to calculate the entropy with the g(R) and the g(angles).

**usage**: `python water_angles.py *.xyz sep len iter`

**Outputs**: A csv file that has the distance and 5 descriptive angles
for pairs of atoms in the system.


## water_angles_util.py

This file contains the methods for calculating pairwise distances and angles
between waters in the system. Has equations for angle calculation, calculating 
distances between pairs, finding n closest pairs. Some methods here are used for
other programs as well.


## volfile.py

A class for making/printing/reading volume-type files.

**usage**: `python volfile.py lmp_run.out`

Will generate a lmp_run.vol file that has every 5th snapshot time, step #,
xmin, xmax, ymin, ymax, zmin, zmax on each line.

**usage** as class: `vC = VolFile(volname)`

**Outputs**: Will create an instance of a VolFile class where vC reads in a \*.out or \*.vol
file and can operate on it.

## write_pdb.py

Given an XYZ file, convert to pdb formatted trajectory. Each snapshot delimited
by an END line.

**usage**: `python write_pdb.py *.xyz sep len iter`
Also expects a \*.vol file.

**Outputs**: will create an output \*.pdb file with same start as \*.xyz

## get_density.py

Will calculate the density of a system or plot it, depending on the input files.

If invoked with file \*.xyz extension, will print out CSV files with a histogram
for the length of the simulation box in X and one column for each atom type.

**usage**: `python get_density.py *.xyz sep len iter`
Also expects a \*.vol file

**Outputs**: will print out a CSV file with the length of the simulation box in X 
and one column for each atom type.

Else, will assume that the input file is a csv with histograms and will print plots.

**usage**: `python get_density.py *.dens_hist sep len iter`
Also expects a \*.vol file

**Outputs**: will print out one figure for each type: `dens_#.png`.

## xyzfile.py

A class for the xyzfiles. When using constructor, will save the filename as
a member `xyzfname`. Contains several methods for getting lists of different atom 
types etc.

**usage**: `xC = XYZFile(fname, VolFile)`

Where `fname` is the name of the xyz-style file and VolFile is a class of 
volume style.



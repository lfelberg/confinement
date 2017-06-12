# Python analysis scripts for confinement project

This file contains an overview of files used for the confinement project.

2.  compute_mean.py
3.  csvfile.py
3.  dics.py
3.  distance_graphene.py
3.  energy_lmp.py
4.  entropy_calc.py
0.  entropy_calc_util.py
3.  fourier.py
5.  g_of_r.py
4.  gauss_fits.py
3.  get_density.py
2.  get_distances.py
3.  hbond_stats.py
3.  msd_water.py
3.  ntot_make.py
3.  order_params_2D.py
3.  order_params_3D.py
4.  plot_2d_heat.py
5.  plot_3B_angles.py
0.  plot_3B_angles_multi.py
3.  plot_diffusion.py
5.  plot_dist_gg.py
4.  plot_dist_gg_all.py
3.  plot_gr.py
2.  plot_msd.py
3.  plot_order.py
3.  plot_vs_time.py
3.  rmsd.py
3.  util.py
3.  volfile.py
4.  water_3B_angles.py
5.  water_angles.py
0.  water_angles2D.py
3.  water_angles_util.py
5.  water_deviation.py
4.  write_pdb.py
3.  xyzfile.py


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

**usage**: `python distance_graphene.py \*.xyz sep len iter`
Will open the file `\*.xyz` and also expects a vol file named `run{sep}_{len}_{iter}.vol`.

**Outputs**: Prints out a file: `run{sep}_{len}_{iter}.distgg`, which is a csv file
with columns: 

1.  dg1, the distance between the given water molecule and graphene sheet 1
1.  dg2, the distance between the given water molecule and graphene sheet 2
2.  dgg, the distance between the two graphene patches that enclose the water 

Note: There is a commented section to be used with the smallest flexible system
to count regions with no water, where the graphene walls are in contact.

## energy_lmp.py

**usage**: `python \*.out` or `python \*.ext`

**Outputs**: Will read either a lammps output file (with the custom output dump
displayed in the beginning of this python file) or a simple `\*.vol` style file and
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

3.  fourier.py
5.  g_of_r.py
4.  gauss_fits.py
3.  get_density.py
2.  get_distances.py
3.  hbond_stats.py
3.  msd_water.py
3.  ntot_make.py
3.  order_params_2D.py
3.  order_params_3D.py
4.  plot_2d_heat.py
5.  plot_3B_angles.py
0.  plot_3B_angles_multi.py
3.  plot_diffusion.py
5.  plot_dist_gg.py
4.  plot_dist_gg_all.py
3.  plot_gr.py
2.  plot_msd.py
3.  plot_order.py
3.  plot_vs_time.py
3.  rmsd.py
3.  util.py
4.  water_3B_angles.py
5.  water_angles.py
0.  water_angles2D.py
3.  water_angles_util.py
5.  water_deviation.py

## fourier.py

**usage**: `python \*.xyz sep len iter`

**Outputs**: 

## 

**usage**: `python \*.xyz sep len iter`

**Outputs**: 

## 

**usage**: `python \*.xyz sep len iter`

**Outputs**: 

## 

**usage**: `python \*.xyz sep len iter`

**Outputs**: 

## 

**usage**: `python \*.xyz sep len iter`

**Outputs**: 

## 

**usage**: `python \*.xyz sep len iter`

**Outputs**: 

## 

**usage**: `python \*.xyz sep len iter`

**Outputs**: 

## 

**usage**: `python \*.xyz sep len iter`

**Outputs**: 

## 

**usage**: `python \*.xyz sep len iter`

**Outputs**: 


## get_density.py

Will calculate the density of a system or plot it, depending on the input files.

If invoked with file \*.xyz extension, will print out CSV files with a histogram
for the length of the simulation box in X and one column for each atom type.

**usage**: `python get_density.py \*.xyz sep len iter`
Also expects a \*.vol file

**Outputs**: will print out a CSV file with the length of the simulation box in X 
and one column for each atom type.

Else, will assume that the input file is a csv with histograms and will print plots.

**usage**: `python get_density.py *.dens_hist sep len iter`
Also expects a \*.vol file

**Outputs**: will print out one figure for each type: `dens_#.png`.

## get_distances.py

Will calculate the distance of non-graphene atoms to the graphene wall. Will
also calculate the distance between the graphene sheets for each frame and 
reports the closest graphene atom number.

**usage**: `python get_distances.py *.xyz sep len iter`
Also expects a \*.vol file.

**Outputs**: 

1. runS_L_it.dist - a file with XYZ structure, reports the distance of all atoms
to the closest graphene atom. Will store the # of graphene atom and the distance.
For graphene atoms, reports for the one wall the distance to the other. For the other 
wall, has 0.0 as values.
2. runS_L_it_graph[0-1].dat - a file with XYZ structure, reports the deviation
of each graphene in wall 0 or 1 from the average X of the wall
3. runS_L_it.grap_sep0.dat - a file with XYZ strucute, reports the separation
of each atom in wall 0 to the closest point in wall 1 without PBCS.

## g_of_r.py

A program to calculate the g(r) for a given pair of atom types. 

**usage**: `python g_of_r.py *xyz sep len iter #Pairs Pr0_0 Pr0_1 Pr1_0 Pr1_1 ... `

**Outputs**:
1. g_r_3D_sep_len_iter_PrN0_PrN1.csv - a csv file with histogram as a function
of x distance from graphene wall of 3D g(r) of pair.
2. FUTURE: 2D g(r)

## plot_2d_heat.py

This will plot a 2D heatmap of either: 

1. The distance between the two graphene plates (filename: runS_L_It_graph[0-1].data)
2. The displacement of each carbon atom from the average x location of the plate
(filename: runS_L_It_graph_sep0.data)

**usage**: `python plot_2d_heat.py runS_L_It_graph*.dat` 

**Outputs**: Will produce a series of plots, one for each snapshot in the dat file.

## plot_gr.py

This will plot a series of curves for g(r). 

**usage**: `python plot_gr.py csvStart nsep nlen niter sep1 sep2... len1 len2... iter1 iter2... ext datLoc`

Should output a plot of the multiple datasets.

## plot_vs_x.py

Plot multiple csv files of data for a versus X coord plot.

**usage**: `python plot_vs_x.py csvStart nsep nlen niter sep1 sep2... len1 len2... iter1 iter2... ext datLoc`

Should output a plot of the multiple datasets.


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

## xyzfile.py

A class for the xyzfiles. When using constructor, will save the filename as
a member `xyzfname`. Contains several methods for getting lists of different atom 
types etc.

**usage**: `xC = XYZFile(fname, VolFile)`

Where `fname` is the name of the xyz-style file and VolFile is a class of 
volume style.



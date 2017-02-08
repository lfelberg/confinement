# Python analysis scripts for confinement project

This file contains an overview of files used for the confinement project.

## csvfile.py

A class for reading in CSV files.

usage: `cC = CSVFile(csvname)`

In the object, saves the filename as `csvfname` and the data as an array, `dat`.


## get_density.py

Will calculate the density of a system or plot it, depending on the input files.

If invoked with file \*.xyz extension, will print out CSV files with a histogram
for the length of the simulation box in X and one column for each atom type.

usage: `python get_density.py \*.xyz sep len iter`
Also expects a \*.vol file

Outputs: will print out a CSV file with the length of the simulation box in X 
and one column for each atom type.

Else, will assume that the input file is a csv with histograms and will print plots.

usage: `python get_density.py \*.dens_hist sep len iter`
Also expects a \*.vol file

Outputs: will print out one figure for each type: `dens_#.png`.

## get_distances.py

Will calculate the distance of non-graphene atoms to the graphene wall. Will
also calculate the distance between the graphene sheets for each frame and 
reports the closest graphene atom number.

usage: `python get_distances.py \*.xyz sep len iter`
Also expects a \*.vol file.

Outputs: 
1. runS_L_it.dist - a file with XYZ structure, reports the distance of all atoms
to the closest graphene atom. Will store the # of graphene atom and the distance.
For graphene atoms, reports for the one wall the distance to the other. For the other 
wall, has 0.0 as values.
2. runS_L_it_graph[0-1].dat - a file with XYZ structure, reports the deviation
of each graphene in wall 0 or 1 from the average X of the wall
3. runS_L_it.grap_sep0.dat - a file with XYZ strucute, reports the separation
of each atom in wall 0 to the closest point in wall 1 without PBCS.

## g_of_r.py


## hbonanza.py

This will calculate the hydrogen bonds in a pdb file. See the file for usage,
because the invocation is long.

Outputs: \*.average_hbonds and \*.frame_by_frame_hbonds.csv and \*.tcl


## hbond_stats.py

This script tries to use the file output by hbonanza.py, \*.frame_by_frame_hbonds.csv
and syncs up that with the distance from the wall in a \*.dist file to calculate the 
number of hydrogen bonds vs distance from the wall.

usage: `python hbond_stats.py xyzfname sep len iter`
Also expects to have a \*.vol file and a \*.dist file in the directory

Will also assume that there is a file with same start as \*.xyz that has the
extension: \*.frame_by_frame_hbonds.csv.

Outputs: 
1. \*.hbond_dat - file with bins from 0 - half carbon sep with hydrogen bonds per
water molecule.

## plot_2d_heat.py

This will plot a 2D heatmap of either: 
1. The distance between the two graphene plates (filename: runS_L_It_graph[0-1].data)
2. The displacement of each carbon atom from the average x location of the plate
(filename: runS_L_It_graph_sep0.data)

usage: `python plot_2d_heat.py runS_L_It_graph*.dat` 

Outputs: Will produce a series of plots, one for each snapshot in the dat file.

## plot_vs_x.py


## volfile.py

A class for making/printing/reading volume-type files.

usage: `python volfile.py lmp_run.out`

Will generate a lmp_run.vol file that has every 5th snapshot time, step #,
xmin, xmax, ymin, ymax, zmin, zmax on each line.

usage as class: `vC = VolFile(volname)`

Outputs: Will create an instance of a VolFile class where vC reads in a \*.out or \*.vol
file and can operate on it.

## write_pdb.py

Given an XYZ file, convert to pdb formatted trajectory. Each snapshot delimited
by an END line.

usage: `python write_pdb.py \*.xyz sep len iter`
Also expects a \*.vol file.

Outputs: will create an output \*.pdb file with same start as \*.xyz

## xyzfile.py

A class for the xyzfiles. When using constructor, will save the filename as
a member `xyzfname`.

usage: `xC = XYZFile(fname, VolFile)`

Where `fname` is the name of the xyz-style file and VolFile is a class of 
volume style.



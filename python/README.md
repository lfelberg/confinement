# Python analysis scripts for confinement project

This file contains an overview of files used for the confinement project.

## csvfile.py

A class for reading in CSV files.

usage: `cC = CSVFile(csvname)`

Saves the filename as `csvfname` and the data as an array, `dat`.


## get_density.py


## get_distances.py


## g_of_r.py


## hbonanza.py


## hbond_stats.py


## plot_2d_heat.py

## plot_vs_x.py


## volfile.py

A class for making/printing/reading volume-type files.

usage: `python volfile.py lmp_run.out`

Will generate a lmp_run.vol file that has every 5th snapshot time, step #,
xmin, xmax, ymin, ymax, zmin, zmax on each line.

usage as class: `vC = VolFile(volname)`

Will create an instance of a VolFile class where vC reads in a \*.out or \*.vol
file and can operate on it.

## write_pdb.py


## xyzfile.py

A class for the xyzfiles. When using constructor, will save the filename as
a member `xyzfname`.

usage: `xC = XYZFile(fname, VolFile)`

Where `fname` is the name of the xyz-style file and VolFile is a class of 
volume style.



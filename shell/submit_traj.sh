#!/bin/bash

# script to move files for each replica of traj

for sep in 16 ; do #16 20 ; do
    cd ${sep}
    for iter in 1 2 3 4 5 6 7 8 9 10; do
        cp ../sub_lammps.qsub sub_lmp.${iter}.qsub
        cp ../water_graph.inp run${sep}_${iter}.inp
        sed -i "s/sep/${sep}/g" sub_lmp.${iter}.qsub run${sep}_${iter}.inp
        sed -i "s/iter/${iter}/g" sub_lmp.${iter}.qsub run${sep}_${iter}.inp
        sed -i "s/seed/${iter}/g" run${sep}_${iter}.inp

        qsub -q all.q sub_lmp.${iter}.qsub
    done
    cd ../
done

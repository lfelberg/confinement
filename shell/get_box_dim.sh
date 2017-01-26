#!/bin/bash

# script to move files for each replica of traj

for sep in 20 ; do #16 20 ; do
    cd ${sep}
    for iter in 1 2 3 4 5 ; do
        python ../../analysis/get_box_dim.py run${sep}_${iter}.out
    done
    cd ../
done

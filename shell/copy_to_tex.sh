#!/bin/bash
# Script to copy images to latex file: April 21, 2017

cp msd/diffusion_coeff_2D.png ~/confinement/analysis/paper/flexible/diffusion_coeff_2D_f_r.png
cp d_gg.png ~/confinement/analysis/paper/flexible/
cp compressibility.png ~/confinement/analysis/paper/flexible/

 cd g_r/
 python3 ~/confinement/analysis/python/plot_gr.py g_r_2D_ 4 1 1 37 6 7 8 46 _1_1.csv
 python3 ~/confinement/analysis/python/plot_gr.py g_r_2D_ 5 1 1 37 91 92 10 11 46 _1_1.csv
 python3 ~/confinement/analysis/python/plot_gr.py g_r_2D_ 6 1 1 37 12 141 142 161 162 46 _1_1.csv
 cd ../rigid/g_r
 python3 ~/confinement/analysis/python/plot_gr.py g_r_2D_ 4 1 1 37 6 7 8 46 _1_1.csv
 python3 ~/confinement/analysis/python/plot_gr.py g_r_2D_ 4 1 1 37 9 10 11 46 _1_1.csv
 python3 ~/confinement/analysis/python/plot_gr.py g_r_2D_ 6 1 1 37 12 141 142 161 162 46 _1_1.csv
 cd ../../ 

 cd angles/
 python3 ~/confinement/analysis/python/plot_3B_angles_multi.py run 4 1 1 37 6 7 8 46 _3B_angles_hist.csv
 python3 ~/confinement/analysis/python/plot_3B_angles_multi.py run 4 1 1 37 9 10 11 46 _3B_angles_hist.csv
 python3 ~/confinement/analysis/python/plot_3B_angles_multi.py run 4 1 1 37 12 14 16 46 _3B_angles_hist.csv
 cd ../rigid/angles
 python3 ~/confinement/analysis/python/plot_3B_angles_multi.py run 4 1 1 37 6 7 8 46 _3B_angles_hist.csv
 python3 ~/confinement/analysis/python/plot_3B_angles_multi.py run 4 1 1 37 9 10 11 46 _3B_angles_hist.csv
 python3 ~/confinement/analysis/python/plot_3B_angles_multi.py run 4 1 1 37 12 14 16 46 _3B_angles_hist.csv
 cd ../../ 

cp g_r/g_r_2D_37_6_7_8_1_1.png ~/confinement/analysis/paper/flexible/g_r_2D_6_7_8_f.png
cp g_r/g_r_2D_37_92_91_10_11_1_1.png ~/confinement/analysis/paper/flexible/g_r_2D_9_10_11_f.png
cp g_r/g_r_2D_37_12_14_16_1_1.png ~/confinement/analysis/paper/flexible/g_r_2D_12_14_16_f.png

cp rigid/g_r/g_r_2D_37_6_7_8_1_1.png ~/confinement/analysis/paper/rigid/g_r_2D_6_7_8_r.png
cp rigid/g_r/g_r_2D_37_9_10_11_1_1.png ~/confinement/analysis/paper/rigid/g_r_2D_9_10_11_r.png
cp rigid/g_r/g_r_2D_37_12_14_16_1_1.png ~/confinement/analysis/paper/rigid/g_r_2D_12_14_16_r.png

cp angles/ang_3B_37_6_7_8.png ~/confinement/analysis/paper/flexible/ang_3B_6_7_8_f.png
cp angles/ang_3B_37_9_10_11.png ~/confinement/analysis/paper/flexible/ang_3B_9_10_11_f.png
cp angles/ang_3B_37_12_14_16.png ~/confinement/analysis/paper/flexible/ang_3B_12_14_16_f.png

cp rigid/angles/ang_3B_37_6_7_8.png ~/confinement/analysis/paper/rigid/ang_3B_6_7_8_r.png
cp rigid/angles/ang_3B_37_9_10_11.png ~/confinement/analysis/paper/rigid/ang_3B_9_10_11_r.png
cp rigid/angles/ang_3B_37_12_14_16.png ~/confinement/analysis/paper/rigid/ang_3B_12_14_16_r.png


for i in 6 7 8 9 10 11 12 13 14 16 ; do 
#for i in 16 ; do 
    for j in 46 ; do
        cd ${i}_${j}_1/
        echo $i $j "flexible!"
       #if [ $i -gt 8 ] ; then
       #    python3 ~/confinement/analysis/python/plot_3B_angles.py run${i}_${j}_1_3B_angles_layers.csv ${i} ${j} 1
       #else
       #    python3 ~/confinement/analysis/python/plot_3B_angles.py run${i}_${j}_1_3B_angles.csv ${i} ${j} 1
       #fi
       #python3 ~/confinement/analysis/python/plot_dist_gg.py run${i}_${j}_1.distgg ${i} ${j} 1
       #python3 ~/confinement/analysis/python/plot_gr.py g_r_2D_ 1 1 1 ${i} ${j} _1_1.csv 2
       #python3 ~/confinement/analysis/python/plot_gr.py g_r_2D_ 1 1 1 ${i} ${j} _1_2.csv 2
       #python3 ~/confinement/analysis/python/plot_rmsd.py msd_ 1 1 1 ${i} ${j} 1 .csv 1
       #python3 ~/confinement/analysis/python/plot_rmsd.py msd_ 1 1 1 ${i} ${j} 1 .csv 7
        cd ../

        cd rigid/${i}_${j}_1/
        echo $i $j "rigid!"
       #python3 ~/confinement/analysis/python/plot_3B_angles.py run${i}_${j}_1_3B_angles_layers.csv ${i} ${j} 1
       #python3 ~/confinement/analysis/python/plot_dist_gg.py run${i}_${j}_1.distgg ${i} ${j} r
       #python3 ~/confinement/analysis/python/plot_gr.py g_r_2D_ 1 1 1 ${i} ${j} _1_1.csv 2
       #python3 ~/confinement/analysis/python/plot_gr.py g_r_2D_ 1 1 1 ${i} ${j} _1_2.csv 2
       #python3 ~/confinement/analysis/python/plot_rmsd.py msd_ 1 1 1 ${i} ${j} 1 .csv 1
       #python3 ~/confinement/analysis/python/plot_rmsd.py msd_ 1 1 1 ${i} ${j} 1 .csv 7
        cd ../../

        # oxygen-graphene and gg dist plots
       #cp rigid/${i}_${j}_1/run${i}_${j}_1.disg_sep_fit.png ~/confinement/analysis/paper/rigid/dgg_fit_${i}_${j}_r.png
       #cp ${i}_${j}_1/run${i}_${j}_1.disg_sep_fit.png ~/confinement/analysis/paper/flexible/dgg_fit_${i}_${j}_f.png

       #cp rigid/${i}_${j}_1/run${i}_${j}_1.dis.png ~/confinement/analysis/paper/rigid/dgg_v_x_${i}_${j}_r.png
       #cp ${i}_${j}_1/run${i}_${j}_1.dis.png ~/confinement/analysis/paper/flexible/dgg_v_x_${i}_${j}_f.png

        # gr plots
       #cp rigid/${i}_${j}_1/g_r_2D_1_1.png ~/confinement/analysis/paper/rigid/g_r_${i}_${j}_r_o_o.png
       #cp rigid/${i}_${j}_1/g_r_2D_1_2.png ~/confinement/analysis/paper/rigid/g_r_${i}_${j}_r_o_h.png
       #cp ${i}_${j}_1/g_r_2D_1_1.png ~/confinement/analysis/paper/flexible/g_r_${i}_${j}_f_o_o.png
       #cp ${i}_${j}_1/g_r_2D_1_2.png ~/confinement/analysis/paper/flexible/g_r_${i}_${j}_f_o_h.png

        # 3Body angle plots
       #if [ $i -gt 8 ] ; then
       #    cp rigid/${i}_${j}_1/run${i}_${j}_1_3B_angles_layers.png ~/confinement/analysis/paper/rigid/ang_3B_${i}_${j}_r.png
       #    cp rigid/${i}_${j}_1/run${i}_${j}_1_3B_angles_layers_fit.png ~/confinement/analysis/paper/rigid/ang_3B_${i}_${j}_r.png
       #    cp ${i}_${j}_1/run${i}_${j}_1_3B_angles_layers.png ~/confinement/analysis/paper/flexible/ang_3B_${i}_${j}_f.png
       #    cp ${i}_${j}_1/run${i}_${j}_1_3B_angles_layers_fit.png ~/confinement/analysis/paper/flexible/ang_3B_${i}_${j}_f.png
       #else
       #    cp rigid/${i}_${j}_1/run${i}_${j}_1_3B_angles.png ~/confinement/analysis/paper/rigid/ang_3B_${i}_${j}_r.png
       #    cp ${i}_${j}_1/run${i}_${j}_1_3B_angles.png ~/confinement/analysis/paper/flexible/ang_3B_${i}_${j}_f.png
       #    cp rigid/${i}_${j}_1/run${i}_${j}_1_3B_angles_fit.png ~/confinement/analysis/paper/rigid/ang_3B_${i}_${j}_r.png
       #    cp ${i}_${j}_1/run${i}_${j}_1_3B_angles_fit.png ~/confinement/analysis/paper/flexible/ang_3B_${i}_${j}_f.png
       #fi

        # msd plots
       #cp rigid/${i}_${j}_1/msd_${i}_${j}_1_MSD3D.png ~/confinement/analysis/paper/rigid/msd_3D_${i}_${j}_r.png
       #cp rigid/${i}_${j}_1/msd__MSDYZ.png ~/Desktop/confinement_test/fig/msd_2D_${i}_${j}_r.png
       #cp ${i}_${j}_1/msd__MSDYZ.png ~/Desktop/confinement_test/fig/msd_2D_${i}_${j}_f.png
       #cp ${i}_${j}_1/msd_${i}_${j}_1_MSD3D.png ~/confinement/analysis/paper/flexible/msd_3D_${i}_${j}_f.png
       #cp ${i}_${j}_1/msd_${i}_${j}_1_MSDYZ.png ~/confinement/analysis/paper/flexible/msd_2D_${i}_${j}_f.png
        
    done
done

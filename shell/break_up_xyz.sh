#!/bin/bash

# script to move files for each replica of traj

sep=6
l=46
iter=2

ext=".velxyz"
ext=".xyz"
xyzf=`echo "traj_${sep}_${l}_${iter}${ext}"`
xyz=`echo "traj_${sep}_${l}"`

if [ $ext = ".xyz" ] ; then
    nat=`head -n1 $xyzf`
    nat=$((nat + 2)) ## for xyz file
    nsnaps=`grep -c Atom $xyzf`
elif [ $ext = ".velxyz" ]  ; then
    nat=`head -n4 $xyzf | tail -n 1`
    nat=$((nat + 9)) ## for vel file
    nsnaps=`grep -c TIMESTEP $xyzf`
fi

echo "$nsnaps"
echo "$nat"
head -n $nat $xyzf > first.xyz

incr=$((nat * 2500))
fst=$((nat * 5000))
scnd=$((fst + incr))
trd=$((scnd + incr))
fou=$((trd + incr))
lst=$((fou + incr))
fst=$((fst + 1))

echo "ext: ${ext} and range ${fst} ${scnd} ${trd} ${fou} ${lst}"
 sed -n ${fst},${scnd}p $xyzf  > q1.xyz
 sed -n $((scnd + 1)),${trd}p $xyzf  > q2.xyz
 sed -n $((trd + 1)),${fou}p $xyzf  > q3.xyz
 sed -n $((fou + 1)),${lst}p $xyzf  > q4.xyz

cat first.xyz q1.xyz > ${xyz}_q1${ext}
cat first.xyz q2.xyz > ${xyz}_q2${ext}
cat first.xyz q3.xyz > ${xyz}_q3${ext}
cat first.xyz q4.xyz > ${xyz}_q4${ext}

rm q1.xyz q2.xyz q3.xyz q4.xyz



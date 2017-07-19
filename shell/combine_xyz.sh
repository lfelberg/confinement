#!/bin/bash

# script to move files for each replica of traj

nat=`head -n1 *q1.xyz`
nat=$((nat + 2)) ## for xyz file
nsnaps=`grep -c Atom *q1.xyz`
nl_q3=`wc -l *q3.xyz | awk '{print $1}'`
nl_v3=`wc -l *q3.vol | awk '{print $1}'`
echo "$nat"

# cat together q2 and q3
q2=`ls *q2.xyz`
q2_vl=`ls *q2.vol`
q3=`ls *q3.xyz`
q3_vl=`ls *q3.vol`
nm=`echo ${q3/q3/h}`
nm_vl=`echo ${q3_vl/q3/h}`

sed -n $((nat + 1)),${nl_q3}p $q3  > q3.xyz
sed -n 2,${nl_v3}p $q3_vl  > q3.vol
cat ${q2} q3.xyz > ${nm}
cat ${q2_vl} q3.vol > ${nm_vl}

rm q3.xyz q3.vol

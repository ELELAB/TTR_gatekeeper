#!/bin/bash

j=0
z=10000
num_times=15
traj=../traj_centered_150ns.xtc
tpr=../md.tpr
ndx=../chains.ndx

# compute the RMSF over windows of width z (ps)
for ((i=1; i<=$num_times; i++)); do
g_rmsf -f $traj -s $tpr -n $ndx -res -o rmsf$i.xvg -b $j -e $z << eof
3
eof

xvg2octave rmsf$i.xvg 
j=$(($j+10000))
z=$(($z+10000))
done

# compute the average
paste rmsf1.oct rmsf2.oct rmsf3.oct rmsf4.oct rmsf5.oct rmsf6.oct rmsf7.oct rmsf8.oct rmsf9.oct rmsf10.oct rmsf11.oct rmsf12.oct rmsf13.oct rmsf14.oct rmsf15.oct | awk '{print $1, " ",$2, " ",$4, " ",$6, " ", $8, " ", $10, " ", $12, " ", $14, " ", $16, " ", $18, " ", $20}' > all.rmsf.oct
cat all.rmsf.oct | awk '{sum=0; for(i=2; i<=NF; i++){sum+=$i}; sum/=NF; print $1, " ", " ", sum}' > average_rmsf.oct
rm all.rmsf.oct

#!/bin/bash

#PBS -l nodes=1:ppn=32:xe
#PBS -l walltime=00:30:00
#PBS -q high

#cd $PBS_O_WORKDIR

export OMP_NUM_THREADS=1
aprun -n 1 -N 1 -j 1 -d $OMP_NUM_THREADS ./mom 

#for i in 1 2 4 8 16;
#do
#  for j in `seq 16 16`;
#  do
#    echo node_$i thread_$j
#    export OMP_NUM_THREADS=$j
#    aprun -n $i -N 1 -d $OMP_NUM_THREADS -j 1 ./mom > outt_cpu_${i}_${j}_j1
#    mv times.txt times_${i}_${j}_j1.txt
#  done
#done

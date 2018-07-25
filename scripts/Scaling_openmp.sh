#!/bin/bash


# Obtain the polynomial degree
echo $(date) >> scaling_openmp.txt

# Find the line number where the number of processes is specified

for nP in $( seq  1 8)
do
   echo 'Running with ' $nP ' threads'
   export OMP_NUM_THREADS=$nP

   # Run the integrator
   ./integrate
   walltime=$( grep 'Forward_Step_RK3' Timing.stats | awk '{ print $2 }' )
   echo $nP $walltime >> scaling_openmp.txt
   
done

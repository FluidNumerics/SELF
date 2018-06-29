#!/bin/bash


# Obtain the polynomial degree
pDeg=$( grep polyDeg runtime.params | awk '{gsub(/,$/,""); print $3}' )
nts=$( grep nTimeSteps runtime.params | awk '{gsub(/,$/,""); print $3}' )
echo $(date) >> scaling_openmp.txt

# Find the line number where the number of processes is specified

for nP in $( seq  1 8)
do
   echo 'Running with ' $nP ' threads'
   export OMP_NUM_THREADS=$nP

   # Run the integrator
   ./integrate_eu
   walltime=$( grep 'Accumulated Time' Timing.stats | awk 'BEGIN { FS = ":" } ; { print $2 }' )
   echo $nts $pDeg $nP $walltime >> scaling_openmp.txt
   
done

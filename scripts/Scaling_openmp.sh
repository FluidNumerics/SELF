#!/bin/bash

# Set the command used to invoke mpirun
mpirunner='/home/joe/Software/openmpi/bin/mpirun'

# Obtain the polynomial degree
pDeg=$( grep polyDeg runtime.params | awk '{gsub(/,$/,""); print $3}' )
nts=$( grep nTimeSteps runtime.params | awk '{gsub(/,$/,""); print $3}' )
echo $(date) >> scaling_openmp.txt

# Find the line number where the number of processes is specified

for nP in $( seq  1 4)
do
   echo 'Running with ' $nP ' threads'
   # modify geometry.params

   # commands for starting a fresh run
   # Decompose a mesh specified by either a SpecMeshFile or a Default Mesh
   ./GenerateFromSpecMesh
   export OMP_NUM_THREADS=$nP
   # Run the initializer
   ./initialize_sw

   # Run the integrator
   (time ./integrate_sw) 2> timing.txt
   walltime=$(head -n 2 timing.txt | awk '{print $2}')
   rm timing.txt

   echo $nts $pDeg $nP $walltime >> scaling_openmp.txt
   
   # Clean up old run files
   ./DeepClean.sh
done

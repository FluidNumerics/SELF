#!/bin/bash

# Set the command used to invoke mpirun
mpirunner='mpirun_pgi'

# Obtain the polynomial degree
pDeg=$( grep polyDeg runtime.params | awk '{gsub(/,$/,""); print $3}' )
nts=$( grep nTimeSteps runtime.params | awk '{gsub(/,$/,""); print $3}' )
echo $(date) >> scaling_mpi.txt

for j in {0..3}
  do  
   nP=$((2**j))
   echo 'Running with ' $nP 'processes'
   # modify geometry.params
   sed -i '/nProcZ/c\nProcZ = '$nP',' runtime.params

   mpicommand=$mpirunner' -np '$nP
   echo $mpicommand
   # commands for starting a fresh run
   # Decompose a mesh specified by either a SpecMeshFile or a Default Mesh
   echo '=============================================================='
   echo ' Generating mesh and domain decomposition '
   ./DecomposeStructuredHexMesh
   echo '=============================================================='
   echo ' Generating Initial Conditions '
   # Run the initializer
   $mpicommand ./initialize
   echo '=============================================================='
   
   echo ' Running Integrator !'
   # Run the integrator
   $mpicommand ./integrate
   walltime=$(more Timing.stats | grep "Accumulated Time" | awk 'BEGIN { FS = ":" }; {print $2}')

   echo $nts $pDeg $nP $walltime >> scaling_mpi.txt
   
   # Clean up old run files
   rm State.* box.* ExtComm.* mesh.*
   echo '=============================================================='

  done

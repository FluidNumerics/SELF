#!/bin/bash

# Set the command used to invoke mpirun
#mpirunner='/opt/pgi/linux86-64/2017/mpi/openmpi-1.10.2/bin/mpirun'
mpirunner='/usr/local/pgi/linux86-64/2016/mpi/openmpi-1.10.2/bin/mpirun'

# Obtain the polynomial degree
pDeg=$( grep polyDeg runtime.params | awk '{gsub(/,$/,""); print $3}' )
nts=$( grep nTimeSteps runtime.params | awk '{gsub(/,$/,""); print $3}' )
echo $(date) >> scaling_mpi.txt

   echo 'Running with 1 process'
   # modify geometry.params
   sed -i '/nProcX/c\nProcX = '1',' runtime.params
   sed -i '/nProcY/c\nProcY = '1',' runtime.params
   sed -i '/nProcZ/c\nProcZ = '1',' runtime.params
   sed -i '/nZElem/c\nZElem = '8',' runtime.params
   sed -i '/nYElem/c\nYElem = '8',' runtime.params
   sed -i '/nXElem/c\nXElem = '8',' runtime.params

   mpicommand=$mpirunner' -np 1'
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

   echo $nts $pDeg 1 $walltime >> scaling_mpi.txt
   
   # Clean up old run files
   rm State.* box.* ExtComm.* mesh.*
   echo '=============================================================='

   echo 'Running with 2 process'
   # modify geometry.params
   sed -i '/nProcX/c\nProcX = '2',' runtime.params
   sed -i '/nProcY/c\nProcY = '1',' runtime.params
   sed -i '/nProcZ/c\nProcZ = '1',' runtime.params
   sed -i '/nZElem/c\nZElem = '8',' runtime.params
   sed -i '/nYElem/c\nYElem = '8',' runtime.params
   sed -i '/nXElem/c\nXElem = '16',' runtime.params

   mpicommand=$mpirunner' -np 2'
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

   echo $nts $pDeg 2 $walltime >> scaling_mpi.txt
   
   echo 'Running with 4 process'
   # modify geometry.params
   sed -i '/nProcX/c\nProcX = '2',' runtime.params
   sed -i '/nProcY/c\nProcY = '2',' runtime.params
   sed -i '/nProcZ/c\nProcZ = '1',' runtime.params
   sed -i '/nZElem/c\nZElem = '8',' runtime.params
   sed -i '/nYElem/c\nYElem = '16',' runtime.params
   sed -i '/nXElem/c\nXElem = '16',' runtime.params

   mpicommand=$mpirunner' -np 4'
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

   echo $nts $pDeg 4 $walltime >> scaling_mpi.txt

   echo 'Running with 8 process'
   # modify geometry.params
   sed -i '/nProcX/c\nProcX = '2',' runtime.params
   sed -i '/nProcY/c\nProcY = '2',' runtime.params
   sed -i '/nProcZ/c\nProcZ = '2',' runtime.params
   sed -i '/nZElem/c\nZElem = '16',' runtime.params
   sed -i '/nYElem/c\nYElem = '16',' runtime.params
   sed -i '/nXElem/c\nXElem = '16',' runtime.params

   mpicommand=$mpirunner' -np 8'
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

   echo $nts $pDeg 8 $walltime >> scaling_mpi.txt
   rm State.* box.* ExtComm.* mesh.*

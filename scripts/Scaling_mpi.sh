#!/bin/bash

# Set the command used to invoke mpirun
mpirunner='/home/joe/Software/openmpi/bin/mpirun'

# Obtain the polynomial degree
pDeg=$( grep polyDeg runtime.params | awk '{gsub(/,$/,""); print $3}' )
nts=$( grep nTimeSteps runtime.params | awk '{gsub(/,$/,""); print $3}' )
echo $(date) >> scaling_mpi.txt

# Find the line number where the number of processes is specified

   echo 'Running with 1 processes'
   # modify geometry.params

   # commands for starting a fresh run
   # Decompose a mesh specified by either a SpecMeshFile or a Default Mesh
   ./GenerateFromSpecMesh

   # Run the initializer
   ./initialize_sw

   # Run the integrator
   (time ./integrate_sw) 2> timing.txt
   walltime=$(head -n 2 timing.txt | awk '{print $2}')
   rm timing.txt

   echo $nts $pDeg '1' $walltime >> scaling.txt
   
   # Clean up old run files
   ./DeepClean.sh
   
for nP in $( seq  2 4)
do
   echo 'Running with ' $nP 'processes'
   # modify geometry.params
   sed -i '/nProc/c\nProc = '$nP',' geometry.params

   mpicommand=$mpirunner' -np '$nP
   echo $mpicommand
   # commands for starting a fresh run
   # Decompose a mesh specified by either a SpecMeshFile or a Default Mesh
   ./DecomposeQuadMesh

   # Run the initializer
   $mpicommand ./initialize_sw_mpi

   # Run the integrator
   (time $mpicommand ./integrate_sw_mpi) 2> timing.txt
   walltime=$(head -n 2 timing.txt | awk '{print $2}')
   rm timing.txt

   # Convert pickup files to tecplot
  # $mpicommand ./PickupToTecplot_mpi

   # "Glue" the mpi output to one file per time step
  # ./GlueMPITecplot.sh

   echo $nts $pDeg $nP $walltime >> scaling_mpi.txt
   
   # Clean up old run files
   ./DeepClean.sh
done

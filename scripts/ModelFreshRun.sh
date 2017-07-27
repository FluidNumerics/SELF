#!/bin/bash

# Set the command used to invoke mpirun
mpirunner='/home/joe/Software/openmpi/bin/mpirun'

#obtain the number of MPI processes consistent with the domain decomposition
nP=$( grep nProc geometry.params | awk '{gsub(/,$/,""); print $3}' )

mpicommand=$mpirunner' -np '$nP

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
$mpicommand ./PickupToTecplot_mpi

# "Glue" the mpi output to one file per time step
./GlueMPITecplot.sh

echo 'Wall Time : '$walltime

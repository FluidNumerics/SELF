#!/bin/bash

# Set the command used to invoke mpirun
mpirunner='/home/joe/Software/openmpi/bin/mpirun'

#obtain the number of MPI processes consistent with the domain decomposition
nP=$( grep nProc geometry.params | awk '{gsub(/,$/,""); print $3}' )

mpicommand=$mpirunner' -np '$nP

# Run the integrator
time $mpicommand ./integrate_sw_mpi

# Convert pickup files to tecplot
$mpicommand ./PickupToTecplot_mpi

# "Glue" the mpi output to one file per time step
./GlueMPITecplot.sh


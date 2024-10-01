#!/bin/bash


for f in *.f90; do
	echo $f
	sed -i 's/SELF_MPI/SELF_DomainDecomposition/g' $f
done

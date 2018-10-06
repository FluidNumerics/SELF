#!/bin/bash

WORKDIR=/home/joe_fluidnumerics_com/SELF-Fluids
INSTDIR=/home/joe_fluidnumerics_com/apps/self-fluids/dev

cd $WORKDIR


# PGI 18.4
module load pgi/18.4 hdf5/1.10.2/serial/pgi/18.4
${WORKDIR}/configure --enable-tecplot \
                     --enable-timing \
                     --enable-diagnostics \
                     --prefix=${INSTDIR}/serial/pgi/18.4
make 
make install
make distclean
module purge

# PGI 18.4 + OpenMP
module load pgi/18.4 hdf5/1.10.2/serial/pgi/18.4
${WORKDIR}/configure --enable-tecplot \
                     --enable-timing \
                     --enable-diagnostics \
                     --enable-openmp \
                     --prefix=${INSTDIR}/openmp/pgi/18.4
make 
make install
make distclean
module purge

# PGI 18.4 + OpenMPI 3.1.2
module load pgi/18.4 openmpi/3.1.2/pgi/18.4 hdf5/1.10.2/parallel/pgi/18.4/openmpi/3.1.2
${WORKDIR}/configure --enable-tecplot \
                     --enable-mpi \
                     --enable-timing \
                     --enable-diagnostics \
                     --prefix=${INSTDIR}/mpi/pgi/18.4/openmpi/3.1.2 \
make 
make install
make distclean
module purge

# PGI 18.4 + OpenMP + OpenMPI 3.1.2
module load pgi/18.4 openmpi/3.1.2/pgi/18.4 hdf5/1.10.2/parallel/pgi/18.4/openmpi/3.1.2
${WORKDIR}/configure --enable-tecplot \
                     --enable-mpi \
                     --enable-timing \
                     --enable-diagnostics \
                     --enable-openmp \
                     --prefix=${INSTDIR}/openmp+mpi/pgi/18.4/openmpi/3.1.2 \
make 
make install
make distclean
module purge

# PGI 18.4 (CUDA)
module load pgi/18.4 hdf5/1.10.2/serial/pgi/18.4
${WORKDIR}/configure --enable-tecplot \
                     --enable-cuda \
                     --enable-timing \
                     --enable-diagnostics \
                     --prefix=${INSTDIR}/cuda/pgi/18.4 \
                     GPU_ARCH=cc70
make 
make install
make distclean
module purge

# PGI 18.4 (CUDA) + OpenMPI 3.1.2
module load pgi/18.4 openmpi/3.1.2/pgi/18.4 hdf5/1.10.2/parallel/pgi/18.4/openmpi/3.1.2
${WORKDIR}/configure --enable-tecplot \
                     --enable-cuda \
                     --enable-mpi \
                     --enable-timing \
                     --enable-diagnostics \
                     --prefix=${INSTDIR}/mpi+cuda/pgi/18.4/openmpi/3.1.2 \
                     GPU_ARCH=cc70
make 
make install
make distclean
module purge


#GCC 8.1.0 (Serial)
module load gcc/8.1.0 hdf5/1.10.2/serial/gcc/8.1.0
${WORKDIR}/configure --enable-tecplot \
                     --enable-timing \
                     --enable-diagnostics \
                     --prefix=${INSTDIR}/serial/gcc/8.1.0
make 
make install
make distclean
module purge

#GCC 8.1.0 + OpenMPI 3.1.2
module load gcc/8.1.0 openmpi/3.1.2/gcc/8.1.0 hdf5/1.10.2/parallel/gcc/8.1.0/openmpi/3.1.2
${WORKDIR}/configure --enable-tecplot \
                     --enable-mpi \
                     --enable-timing \
                     --enable-diagnostics \
                     --prefix=${INSTDIR}/mpi/gcc/8.1.0/openmpi/3.1.2
make 
make install
make distclean
module purge

# GCC 8.1.0 + OpenMP
module load gcc/8.1.0 hdf5/1.10.2/serial/gcc/8.1.0
${WORKDIR}/configure --enable-tecplot \
                     --enable-timing \
                     --enable-diagnostics \
                     --enable-openmp \
                     --prefix=${INSTDIR}/openmp/gcc/8.1.0
make 
make install
make distclean
module purge

#GCC 8.1.0 + OpenMP + OpenMPI 3.1.2
module load gcc/8.1.0 openmpi/3.1.2/gcc/8.1.0 hdf5/1.10.2/parallel/gcc/8.1.0/openmpi/3.1.2
${WORKDIR}/configure --enable-tecplot \
                     --enable-mpi \
                     --enable-timing \
                     --enable-diagnostics \
                     --enable-openmp \
                     --prefix=${INSTDIR}/openmp+mpi/gcc/8.1.0/openmpi/3.1.2
make 
make install
make distclean
module purge


# DOUBLE PRECISION

# PGI 18.4
module load pgi/18.4 hdf5/1.10.2/serial/pgi/18.4
${WORKDIR}/configure --enable-tecplot \
                     --enable-timing \
                     --enable-diagnostics \
                     --enable-double \
                     --prefix=${INSTDIR}/serial_double/pgi/18.4
make 
make install
make distclean
module purge

# PGI 18.4 + OpenMP
module load pgi/18.4 hdf5/1.10.2/serial/pgi/18.4
${WORKDIR}/configure --enable-tecplot \
                     --enable-timing \
                     --enable-diagnostics \
                     --enable-double \
                     --enable-openmp \
                     --prefix=${INSTDIR}/openmp_double/pgi/18.4
make 
make install
make distclean
module purge

# PGI 18.4 + OpenMPI 3.1.2
module load pgi/18.4 openmpi/3.1.2/pgi/18.4 hdf5/1.10.2/parallel/pgi/18.4/openmpi/3.1.2
${WORKDIR}/configure --enable-tecplot \
                     --enable-mpi \
                     --enable-timing \
                     --enable-diagnostics \
                     --enable-double \
                     --prefix=${INSTDIR}/mpi_double/pgi/18.4/openmpi/3.1.2 \
make 
make install
make distclean
module purge

# PGI 18.4 + OpenMP + OpenMPI 3.1.2
module load pgi/18.4 openmpi/3.1.2/pgi/18.4 hdf5/1.10.2/parallel/pgi/18.4/openmpi/3.1.2
${WORKDIR}/configure --enable-tecplot \
                     --enable-mpi \
                     --enable-timing \
                     --enable-diagnostics \
                     --enable-double \
                     --enable-openmp \
                     --prefix=${INSTDIR}/openmp+mpi_double/pgi/18.4/openmpi/3.1.2 \
make 
make install
make distclean
module purge

# PGI 18.4 (CUDA)
module load pgi/18.4 hdf5/1.10.2/serial/pgi/18.4
${WORKDIR}/configure --enable-tecplot \
                     --enable-cuda \
                     --enable-timing \
                     --enable-diagnostics \
                     --enable-double \
                     --prefix=${INSTDIR}/cuda_double/pgi/18.4 \
                     GPU_ARCH=cc70
make 
make install
make distclean
module purge

# PGI 18.4 (CUDA) + OpenMPI 3.1.2
module load pgi/18.4 openmpi/3.1.2/pgi/18.4 hdf5/1.10.2/parallel/pgi/18.4/openmpi/3.1.2
${WORKDIR}/configure --enable-tecplot \
                     --enable-cuda \
                     --enable-mpi \
                     --enable-double \
                     --enable-timing \
                     --enable-diagnostics \
                     --prefix=${INSTDIR}/mpi+cuda_double/pgi/18.4/openmpi/3.1.2 \
                     GPU_ARCH=cc70
make 
make install
make distclean
module purge


#GCC 8.1.0 (Serial)
module load gcc/8.1.0 hdf5/1.10.2/serial/gcc/8.1.0
${WORKDIR}/configure --enable-tecplot \
                     --enable-timing \
                     --enable-diagnostics \
                     --enable-double \
                     --prefix=${INSTDIR}/serial_double/gcc/8.1.0
make 
make install
make distclean
module purge

#GCC 8.1.0 + OpenMPI 3.1.2
module load gcc/8.1.0 openmpi/3.1.2/gcc/8.1.0 hdf5/1.10.2/parallel/gcc/8.1.0/openmpi/3.1.2
${WORKDIR}/configure --enable-tecplot \
                     --enable-mpi \
                     --enable-timing \
                     --enable-diagnostics \
                     --enable-double \
                     --prefix=${INSTDIR}/mpi_double/gcc/8.1.0/openmpi/3.1.2
make 
make install
make distclean
module purge

# GCC 8.1.0 + OpenMP
module load gcc/8.1.0 hdf5/1.10.2/serial/gcc/8.1.0
${WORKDIR}/configure --enable-tecplot \
                     --enable-timing \
                     --enable-diagnostics \
                     --enable-openmp \
                     --enable-double \
                     --prefix=${INSTDIR}/openmp_double/gcc/8.1.0
make 
make install
make distclean
module purge

#GCC 8.1.0 + OpenMP + OpenMPI 3.1.2
module load gcc/8.1.0 openmpi/3.1.2/gcc/8.1.0 hdf5/1.10.2/parallel/gcc/8.1.0/openmpi/3.1.2
${WORKDIR}/configure --enable-tecplot \
                     --enable-mpi \
                     --enable-timing \
                     --enable-diagnostics \
                     --enable-openmp \
                     --enable-double \
                     --prefix=${INSTDIR}/openmp+mpi_double/gcc/8.1.0/openmpi/3.1.2
make 
make install
make distclean
module purge

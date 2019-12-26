# PGI+OpenMPI Dependency Container

To rebuild Docker file, use HPC container maker
```
hpccm --format=docker --recipe=pgi-openmpi-hdf5.py  --userarg pgi_eula_accept=yes > Dockerfile
```

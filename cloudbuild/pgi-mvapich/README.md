# PGI+MVAPICH2 Dependency Container

To rebuild Docker file, use HPC container maker
```
hpccm --format=docker --recipe=pgi-mvapich2-hdf5.py  --userarg centos=true pgi_eula_accept=yes > Dockerfile
```

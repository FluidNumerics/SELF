# PGI Dependency Container

To rebuild Docker file, use HPC container maker
```
hpccm --format=docker --recipe=pgi-hdf5.py  --userarg centos=true pgi_eula_accept=yes > Dockerfile
```

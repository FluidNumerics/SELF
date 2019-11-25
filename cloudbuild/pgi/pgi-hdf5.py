"""
HPC Base image

Contents:
  CUDA version 10.0
  HDF5 version 1.10.5
  PGI compilers version 19.4
  Python 2 and 3 (upstream)
"""
# pylint: disable=invalid-name, undefined-variable, used-before-assignment

# The PGI End-User License Agreement (https://www.pgroup.com/doc/LICENSE)
# must be accepted.
pgi_eula=False
if USERARG.get('pgi_eula_accept', False):
  pgi_eula=True
else:
  raise RuntimeError('PGI EULA not accepted. To accept, use "--userarg pgi_eula_accept=yes"\nSee PGI EULA at https://www.pgroup.com/doc/LICENSE')

# Choose between either Ubuntu 16.04 (default) or CentOS 7
# Add '--userarg centos=true' to the command line to select CentOS
devel_image = 'nvidia/cuda:10.1-devel-ubuntu18.04'
runtime_image = 'nvidia/cuda:10.1-runtime-ubuntu18.04'
if USERARG.get('centos', False):
    devel_image = 'nvidia/cuda:10.1-devel-centos7'
    runtime_image = 'nvidia/cuda:10.1-runtime-centos7'

######
# Devel stage
######

Stage0 += comment(__doc__, reformat=False)

Stage0 += baseimage(image=devel_image, _as='devel')

# Python
Stage0 += python()

# PGI compilers
compiler = pgi(eula=pgi_eula, version='19.4')
Stage0 += compiler

# HDF5
Stage0 += hdf5(version='1.10.5', mpi=False, toolchain=compiler.toolchain)

# Metis
Stage0 += shell(commands=['yum install -y cmake'])
Stage0 += shell(commands=['mkdir -p /var/tmp',
                          'wget -q -nc --no-check-certificate -P /var/tmp http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz',
                          'tar -xzf /var/tmp/metis-5.1.0.tar.gz -C /var/tmp',
                          'cd /var/tmp/metis-5.1.0',
                          'make config prefix=/usr/local/metis',
                          'make install'])
Stage0 += environment(variables={'LIB_METIS':'/usr/local/metis/libmetis.a'})

######
# Runtime image
######

Stage1 += baseimage(image=runtime_image)

Stage1 += Stage0.runtime(_from='devel')

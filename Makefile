# SELF Makefile
#
# **Make system inspired by and adopted from https://github.com/LKedward/focal**
#
#-------------------
# SELF Dependencies
#-------------------
#   
#   * FLAP (https://github.com/szaghi/FLAP) ( For dev builds )
#   * feq-parse (https://github.com/fluidnumerics/feq-parse) ( For dev builds )
#   * HDF5
#
#-----------------------
# Environment variables
#-----------------------
#
#   Installation Options
#     SELF_PREFIX            Set the path to install SELF (Default: /opt/self)
#     BUILD                  Set the type of build. (Default: dev | Options : dev, release)
#     SELF_DIR               Set the path to the repository root ( Default $(shell pwd) )
#
#   Precision Options
#     PREC                   Set the floating point precision in SELF. (Default: single | Options : single, double)
#
#   Compiler and Target Hardware Options
#     MPIFC                  Set the full path to a MPI fortran compiler. (Default: mpifort)
#     GPU_PLATFORM           Set the target GPU platform. (Default: gfx900 | Options gfx900, sm35, sm50, sm70) 
#     HIP_PLATFORM           Set the target vendor. (Default: amd | Options amd, nvidia) 
#
#   Dependency Options
#     ROCM_DIR               Set the path to ROCm installation (Default: /opt/rocm)
#     CUDA_PATH              Set the path to CUDA installation (Default: /usr/local/cuda) (Needed if GPU_TARGET is an Nvidia GPU)
#     SELF_JSONF_LIBS        Set the linker flags for json-fortran (Default: -L/opt/view/lib -ljsonfortran)
#     SELF_JSONF_INC         Set the includes flags for json-fortran (Default: -I/opt/view/include)
#     SELF_FEQPARSE_LIBS     Set the linker flags for feq-parse (Default: -L/opt/view/lib -lfeqparse)
#     SELF_FEQPARSE_INC      Set the includes flags for feq-parse (Default: -I/opt/view/include)
#     SELF_HDF5_LIBS         Set the linker flags for hdf5 (Default: -L/opt/view/lib -lhdf5_fortran -lhdf5 -lz -lm)
#     SELF_HDF5_INC          Set the includes flags for hdf5 (Default: -I/opt/view/include/shared)
#     SELF_GPU_LIBS          Set the linker flags for either HIP(Default) or CUDA (Default: -L/opt/rocm/lib -lamdhip64)
#     SELF_GPU_INC           Set the includes flags for either HIP(Default) or CUDA (Default: -I/opt/rocm/include)
#
# ================================================================================================================================================= #

SELF_PREFIX ?= /opt/self
SELF_DIR ?= .

# Build Target
install: all
	mkdir -p $(SELF_PREFIX)/bin
	mkdir -p $(SELF_PREFIX)/lib
	mkdir -p $(SELF_PREFIX)/include
	mkdir -p $(SELF_PREFIX)/test
	mkdir -p $(SELF_PREFIX)/util
	mv -f $(SELF_DIR)/build/include/* $(SELF_PREFIX)/include/
	mv -f $(SELF_DIR)/build/lib/* $(SELF_PREFIX)/lib/
	mv -f $(SELF_DIR)/build/bin/* $(SELF_PREFIX)/bin/
	mv -f $(SELF_DIR)/util/* $(SELF_PREFIX)/util/
	cp -r $(SELF_DIR)/etc $(SELF_PREFIX)/etc

all: self libself

include ${SELF_DIR}/make.include

clean: self_clean


.PHONY: all clean

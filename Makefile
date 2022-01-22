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
#
#   Precision Options
#     PREC                   Set the floating point precision in SELF. (Default: single | Options : single, double)
#
#   Compiler and Target Hardware Options
#     HIPFC                  Set the full path to hipfc. (Default: hipfc)
#     HIPFORT_COMPILER       Set the Fortran compiler used by hipfort. (Default: gfortran | Options: gfortran)
#     HIPFORT_GPU            Set the target GPU platform. (Default: gfx900 | Options gfx900, sm35, sm50, sm70) 
#     ROCM_DIR               Set the path to ROCm installation (Default: /opt/rocm)
#     CUDA_PATH              Set the path to CUDA installation (Default: /usr/local/cuda) (Needed if GPU_TARGET is an Nvidia GPU)
#
#   Dependency Options
#     SELF_JSONF_LIBS        Set the linker flags for json-fortran (Default: -L/opt/view/lib -ljsonfortran)
#     SELF_JSONF_INC         Set the includes flags for json-fortran (Default: -I/opt/view/include)
#     SELF_FEQPARSE_LIBS     Set the linker flags for feq-parse (Default: -L/opt/view/lib -lfeqparse)
#     SELF_FEQPARSE_INC      Set the includes flags for feq-parse (Default: -I/opt/view/include)
#     SELF_FLAP_LIBS         Set the linker flags for FLAP (Default: -L/opt/view/lib/ -lFLAP -lFACE -lPENF) 
#     SELF_FLAP_INC          Set the includes flags for FLAP (Default: -I/opt/view/include/FLAP -I/opt/view/include/PENF -I/opt/view/include/FACE)
#     SELF_HDF5_LIBS         Set the linker flags for hdf5 (Default: -L/opt/view/lib -lhdf5_fortran -lhdf5 -lz -lm)
#     SELF_HDF5_INC          Set the includes flags for hdf5 (Default: -I/opt/view/include/shared)
#
# ================================================================================================================================================= #

SELF_PREFIX ?= /opt/self
SELF_DIR ?= .

# Build Target
install: all
	mkdir -p $(SELF_PREFIX)
	mv $(SELF_DIR)/*.mod $(SELF_DIR)/include/
	mv $(SELF_DIR)/include $(SELF_DIR)/lib $(SELF_DIR)/obj $(SELF_DIR)/bin $(SELF_PREFIX)

all: self

include ${SELF_DIR}/make.include

clean: self_clean

# test:


.PHONY: all clean

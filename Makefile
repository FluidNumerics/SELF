# SELF Makefile
#
# SELF Dependencies
#   
#   * FLAP (https://github.com/szaghi/FLAP) ( For dev builds )
#   * feq-parse (https://github.com/fluidnumerics/feq-parse) ( For dev builds )
# 
#
# Make system inspired by and adopted from https://github.com/LKedward/focal
#
# Environment variables
#
#   SELF_PREFIX            Set the path to install SELF (make install)
#   BUILD                  Set the type of build. One of "dev" or "release"
#
#   SELF_FEQPARSE_LIBS     Set the linker flags for feq-parse
#   SELF_FEQPARSE_INC      Set the includes flags for feq-parse
#   SELF_FLAP_LIBS         Set the linker flags for FLAP
#   SELF_FLAP_INC          Set the includes flags for FLAP
#

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

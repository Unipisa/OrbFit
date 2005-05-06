# Makefile options for HP PA RISC, debugging

# Fortran compiler
FC=f77
# Options for Fortran compiler:
FFLAGS=-C -g -K -I../include
# "ranlib" command: if it is not needed, use "RANLIB=touch"
RANLIB=ranlib
VPATH=../include

# Makefile options for NAG FORTRAN90 compiler, debugging

# Fortran compiler
FC=f95
# Options for Fortran compiler:
FFLAGS=-C -g90  -I../include
# "ranlib" command: if it is not needed, use "RANLIB=touch"
RANLIB=ranlib
VPATH=../include

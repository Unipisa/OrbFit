# Makefile options for NAG FORTRAN90 compiler, debugging

# Fortran compiler
FC=f95
# Options for Fortran compiler:
FFLAGS=-O -pg  -I../include
# "ranlib" command: if it is not needed, use "RANLIB=touch"
RANLIB=ranlib
VPATH=../include

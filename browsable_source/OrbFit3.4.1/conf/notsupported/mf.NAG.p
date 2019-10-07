# Makefile options for NAG f95 compiler, profiling

# Fortran compiler
FC=f95
# Options for Fortran compiler:
FFLAGS=-O -pg -dcfuns -I../include
# "ranlib" command: if it is not needed, use "RANLIB=touch"
RANLIB=ranlib
VPATH=../include

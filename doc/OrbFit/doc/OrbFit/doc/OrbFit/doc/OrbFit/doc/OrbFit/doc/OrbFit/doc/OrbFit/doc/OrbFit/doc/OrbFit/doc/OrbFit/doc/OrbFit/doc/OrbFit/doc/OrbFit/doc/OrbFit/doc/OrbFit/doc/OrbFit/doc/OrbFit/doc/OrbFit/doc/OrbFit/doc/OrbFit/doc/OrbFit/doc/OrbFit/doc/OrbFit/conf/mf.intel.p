# Makefile options for the INTEL compiler, optimized and profiling

# Fortran compiler
FC=ifc
# Options for Fortran compiler:
FFLAGS= -cm -O3 -mp1 -pg -Vaxlib -I../include 
# "ranlib" command: if it is not needed, use "RANLIB=touch"
RANLIB=ranlib
VPATH=../include
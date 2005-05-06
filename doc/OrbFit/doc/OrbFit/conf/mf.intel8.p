# Makefile options for the INTEL compiler, optimized and profiling

# Fortran compiler
FC=ifort
# Options for Fortran compiler:
FFLAGS= -cm -O -mp1 -pg -save -assume byterecl -Vaxlib -I../include 
# "ranlib" command: if it is not needed, use "RANLIB=touch"
RANLIB=ranlib
VPATH=../include

# Makefile options for the INTEL compiler, not optimized, high accuracy

# Fortran compiler
FC=ifort
# Options for Fortran compiler for debugging:
FFLAGS=  -cm -g -pc80 -mp -CB -traceback -save -assume byterecl -Vaxlib -I../include 
# "ranlib" command: if it is not needed, use "RANLIB=touch"
RANLIB=ranlib
VPATH=../include

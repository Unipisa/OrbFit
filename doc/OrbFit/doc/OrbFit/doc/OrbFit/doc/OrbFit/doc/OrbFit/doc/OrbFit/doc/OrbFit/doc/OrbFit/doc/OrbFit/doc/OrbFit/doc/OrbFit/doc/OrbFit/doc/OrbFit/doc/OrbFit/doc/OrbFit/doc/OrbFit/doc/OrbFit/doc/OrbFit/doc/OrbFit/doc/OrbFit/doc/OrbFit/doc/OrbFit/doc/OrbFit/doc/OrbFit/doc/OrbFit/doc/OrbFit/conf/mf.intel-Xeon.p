# Makefile options for the INTEL compiler, optimized

# Fortran compiler
FC=ifc
# Options for Fortran compiler:
FFLAGS= -cm -O3 -mp1 -tpp7 -xW -pg -Vaxlib -I../include
# "ranlib" command: if it is not needed, use "RANLIB=touch"
RANLIB=ranlib
VPATH=../include

# Makefile options for Portland PGI compiler, profiling

# Fortran compiler
FC=pgf90
# Options for Fortran compiler:
FFLAGS= -fast -Mprof=func  -I../include
# "ranlib" command: if it is not needed, use "RANLIB=touch"
RANLIB=ranlib
VPATH=../include

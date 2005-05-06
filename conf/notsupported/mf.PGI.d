# Makefile options for Portland PGI compiler, debugging

# Fortran compiler
FC=pgf90
# Options for Fortran compiler:
FFLAGS= -g -C -Mstandard  -I../include
# "ranlib" command: if it is not needed, use "RANLIB=touch"
RANLIB=ranlib
VPATH=../include

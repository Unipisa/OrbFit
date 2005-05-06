# Makefile options for GNU gfortran compiler, debugging

# Fortran compiler
FC=/usr/local/gfortran/irun/bin/gfortran 
# Options for Fortran compiler:
FFLAGS= -static -g   -I../include
# "ranlib" command: if it is not needed, use "RANLIB=touch"
RANLIB=ranlib
VPATH=../include

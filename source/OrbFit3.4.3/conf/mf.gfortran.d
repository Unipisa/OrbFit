# Makefile options for GNU gfortran compiler, debugging

# Fortran compiler
FC=/home/bedini/gcc-4_1-local/bin/gfortran-4.1
# Options for Fortran compiler:
FFLAGS= -static -g   -I../include
# "ranlib" command: if it is not needed, use "RANLIB=touch"
RANLIB=ranlib
VPATH=../include
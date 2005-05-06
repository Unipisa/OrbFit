# Makefile options for NAG f95 compiler, debugging

# Fortran compiler
FC=f95
# Options for Fortran compiler:
FFLAGS=-C -w=obs -dcfuns -g -gline -I../include
# "ranlib" command: if it is not needed, use "RANLIB=touch"
RANLIB=ranlib
VPATH=../include

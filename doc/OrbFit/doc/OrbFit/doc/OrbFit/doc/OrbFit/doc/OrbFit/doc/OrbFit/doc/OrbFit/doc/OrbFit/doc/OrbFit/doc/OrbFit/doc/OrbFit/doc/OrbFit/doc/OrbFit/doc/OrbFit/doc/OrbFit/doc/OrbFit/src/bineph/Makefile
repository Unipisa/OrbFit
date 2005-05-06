include ../make.flags
include modules.flags

all: bineph.x

# Executables

BINEPH=bineph.o rdoptb.o ra15b.o forceb.o

bineph.x: $(BINEPH) ../lib/libsuit.a 
	$(FC) $(FFLAGS) -o bineph.x $(BINEPH) ../lib/libsuit.a 

depend:
	../../bin/mkdep90 -I../include -L../suit *.f *.f90 *.h *.h90 > make.dep
	../../bin/mkmoddep -f ../.compiler ../suit > modules.flags

clean:
	rm -f *.x *.o *.mod lib*.a *% *~ core work.pc*

distclean: clean
	../lib/deinst_src

# Windows targets
win: ..\..\bin\bineph.exe

..\..\bin\bineph.exe: $(BINEPH) ..\lib\suit.lib
	$(FC) $(FFLAGS) /link /out:$@ $**

windist:
	@del *.o
	@del *.pdb
	@del ..\..\bin\*.pdb
	@del ..\..\bin\*.ilk

winclean:windist
	@del ..\..\bin\bineph.*

# Dependencies
include make.dep

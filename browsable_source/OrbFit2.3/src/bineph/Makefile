.f.o:
	$(FC) $(FFLAGS) -c $*.f -o $*.o 

include ../make.flags

all: bineph.x

# Executables

BINEPH=bineph.o rdoptb.o ra15b.o forceb.o

bineph.x: $(BINEPH) ../lib/libsuit.a 
	$(FC) $(FFLAGS) -o bineph.x $(BINEPH) ../lib/libsuit.a 

depend:
	./mkdep *.f *.h ../include/*.h > make.dep

clean:
	rm -f *.x *.o lib*.a *% *~ core work.pc*

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

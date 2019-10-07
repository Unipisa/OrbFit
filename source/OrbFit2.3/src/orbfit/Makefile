.f.o:
	$(FC) $(FFLAGS) -c $*.f -o $*.o 

include ../make.flags

all: orbfit.x

# Auxiliary routines (internal to orbfit)
ORBAUX=rdopto.o rdoptf.o rdopti.o rdopte.o srtobs.o seltid.o ustsp.o \
       ofiobs.o ofiorb.o ofinod.o ofofit.o ofiden.o ofephe.o ofprop.o \
       ofinip.o ofclrf.o

ORBFIT=orbfit.o $(ORBAUX)

ORBFITDEPS=$(ORBFIT) ../lib/libprop.a ../lib/libmoid.a ../lib/libgauss.a ../lib/libsuit.a

orbfit.x: $(ORBFITDEPS) 
	$(FC) $(FFLAGS) -o orbfit.x $(ORBFITDEPS)

depend:
	./mkdep *.f *.h ../include/*.h > make.dep

clean:
	rm -f *.x *.o lib*.a *% *~ core work.pc*

distclean: clean
	../lib/deinst_src

# Windows targets
win: ..\..\bin\orbfit.exe

..\..\bin\orbfit.exe: $(ORBFIT) ..\lib\prop.lib ..\lib\moid.lib ..\lib\suit.lib ..\lib\gauss.lib
	$(FC) $(FFLAGS) /link /out:$@ $**

windist:
	@del *.o
	@del *.pdb
	@del ..\..\bin\*.pdb
	@del ..\..\bin\*.ilk

winclean:windist
	@del ..\..\bin\orbfit.*

# Dependencies
include make.dep

.f.o:
	$(FC) $(FFLAGS) -c $*.f -o $*.o 

include ../make.flags

all: fitobs.x catpro.x

#packages

FIT=fgauss.o stacop.o fobpre.o

FIDENT=fident.o identsub.o

FINOUT=finopt.o finobs.o finele.o

FMULTI=mult_subs.o

FITOBS= fitobs.o $(FIT) $(FIDENT) $(FINOUT) $(FMULTI)

LIBDEPS= ../lib/libprop.a ../lib/libmoid.a ../lib/libgauss.a ../lib/libsuit.a

# Executables

fitobs.x: $(FITOBS) $(LIBDEPS)
	$(FC) $(FFLAGS) -o fitobs.x  $(FITOBS) $(LIBDEPS)

catpro.x: catpro.o cinopt.o  $(LIBDEPS)
	$(FC) $(FFLAGS) catpro.o cinopt.o $(LIBDEPS) -o catpro.x

depend:
	./mkdep *.f *.h ../include/*.h > make.dep

clean:
	-(rm -f *.x *.o *% *~ core work.pc*)

distclean: clean
	../lib/deinst_src

# Windows targets
win: ..\..\bin\fitobs.exe

..\..\bin\fitobs.exe:$(FITOBS) ..\lib\prop.lib ..\lib\moid.lib ..\lib\suit.lib ..\lib\gauss.lib
	$(FC) $(FFLAGS) /link /out:$@ $**

windist:
	@del *.o
	@del *.pdb
	@del ..\..\bin\*.pdb
	@del ..\..\bin\*.ilk

winclean:windist
	@del ..\..\bin\fitobs.*

# Dependencies
include make.dep

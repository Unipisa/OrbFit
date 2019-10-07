include ../make.flags

all: libmoid.a statcomp.x mincomp.x mintest.x

STATCOMP= compute_minima.o mutualrefcha.o aical.o comp_heart.o matrixdat.o \
	clevcompsylv.o rvfft.o decode_out.o code_inp.o irvfft.o solv20.o \
	solvesystem.o hess.o d2eval.o choosedeg.o sortd2.o 

CDET= cdetcomp.o zgetf2.o izamax.o dcabs1.o zscal.o \
	zswap.o zgeru.o xerbla.o

# library to compute stationary points of D2
libmoid.a: $(STATCOMP) $(CDET)
	ar r ../lib/libmoid.a $(STATCOMP) $(CDET)
	$(RANLIB) ../lib/libmoid.a

LIBDEPS= ../lib/libprop.a ../lib/libmoid.a ../lib/libgauss.a ../lib/libsuit.a

# main program
statcomp.x: statcomp.o optpro.o writestat.o compute_statpts.o $(LIBDEPS)
	$(FC) $(FFLAGS) statcomp.o optpro.o writestat.o  compute_statpts.o \
	$(LIBDEPS)  -o statcomp.x

mincomp.x:mincomp.o optpro.o writemin.o $(LIBDEPS)
	$(FC) $(FFLAGS) mincomp.o optpro.o writemin.o $(LIBDEPS) \
	-o mincomp.x

mintest.x:mintest.o $(LIBDEPS)
	$(FC) $(FFLAGS) mintest.o optpro.o writemin.o $(LIBDEPS) \
	-o mintest.x


depend:
	./mkdep *.f *.h ../include/*.h > make.dep

clean:
	-(rm -f *.x *.o ../lib/libmoid.a *% *~  core work.pc*)

distclean: clean
	../lib/deinst_src

neotar:	clean
	(cd ../ ; tar -cf neodys.tar neodys/)

# Windows targets
win: ..\lib\moid.lib 

..\lib\moid.lib:$(STATCOMP) $(CDET)
	lib /out:$@ $**

winclean:
	@del *.o
	@del *.pdb
	@del ..\lib\moid.lib


# Dependencies
include make.dep


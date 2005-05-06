include ../make.flags
include modules.flags

all: libmoid.a min_comp.x stat_comp.x stat_search.x

STATCOMP= writemin.o writestat3.o compute_minima.o compute_statpts.o \
	mutualrefcha.o aical.o comp_heart_rot.o matrixdat_rot.o \
	clevcompsylv.o rvfft.o decode_out.o code_inp.o irvfft.o \
	solvpoly.o solvesystem_rot.o hess.o int_eval.o d2eval.o \
	choosedeg.o sortd2.o 

CDET= cdetcomp.o zgetf2.o izamax.o dcabs1.o zscal.o \
	zswap.o zgeru.o xerbla.o

# library to compute stationary points of D2
libmoid.a: $(STATCOMP) $(CDET)
	ar r ../lib/libmoid.a $(STATCOMP) $(CDET)
	$(RANLIB) ../lib/libmoid.a

LIBDEPS= ../lib/libprop.a ../lib/libmoid.a ../lib/libgauss.a ../lib/libsuit.a

# test program

min_comp.x:min_comp.o  $(LIBDEPS)
	$(FC) $(FFLAGS) min_comp.o  $(LIBDEPS) \
	-o min_comp.x

stat_comp.x:stat_comp.o  $(LIBDEPS)
	$(FC) $(FFLAGS) stat_comp.o  $(LIBDEPS) \
	-o stat_comp.x

#stat_search.x: statsearch.o $(LIBDEPS) optpro.o 
#	$(FC) $(FFLAGS) statsearch.o $(LIBDEPS) optpro.o -o stat_search.x
stat_search.x: stat_search.o $(LIBDEPS) 
	$(FC) $(FFLAGS) stat_search.o $(LIBDEPS) -o stat_search.x


depend:
	../../bin/mkdep90 -I../include -L../suit *.f *.f90 *.h *.h90 > make.dep
	../../bin/mkmoddep -f ../.compiler ../suit > modules.flags

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


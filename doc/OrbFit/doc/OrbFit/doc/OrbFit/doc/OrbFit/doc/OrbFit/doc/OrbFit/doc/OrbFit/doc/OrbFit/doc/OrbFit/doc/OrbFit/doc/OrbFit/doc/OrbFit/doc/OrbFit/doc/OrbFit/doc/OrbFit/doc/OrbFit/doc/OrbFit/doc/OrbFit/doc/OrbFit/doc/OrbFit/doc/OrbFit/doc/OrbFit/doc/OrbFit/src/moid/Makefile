include ../make.flags
include modules.flags

all: libmoid.a  min_comp.x statcomp_ta_shift.x stat_search_ta_shift.x ast2ast_moid.x 

STATCOMP= compute_minima_ta.o \
	compute_critical_points_shift.o solvpoly_ta.o critical_points_shift.o \
	rvfft.o decode_out.o code_input.o irvfft.o choosedeg.o

CDET= cdetcomp.o zgetf2.o izamax.o dcabs1.o zscal.o \
	zswap.o zgeru.o xerbla.o

# library to compute stationary points of D2
libmoid.a: $(STATCOMP) $(CDET)
	ar r ../lib/libmoid.a $(STATCOMP) $(CDET)
	$(RANLIB) ../lib/libmoid.a

LIBDEPS= ../lib/libprop.a ../lib/libmoid.a ../lib/libgauss.a ../lib/libsuit.a

# programs

min_comp.x:mincomp_ta_shift.o  $(LIBDEPS)
	$(FC) $(FFLAGS) mincomp_ta_shift.o  $(LIBDEPS) \
	-o min_comp.x

statcomp_ta_shift.x:statcomp_ta_shift.o  $(LIBDEPS)
	$(FC) $(FFLAGS) statcomp_ta_shift.o  $(LIBDEPS) \
	-o statcomp_ta_shift.x

stat_search_ta_shift.x:stat_search_ta_shift.o $(LIBDEPS)
	$(FC) $(FFLAGS) stat_search_ta_shift.o $(LIBDEPS) \
	-o stat_search_ta_shift.x

ast2ast_moid.x:ast2ast_moid.o $(LIBDEPS)
	$(FC) $(FFLAGS) ast2ast_moid.o $(LIBDEPS) \
	-o ast2ast_moid.x

#statcomp3.x:statcomp3.o  $(LIBDEPS)
#	$(FC) $(FFLAGS) statcomp3.o  $(LIBDEPS) \
#	-o statcomp3.x

depend:
	../../bin/mkdep90 -I../include -L../suit -L../propag *.f *.f90 *.h *.h90 > make.dep
	../../bin/mkmoddep -f ../.compiler ../suit ../propag > modules.flags

clean:
	-(rm -f *.x *.o *.mod ../lib/libmoid.a *% *~ core work.pc*)

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

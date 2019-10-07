include ../make.flags
include modules.flags

all: ../lib/libprop.a 

PROP= propag_state.o force_model.o rmodel.o yark_pert.o moid_compute.o runge_kutta_gauss.o ra15_mod.o
      
CLOAP=  close_app.o tp_trace.o semi_linear.o
# target_plane.o close_app_data.o
PUB= ephem_prop.o plocbd.o fitsubs.o pred_obs.o

DIFCOR= multiple_sol.o least_squares.o obs_correl.o count_opp.o count_opp2.o



../lib/libprop.a: $(PROP) $(DIFCOR) $(PUB) $(CLOAP)
	ar r ../lib/libprop.a $(PROP) $(DIFCOR) $(PUB) $(CLOAP) 
	$(RANLIB) ../lib/libprop.a

depend:
	../../bin/mkdep90 -I../include -L../suit *.f *.f90 *.h *.h90 > make.dep
	../../bin/mkmoddep -f ../.compiler ../suit > modules.flags

clean:
	rm -f *.x *.o *.mod lib*.a *% *~ core work.pc*
	rm -f ../lib/libprop*.a

distclean: clean
	../lib/deinst_src

# Windows targets
win: ..\lib\prop.lib 

..\lib\prop.lib:$(PROP) $(DIFCOR) $(PUB) $(CLOAP)
	lib /out:$@ $**

winclean:
	@del *.o
	@del *.pdb
	@del ..\lib\prop.lib

# Dependencies
include make.dep

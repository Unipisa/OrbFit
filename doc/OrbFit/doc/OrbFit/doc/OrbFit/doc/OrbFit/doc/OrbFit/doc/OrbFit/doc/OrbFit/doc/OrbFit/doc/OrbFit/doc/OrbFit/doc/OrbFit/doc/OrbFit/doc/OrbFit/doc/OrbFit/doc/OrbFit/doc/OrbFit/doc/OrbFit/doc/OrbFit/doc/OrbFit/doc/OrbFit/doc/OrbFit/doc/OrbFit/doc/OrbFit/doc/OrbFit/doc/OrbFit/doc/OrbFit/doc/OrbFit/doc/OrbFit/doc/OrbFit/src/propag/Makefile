include ../make.flags
include modules.flags


all: ../lib/libprop.a 

PROP= propag_state.o force_model.o rmodel.o yark_pert.o \
        moid_compute.o runge_kutta_gauss.o ra15_mod.o
      
CLOAP=  close_app.o tp_trace.o semi_linear.o shower_analysistp.o \
        ret_analysistp.o offlov_checktp.o eval_risk.o obssto.o virtual_impactor.o

PUB= ephem_prop.o plocbd.o fitsubs.o pred_obs.o

DIFCOR= multiple_sol.o least_squares.o obs_correl.o count_opp.o count_opp2.o

QUADP = quadp/dqags.o quadp/dqagse.o quadp/dqagp.o quadp/dqagpe.o \
        quadp/dqelg.o quadp/dqk21.o quadp/dqpsrt.o quadp/dqagsc.o \
        quadp/dqagsec.o quadp/dqagpc.o quadp/dqagpec.o quadp/dqelgc.o \
        quadp/dqk21c.o quadp/dqpsrtc.o 

libquamia:
	cd quadp; make

../lib/libprop.a: $(PROP) $(DIFCOR) $(PUB) $(CLOAP) libquamia
	ar r ../lib/libprop.a $(PROP) $(DIFCOR) $(PUB) $(CLOAP) $(QUADP) 
	$(RANLIB) ../lib/libprop.a

depend:
	../../bin/mkdep90 -I../include -L../suit *.f *.f90 *.h *.h90 > make.dep
	../../bin/mkmoddep -f ../.compiler ../suit > modules.flags

clean:
	rm -f *.x *.o *.mod lib*.a *% *~ core work.pc*
	rm -f ../lib/libprop*.a
	cd quadp ; make clean  

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

.f.o:
	$(FC) $(FFLAGS) -c $*.f -o $*.o 

include ../make.flags

all: ../lib/libprop.a 

PROP= propag_state.o force_model.o moid_compute.o 
      
CLOAP= close_app.o target_plane.o 

PUB= predict_obs.o ephem_prop.o plocbd.o fitsubs.o

DIFCOR=obs_compute.o least_sq.o mult_sol.o 
# experimental zmulti.o

../lib/libprop.a: $(PROP) $(DIFCOR) $(PUB) $(CLOAP) 
	ar r ../lib/libprop.a $(PROP) $(DIFCOR) $(PUB) $(CLOAP)
	$(RANLIB) ../lib/libprop.a

depend:
	./mkdep *.f *.h ../include/*.h > make.dep

clean:
	-(rm -f *.x *.o ../lib/libprop*.a *% *~ core work.pc*); 

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

.f.o:
	$(FC) $(FFLAGS) -c $*.f -o $*.o 

include ../make.flags

all: ../lib/libgauss.a

# PUBLIC LIBRARY for Vaisala method
# (the only public routines are: vaisala, vaisub sel3g and alfdel)

VAIS=vaisala.o vaisub.o

LIB=iodini.o sel3g.o $(VAIS) iodet.o gaussn.o gibbs.o solv8.o gauseq.o \
    gauini.o ccek1.o modgeq.o anecc.o hypan.o ekcc1.o elekdt.o \
    sel3mc.o adds3p.o s3dtw.o s3srt.o sel3mg.o iodrms.o iodsdt.o iodsbs.o \
    random.o

../lib/libgauss.a: $(LIB)
	ar r ../lib/libgauss.a $(LIB)
	$(RANLIB) ../lib/libgauss.a

depend:
	./mkdep *.f *.h ../include/*.h > make.dep

clean:
	rm -f *.x *.o lib*.a *% *~ core work.pc*
	rm -f ../lib/libgauss.a

distclean: clean
	../lib/deinst_src

# Windows targets
win: ..\lib\gauss.lib

..\lib\gauss.lib: $(LIB)
	lib /out:$@ $**

winclean:
	@del *.o
	@del *.pdb
	@del ..\lib\gauss.lib

# Dependencies
include make.dep

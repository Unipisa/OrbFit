.f.o:
	$(FC) $(FFLAGS) -c $*.f -o $*.o 

include ../make.flags

all: ../lib/libsuit.a

# LIBRARY OBJECT FILES

# Orbital elements
ELEM_OBJ=orb_els.o
#	coocha.o prop2b.o carequ.o equcar.o kepequ.o equkep.o eqpequ.o \
#       equeqp.o ekensd.o rdelem.o wro1lh.o wro1lr.o outele.o outel1.o \
#       wromlh.o wromlr.o rdast1.o rdast2.o oefdet.o rdmpca.o mpcpds.o \
#       mpcdat.o rdoef.o fixcnm.o rdmpca2.o oporbf.o clorbf.o rdorb.o

# Observations & Statistics
OBS_OBJ=obs_wei.o
#        obsrms.o rdorms.o wrorms.o wrores.o fitwgt.o iaucod.o iaucod2.o \
#        wrirwo.o wrirms.o heapsort.o rearwo.o addobs.o inobs.o wrirwg.o \
#        astrow.o magrms.o rrmscl.o accstr.o crmscl.o pvobs2.o errel.o \
#        astrwb.o crmscn.o rdcorm.o obscor.o fcorob.o fcsfun.o \
#        mpcin.o mpctr.o mpcrad.o jplin.o jplrad.o sessag.o astrad.o rdanga.o

# Mathematical routines
MATH_OBJ= math_lib.o pol_zeros.o chebi_pol.o
#         heapsort.o princ.o prscal.o prvec.o vsize.o prodmv.o 
#         mulmat.o mulmat_cut.o lincom.o \
#         plegm.o plegd1.o evdis.o tchol.o mult3.o rotmt.o trsp.o \
#         matvet.o prodmm.o vetmat.o snorm.o snormd.o vdiff.o \
#         transp.o matin.o ident.o vcopy.o mcopy.o primea.o pridif.o intlo.o \
#         inthi.o bilin.o lincog.o vsumg.o prscag.o inv22.o mulmav.o tchinv.o \
#         covprs.o norprs.o polar.o roff.o truncat.o \
# these are used only in iers_ser.f
#         assmat.o pdmat.o pd1mat.o pd2mat.o prodvs.o rotmt1.o \
#         rotmt2.o summat.o sumv.o trsp3.o

#removed: inver2, invmat3, tchol2, epsini, comeps.h + all chebichev

# JPL Digital Ephemeris and binary ephemeris
EPH_OBJ=jpl_ephem.o 
#        rdbep.o jpllis.o fszer2.o dpleph.o state.o interp.o const.o split.o \
#        masjpl.o

# File operations
FILE_OBJ=file_oper.o
#	 filopn.o filclo.o filass.o filopl.o libini.o filopf.o dlifex.o \
#        filnam.o splinam.o fidinam.o dircom2.o

# Header namelists
HNAM_OBJ= header_input.o
#         rdfnam.o rdfcha.o rdfint.o rdflog.o rdfrea.o rdftim.o splkvc.o \
#         rdfref.o getrsc.o 

# Main option namelist
OPT_OBJ=option_input.o
#        getkv.o namini.o rdnam.o rdnam1.o rdnam2.o rdnam3.o rdklst.o \
#        chkkey.o rdncha.o rdnint.o rdnlog.o rdnrea.o rdntim.o rdnref.o \
#        chkfln.o rdmint.o rdmrea.o rdmcha.o rdvint.o rdvrea.o rdvcha.o \
#        rdkls1.o chkpdf.o

# Character strings
CHA_OBJ=char_str.o
#	lench.o rmsp.o norstr.o strcnt.o upcase.o stspli.o isnum.o islett.o \
#        nitchs.o locase.o spflds.o titast.o isbnum.o filstr.o

# Reference systems
REF_OBJ= ref_syst.o    
#       rotpn.o obleq.o rnut80.o prec.o nutn80.o equequ.o gmst.o  \
#       obscoo.o chkref.o ch2ref.o posobs.o pvobs.o

# Time and timescales
TIME_OBJ=time_scales.o
#	ch2tim.o intmon.o tjm1.o timnf.o itaiut.o chktsc.o cnvtim.o \
#         deltt.o mjddat.o chmon.o bessep.o chmo2i.o julian.o

# Handling of IERS time series
IERS_OBJ=iers_ser.o
#         delthp.o diurot.o dut1r.o dut1s.o equeqd.o gmsnom.o gmstd.o \
#         ierini.o iersts.o isbadr.o nutarg.o nutnd.o nutwhr.o obleqd.o \
#         precd.o rdbula.o rnutd.o rotpv.o rotsys.o xypol.o sv2int.o

# EISPACK routines (to be substituted?)
EIS_OBJ=eigen_val.o
#	rs.o tred1.o tqlrat.o tred2.o tql2.o

# Various utilities
UTIL_OBJ= util_suit.o 
#       angvcf.o  tee.o menu.o nights.o
# angles I/O routines
# rdanga.o sessag.o angvcf.o

# currently not used
NOT_USED= not_used.o
#	wrfcha.o wrfini.o wrfint.o wrfinv.o wrflog.o \
#       wrfrea.o wrfrev.o wrftim.o \
# 	eildri.o eildrf.o eildrl.o

SUIT_OBJ=$(FILE_OBJ) $(HNAM_OBJ) $(CHA_OBJ) $(UTIL_OBJ) $(TIME_OBJ) \
         $(ELEM_OBJ) $(MATH_OBJ) $(EPH_OBJ) $(REF_OBJ) \
         $(EIS_OBJ) $(OPT_OBJ) $(OBS_OBJ) $(IERS_OBJ)

../lib/libsuit.a: $(SUIT_OBJ)
	ar r ../lib/libsuit.a $(SUIT_OBJ)
	$(RANLIB) ../lib/libsuit.a

depend:
	./mkdep *.f *.h ../include/*.h > make.dep

clean:
	-(rm -f *.x *.o ../lib/libsuit.a *% *~ core work.pc*); 

distclean: clean
	../lib/deinst_src

# Windows targets
win: ..\lib\suit.lib

..\lib\suit.lib:$(SUIT_OBJ)
	lib /out:$@ $**

winclean:
	@del *.o
	@del *.pdb
	@del ..\lib\suit.lib

# Dependencies
include make.dep

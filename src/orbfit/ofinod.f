* Copyright (C) 1997-1999 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: February 11, 1999
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         O F I N O D                           *
*  *                                                               *
*  *                Auxiliary routine for ORBFIT:                  *
*  *                 initial orbit determination                   *
*  *                                                               *
*  *****************************************************************
**
* INPUT:    UNIREP    -  FORTRAN unit for report
*           ELEOUT    -  Name of file for output of orbital elements
*           OPT       -  Initial orbit determination option
*           NAME      -  Object names
*           NAMOF     -  Names of observation files (w/o extension)
*           DEFORB    -  Orbit definition flag
*           DEFCN     -  Tells whether covariance/normal matrices
*                            are defined
*           H,G       -  Magnitude parameters
*           DIR       -  Directory containing observation/residual file
*           NOBJ      -  Number of objects
*           OBJID     -  Object identifier
*           TDT       -  Times of observations (MJD, TDT)
*           TUTM      -  Time (MJD, UTM)
*           ALPHA     -  Right ascension (rad)
*           DELTA     -  Declination (rad)
*           IOBS      -  Observation type
*           OBSCOD    -  Observatory code
*           SMAG      -  Measured magnitude (string)
*           RMSA      -  A-priori RMS of RA (rad)
*           RMSD      -  A-priori RMS of DEC (rad)
*           RMSMAG    -  A-priori RMS of magnitude
*           SEL       -  Selection index (0=don't use; 1=use for fit;
*                                         2=use for fit & Gauss method)
*           N         -  Number of observations for each object
*           NT        -  Total number of observations
*           IP1       -  Pointer to first observation for each object
*
* OUTPUT:   ELEM      -  Orbital elements
*           TELEM     -  Epoch of orbital elements (MJD, TDT)
*           ELTYPE    -  Type of orbital elements (EQU/KEP/CAR)
*           COMELE    -  Comment on orbital elements
*           RESA      -  Residuals in RA (rad)
*           RESD      -  Residuals in DEC (rad)
*
* Input variables DEFORB, DEFCN, SEL are modified by the routine
*
      SUBROUTINE ofinod(unirep,eleout,opt,name,namof,deforb,defcn,h,g,
     +                  dir,nobj,objid,tdt,tutm,alpha,delta,iobs,
     +                  obscod,smag,rmsa,rmsd,rmsmag,sel,n,nt,ip1,elem,
     +                  telem,eltype,comele,resa,resd)
      IMPLICIT NONE

      INTEGER unirep,opt,nobj,n(nobj),nt,ip1(nobj),sel(nt),obscod(nt)
      INTEGER iobs(nt)
      LOGICAL deforb(nobj),defcn(nobj)
      DOUBLE PRECISION tdt(nt),tutm(nt),alpha(nt),delta(nt)
      DOUBLE PRECISION rmsa(nt),rmsd(nt),rmsmag(nt),resa(nt),resd(nt)
      DOUBLE PRECISION elem(6,nobj),telem(nobj),h(nobj),g(nobj)
      CHARACTER*(*) eleout,name(nobj),namof(nobj),eltype(nobj),dir(nobj)
      CHARACTER*(*) comele(nobj),objid(nt),smag(nt)

      INTEGER i,ln,uniele,lc,ld
      DOUBLE PRECISION cove(6,6),nore(6,6)
      CHARACTER*150 rwofil
      LOGICAL opnd,doit,fail
      INCLUDE 'parcmc.h'

      INTEGER lench
      EXTERNAL lench

      opnd=.false.

      DO 1 i=1,nobj
      doit=.false.
      IF(opt.EQ.1) THEN
          doit=(.NOT.deforb(i))
      ELSEIF(opt.EQ.2) THEN
          doit=.true.
      END IF
      IF(.NOT.doit) GOTO 1
      ln=lench(namof(i))
      ld=lench(dir(i))
      IF(ld.GT.0) THEN
          rwofil=dir(i)(1:ld)//namof(i)(1:ln)//'.rwo'
      ELSE
          rwofil=namof(i)(1:ln)//'.rwo'
      END IF
      CALL iodet(unirep,rwofil,name(i),objid(ip1(i)),tdt(ip1(i)),
     +           tutm(ip1(i)),alpha(ip1(i)),delta(ip1(i)),
     +           iobs(ip1(i)),obscod(ip1(i)),smag(ip1(i)),rmsa(ip1(i)),
     +           rmsd(ip1(i)),rmsmag(ip1(i)),sel(ip1(i)),n(i),2,
     +           elem(1,i),telem(i),resa(ip1(i)),resd(ip1(i)),eltype(i),
     +           comele(i),fail)
      deforb(i)=(.NOT.fail)
      defcn(i)=.false.
      IF(.NOT.fail) THEN
          IF(.NOT.opnd) THEN
              CALL filopn(uniele,eleout,'UNKNOWN')
              CALL wromlh(uniele,'ECLM','J2000')
              opnd=.true.
          END IF
          lc=lench(comele(i))
          WRITE(uniele,300) comcha,comele(i)(1:lc)
          CALL wromlr(uniele,name(i),elem(1,i),eltype(i),telem(i),cove,
     +                .false.,nore,.false.,h(i),g(i),0.D0)
      END IF
 300  FORMAT(A,1X,A)

 1    CONTINUE

      IF(opnd) CALL filclo(uniele,' ')

      END

* Copyright (C) 1997-1999 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: February 11, 1999
* Version: November 4, 1999 (new inobs call, SRC)
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         O F I O B S                           *
*  *                                                               *
*  *                Auxiliary routine for ORBFIT:                  *
*  *                    input of observations                      *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    UNIREP    -  FORTRAN unit for report
*           NAME      -  Object file names
*           NAMOF     -  Names of observation files (w/o extension)
*           DIR       -  Directory containing observation/residual files
*           NOBJ      -  Number of objects
*           NOBX      -  Max number of observations
*           OPUPD     -  Option: update weight file
*
* OUTPUT:   N         -  Number of observations for each object
*           NT        -  Total number of observations
*           IP1       -  Pointer to first observation for each object
*           IP2       -  Pointer to first weight for each object
*           OBJID     -  Object identifier
*           TDT       -  Time (MJD, TDT)
*           TUTM      -  Time (MJD, UTM)
*           ALPHA     -  Right ascension (rad)
*           DELTA     -  Declination (rad)
*           IOBS      -  Observation type
*           SMAG      -  Measured magnitude (string)
*           OBSCOD    -  Observatory code


*           SEL       -  Selection index (0=don't use; 1=use for fit;
*                        2=use for fit & Gauss method)
*           RMSA      -  A-priori RMS of RA (rad)
*           RMSD      -  A-priori RMS of DEC (rad)
*           RMSMAG    -  A-priori RMS of magnitude
*
      SUBROUTINE ofiobs(unirep,name,namof,dir,nobj,nobx,opupd,n,nt,
     +                  ip1,ip2,objid,tdt,tutm,alpha,delta,iobs,smag,
     +                  obscod,sel,rmsa,rmsd,rmsmag)
      IMPLICIT NONE

      INTEGER unirep,nobj,nobx,n(nobj),nt,ip1(nobj),ip2(nobj)
      INTEGER obscod(nobx),sel(nobx),iobs(nobx),opupd
      CHARACTER*(*) name(nobj),namof(nobj),dir(nobj)
      CHARACTER*(*) objid(nobx),smag(nobx)
      DOUBLE PRECISION tdt(nobx),tutm(nobx),alpha(nobx),delta(nobx)
      DOUBLE PRECISION rmsa(nobx),rmsd(nobx),rmsmag(nobx)
      LOGICAL precob,ok

      INTEGER i,ln,ld
      LOGICAL change

      INTEGER lench
      EXTERNAL lench

      WRITE(unirep,120)
 120  FORMAT('Input of observations:')

      nt=0
      precob=.false.

      DO 1 i=1,nobj
      ip1(i)=nt+1
      ip2(i)=2*nt+1
      CALL inobs(dir(i),namof(i),precob,objid(ip1(i)),ok,n(i),
     +           iobs(ip1(i)),tdt(ip1(i)),alpha(ip1(i)),delta(ip1(i)),
     +           tutm(ip1(i)),obscod(ip1(i)),sel(ip1(i)),rmsa(ip1(i)),
     +           rmsd(ip1(i)),rmsmag(ip1(i)),smag(ip1(i)),nobx-nt,
     +           unirep,change)
      IF(.NOT.ok) STOP '**** ofiobs: abnormal end ****'
      nt=nt+n(i)
      ld=lench(dir(i))
      ln=lench(name(i))
      WRITE(unirep,200) n(i),name(i)(1:ln)
 200  FORMAT(I9,' observations read for object ',A)
      CALL srtobs(objid(ip1(i)),tdt(ip1(i)),tutm(ip1(i)),alpha(ip1(i)),
     +            delta(ip1(i)),iobs(ip1(i)),smag(ip1(i)),
     +            obscod(ip1(i)),sel(ip1(i)),rmsa(ip1(i)),
     +            rmsd(ip1(i)),rmsmag(ip1(i)),n(i))
 1    CONTINUE

      END

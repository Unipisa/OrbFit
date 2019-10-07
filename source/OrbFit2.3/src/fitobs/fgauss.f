* Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)
* Version OFINOD: December 1, 1997 
* Adapted for FITOBS, AM/ZK March 1998
* Adapted for use of VAISALA, ZK November 1998 
* IOBS added (Feb 10, 1999) MC
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         F G A U S S                           *
*  *                                                               *
*  *                Auxiliary routine for FITOBS:                  *
*  *                 initial orbit determination                   *
*  *                                                               *
*  *****************************************************************
*
* In the present version, only Gauss' method is supported
*
* INPUT:    UNIREP    -  FORTRAN unit for report
*           UNIELE    -  FORTRAN unit for elements
*           NAME      -  Object name
*           DEFORB    -  Orbit definition flag
*           DEFCOV    -  Orbit covariance definition flag
*           H,G       -  Magnitude parameters
*           RWFIL     -  Names of residual/weight files
*           TDT       -  Times of observations (MJD, TDT)
*           TUTM      -  Time (MJD, UTM)
*           ALPHA     -  Right ascension (rad)
*           DELTA     -  Declination (rad)
*           IOBS      -  Observation type (1=RA,DEC; 2=R,RDOT)
*           OBSCOD    -  Observatory code
*           RMSA      -  A-priori RMS of RA (rad)
*           RMSD      -  A-priori RMS of DEC (rad)
*           OBJID     -  identifier for each observation
*           SMAG      -  magnitude observation (string)
*           RMSMAG    -  magnitude weight
*           SEL       -  Selection index (0=don't use; 1=use for fit;
*                        2=use for fit & Gauss method)
*           N         -  Number of observations for each object
*           IMETH     -  Method to be used: 1=Gauss, 2=Vaisala 
*
* OUTPUT:   EQ        -  Equinoctal Orbital elements
*           TELEM     -  Epoch of orbital elements (MJD, TDT)
* NOT OUTPUT, but available:
*           RESA      -  Residuals in RA (rad)
*           RESD      -  Residuals in DEC (rad)
*
* Input variables DEFORB, DEFCOV, SEL are modified by the routine
*
      SUBROUTINE fgauss(unirep,uniele,name,deforb,defcov,h,g,
     +     rwfil,tdt,tutm,alpha,delta,iobs,obscod,
     +     rmsa,rmsd,objid,smag,rmsmag,sel,n,imeth,
     +     eq,telem)
      IMPLICIT NONE

      INTEGER unirep,uniele,n,sel(n),obscod(n),imeth,iobs(n)
      LOGICAL deforb,defcov
      DOUBLE PRECISION tdt(n),tutm(n),alpha(n),delta(n)
      DOUBLE PRECISION rmsa(n),rmsd(n)
      INCLUDE 'parobx.h'
c asteroid identifier, apparent magnitude
      CHARACTER*9 objid(nobx)
      CHARACTER*6 smag(nobx)
      DOUBLE PRECISION rmsmag(nobx)
      DOUBLE PRECISION rms
      DOUBLE PRECISION elem(6),eq(6),telem,h,g,gg(6,6),cc(6,6)
      CHARACTER*18 name
      CHARACTER*60 rwfil
      INCLUDE 'trig.h'
      INTEGER ln,j

      INCLUDE 'parcmc.h'

      INTEGER lench
      EXTERNAL lench
c ======== constant of gravitation ==============
      INCLUDE 'sunmass.h'
c =====initial orbit determination===========================
      INCLUDE 'pariod.h'
      INCLUDE 'comiod.h'
      LOGICAL fail
      CHARACTER*3 eltype
      CHARACTER*60 comele
      DOUBLE PRECISION resa(nobx),resd(nobx),enne
c added for RMS
      DOUBLE PRECISION csi(nob2x),w(nob2x),snormd
      integer nused
      EXTERNAL snormd

      ln=lench(name)
      IF(imeth.eq.1)THEN
         WRITE(unirep,203) name(1:ln)
 203     FORMAT('Automatic method for object ',A,':')
         iodnm=2
         iodmet(1)=1
         iodmen(1)='GAUSS'
         iodmet(2)=2
         iodmen(2)='VAISALA'
      ELSEIF(imeth.eq.2)THEN
         WRITE(unirep,201) name(1:ln)
 201     FORMAT('Gauss method for object ',A,':')
         iodnm=1
         iodmet(1)=1
         iodmen(1)='GAUSS'
      ELSEIF(imeth.eq.3)THEN
         WRITE(unirep,202) name(1:ln)
 202     FORMAT('Vaisala method for object ',A,':')
         iodnm=1
         iodmet(1)=2
         iodmen(1)='VAISALA'
      ENDIF
      CALL iodet(unirep,rwfil,name,objid,tdt,tutm,alpha,delta,
     +     iobs,obscod,smag,rmsa,rmsd,rmsmag,sel,n,1,elem,
     +     telem,resa,resd,eltype,comele,fail)
c *************************************************************************

      deforb=(.NOT.fail)
      defcov=.false.
      IF(fail) THEN
         WRITE(*,204)
 204     FORMAT(2X,'INITIAL ORBIT DETERMINATION FAILED')
      ELSE
c output new elements, multiline format, no header, in .fel file.
         CALL wro1lh(uniele,'ECLM','J2000','KEP')
         CALL wromlr(uniele,name,elem,'KEP',telem,gg,.false.,cc,.false.
     +        ,h,g,0.d0)
c compute residuals RMS 
         call fitwgt(rmsa,rmsd,delta,sel,iobs,w,n,.TRUE.)
         do j=1,n
            csi(2*j-1)=resa(j)
            csi(2*j)=resd(j)
         enddo
         rms = snormd(csi,w,2*n,nused)
         WRITE(unirep,222)rms
         WRITE(*,222)rms
 222     FORMAT(' RMS of residuals, preliminary orbit=',f10.4)
c change to equinoctal
         CALL coocha(elem,eltype,gms,eq,'EQU',enne)
         WRITE(*,220) telem,eq
 220     FORMAT(' preliminary orbit elements for epoch=',f12.4/6f13.7)
      END IF
      RETURN
      END





* Copyright (C) 2000 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: June 5, 2000
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         I O D R M S                           *
*  *                                                               *
*  *      Orbital residuals and RMS for preliminary orbits         *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    ELEM      -  Orbital elements (mean ecliptic
*                                          and equinox J2000)
*           ELTYPE    -  Type of orbital elements
*           TELEM     -  Epoch of orbital elements (MJD, TDT)
*           TOBS      -  Time of observations (MJD, TDT)
*           ALPHA     -  Right ascension (rad)
*           DELTA     -  Declination (rad)
*           IOBS      -  Observation type (1=RA,DEC; 2=R,RDOT)
*           OBSCOD    -  Observatory code
*           SEL       -  Selection index (0: don't use; >=1: use)
*           IO1,IO2   -  First and last observations to be used
*
* OUTPUT:   RESA      -  Residuals in RA (rad)
*           RESD      -  Residuals in DEC (rad)
*           RMS       -  RMS error of residuals (rad)
*
      SUBROUTINE iodrms(elem,eltype,telem,tobs,alpha,delta,iobs,obscod,
     +                  sel,io1,io2,resa,resd,rms)
      IMPLICIT NONE

      INTEGER io1,io2
      INTEGER iobs(io2),obscod(io2),sel(io2)
      DOUBLE PRECISION elem(6),telem,tobs(io2),alpha(io2),delta(io2)
      DOUBLE PRECISION resa(io2),resd(io2),rms
      CHARACTER*(*) eltype

      INCLUDE 'trig.h'

      INTEGER i,k,ns
      DOUBLE PRECISION rot(3,3),gk,gkp,gm,xobs(3),xast(3),vast(3),xv(6)
      DOUBLE PRECISION xrel(3),alphac,deltac,rrel
      LOGICAL first
      DATA first/.true./

      SAVE first,rot,gm

      DOUBLE PRECISION princ
      EXTERNAL princ

      IF(first) THEN
          gk=0.01720209895d0
          gkp=gk
          gm=gkp**2
          CALL rotpn(rot,'ECLM','J2000',0.d0,'EQUM','J2000',0.d0)
          first=.false.
      END IF

      rms=0.d0
      ns=0
      DO 1 k=io1,io2
      IF(iobs(k)/1000.NE.1) GOTO 1
      IF(sel(k).LT.1) GOTO 1
      CALL posobs(tobs(k),obscod(k),1,xobs)
      CALL ekcc1(elem,eltype,xv,gm,tobs(k)-telem)
      CALL prodmv(xast,rot,xv(1))
      CALL prodmv(vast,rot,xv(4))
      DO 2 i=1,3
      xrel(i)=xast(i)-xobs(i)
 2    CONTINUE
      CALL aber1(xrel,vast,xrel)
      CALL polar(xrel,alphac,deltac,rrel)
      resa(k)=princ(alpha(k)-alphac)
      IF(resa(k).GT.pig) resa(k)=resa(k)-dpig
      resd(k)=delta(k)-deltac
      rms=rms+(COS(delta(k))*resa(k))**2+resd(k)**2
      ns=ns+1
 1    CONTINUE

      IF(ns.GT.0) rms=SQRT(rms/(2*ns))

      END

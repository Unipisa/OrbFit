* Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: December 1, 1997
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         O F E P H E                           *
*  *                                                               *
*  *                Auxiliary routine for ORBFIT:                  *
*  *                          ephemerides                          *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    UNIT      -  FORTRAN unit for output
*           NAME      -  Object names
*           DEFORB    -  Orbit definition flag
*           DEFCOV    -  Orbit covariance definition flag
*           ELEM      -  Orbital elements
*           TELEM     -  Epoch of orbital elements (MJD, TDT)
*           ELTYPE    -  Type of orbital elements (EQU/KEP/CAR)
*           COVE      -  Covariance matrix of orbital elements
*           H,G       -  Magnitude parameters
*           MASS      -  Object masses
*           COMELE    -  Comment on orbital elements
*           NOBJ      -  Number of objects
*
      SUBROUTINE ofephe(unit,name,deforb,defcov,elem,telem,eltype,cove,
     +                  h,g,mass,comele,nobj)
      IMPLICIT NONE

      INTEGER unit,nobj
      DOUBLE PRECISION elem(6,nobj),telem(nobj),mass(nobj)
      DOUBLE PRECISION cove(6,6,nobj),h(nobj),g(nobj)
      LOGICAL deforb(nobj),defcov(nobj)
      CHARACTER*(*) name(nobj),eltype(nobj),comele(nobj)

* NEEDED common blocks:
      INCLUDE 'comeph.h'

      INTEGER i,k,ln,lc

      INTEGER lench
      EXTERNAL lench

      IF(iiceph.NE.36) STOP '**** ofephe: internal error (01) ****'

      DO 1 k=1,nepobj
      i=kepobj(k)
      IF(i.LT.1 .OR. i.GT.nobj) GOTO 1
      IF(.NOT.deforb(i)) GOTO 1
      ln=lench(name(i))
      WRITE(unit,100) name(i)(1:ln)
 100  FORMAT(/'Ephemeris for object ',A)
      lc=lench(comele(i))
      IF(lc.GT.0) WRITE(unit,101) comele(i)(1:lc)
 101  FORMAT(5X,'Origin of orbital elements: ',A)
      WRITE(unit,102)
 102  FORMAT(/)
      CALL ephemc(unit,eltype(i),telem(i),elem(1,i),cove(1,1,i),
     +            defcov(i),teph1,teph2,dteph,mass(i),
     +            h(i),g(i),idsta,ephtsc,ephfld)
 1    CONTINUE

      END

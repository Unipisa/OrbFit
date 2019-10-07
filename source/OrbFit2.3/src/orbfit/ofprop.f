* Copyright (C) 1997-2000 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: June 7, 2000
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         O F P R O P                           *
*  *                                                               *
*  *                Auxiliary routine for ORBFIT:                  *
*  *               propagation of orbital elements                 *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    UNIREP    -  FORTRAN unit for report
*           UNIELE    -  FORTRAN unit for output of orbital elements
*           OPELE     -  Flag (need to open UNIELE)
*           ELEOUT    -  Name of file for output of orbital elements
*           NAME      -  Object names
*           DEFORB    -  Orbit definition flag
*           DEFCN     -  Tells whether covariance/normal matrices
*                            are defined
*           ELEM      -  Orbital elements
*           TELEM     -  Epoch of orbital elements (MJD, TDT)
*           ELTYPE    -  Type of orbital elements (EQU/KEP/CAR)
*           COVE      -  Covariance matrix of orbital elements
*           NORE      -  Normal matrix of orbital elements
*           H,G       -  Magnitude parameters
*           MASS      -  Object masses
*           COMELE    -  Comment on orbital elements
*           NOBJ      -  Number of objects
*           GMSUN     -  G*Mass(Sun)
*           OEPTIM    -  Epoch of output elements (MJD, TDT)
*           OEPSET    -  Flag stating that an output epoch is requested
*           OETYPE    -  Type of output elements (CAR/EQU/KEP/EQP)
*
      SUBROUTINE ofprop(unirep,uniele,opele,eleout,name,deforb,defcn,
     +                  elem,telem,eltype,cove,nore,h,g,mass,
     +                  comele,nobj,gmsun,oeptim,oepset,oetype)
      IMPLICIT NONE

      INTEGER unirep,uniele,nobj
      DOUBLE PRECISION elem(6,nobj),telem(nobj),gmsun
      DOUBLE PRECISION cove(6,6,nobj),nore(6,6,nobj)
      DOUBLE PRECISION h(nobj),g(nobj),mass(nobj),oeptim
      LOGICAL opele,deforb(nobj),defcn(nobj),oepset
      CHARACTER*(*) eleout,name(nobj),eltype(nobj),comele(nobj),oetype

      INCLUDE 'parcmc.h'

      INTEGER i,ln,lc
      DOUBLE PRECISION gma1,enne,dxde(6,6)
      DOUBLE PRECISION elemt(6),covet(6,6),noret(6,6)
      DOUBLE PRECISION elemo(6),coveo(6,6),noreo(6,6)
      LOGICAL error,defnro

      INTEGER lench
      EXTERNAL lench

      IF(.NOT.oepset) RETURN
      WRITE(unirep,100)
 100  FORMAT('Propagation of orbital elements:')

      DO 1 i=1,nobj
      IF(.NOT.deforb(i)) GOTO 1
      gma1=gmsun*(1+mass(i))

      IF(defcn(i)) THEN
          CALL proelc(eltype(i),telem(i),elem(1,i),cove(1,1,i),
     +                nore(1,1,i),oeptim,elemt,covet,noret)
          CALL cooder(elemt,eltype(i),gma1,elemo,oetype,enne,dxde)
          CALL covprs(covet,dxde,6,coveo)
          CALL norprs(noret,dxde,6,noreo,error)
          defnro=(.NOT.error)
      ELSE
          CALL proele(eltype(i),telem(i),elem(1,i),oeptim,elemt)
          CALL coocha(elemt,eltype(i),gma1,elemo,oetype,enne)

      END IF
      ln=lench(name(i))
      lc=lench(comele(i))
      IF(lc.LE.0) THEN
          WRITE(unirep,101) name(i)(1:ln)
      ELSE
          WRITE(unirep,102) name(i)(1:ln),comele(i)(1:lc)
      END IF
      CALL outele(unirep,elemo,oetype,oeptim,' ',.true.,.false.)

      IF(opele) THEN
          OPEN(uniele,FILE=eleout,STATUS='UNKNOWN')
          opele=.false.
          WRITE(uniele,301) comcha
          CALL wromlh(uniele,'ECLM','J2000')
      END IF
      IF(lc.LE.0) THEN
          WRITE(uniele,103) comcha,name(i)(1:ln)
      ELSE
          WRITE(uniele,104) comcha,name(i)(1:ln),comele(i)(1:lc)
      END IF
      CALL wromlr(uniele,name(i),elemo,oetype,oeptim,
     +            coveo,defcn(i),noreo,defnro,h(i),g(i),mass(i))
 1    CONTINUE
 101  FORMAT(5X,'Propagated orbital elements for object ',A,':')
 102  FORMAT(5X,'Propagated orbital elements for object ',A,' (',A,'):')
 103  FORMAT(A,' Orbital elements of object ',A)
 104  FORMAT(A,' Orbital elements of object ',A,' (',A,')')
 301  FORMAT(A,' Propagated orbits')

      END

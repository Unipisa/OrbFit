* Copyright (C) 1997-2000 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: June 7, 2000
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         O F I O R B                           *
*  *                                                               *
*  *                Auxiliary routine for ORBFIT:                  *
*  *                  input of orbital elements                    *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    UNIREP    -  FORTRAN unit for report
*           ELFT      -  Input file names (for all objects)
*           NELFT     -  Number of input files (for all objects)
*           ELF1      -  Input file names (for each object)
*           NELF1     -  Number of input files (for each object)
*           NAME      -  Object names (true)
*           NAMEO     -  Names to be searched in orbital element files
*           NOBJ      -  Number of objects
*           NOBJX     -  Max number of objects
*           NFIX      -  First physical dimension of array ELF1
*
* OUTPUT:   DEFORB    -  Orbit definition flag
*           DEFCN     -  Tells whether covariance/normal matrices
*                            are defined
*           ELTYPE    -  Type of orbital elements (EQU/KEP/CAR)
*           ELEM      -  Orbital elements (ECLM J2000)
*           TELEM     -  Epoch of orbital elements (MJD, TDT)
*           COVE      -  Covariance matrix of orbital elements
*           NORE      -  Normal matrix of orbital elements
*           H,G       -  Magnitude parameters
*           MASS      -  Mass
*           COMELE    -  Comment on orbital elements
*
      SUBROUTINE ofiorb(unirep,elft,nelft,elf1,nelf1,name,nameo,
     +                  nobj,nobjx,deforb,defcn,eltype,elem,
     +                  telem,cove,nore,h,g,mass,comele,nfix)
      IMPLICIT NONE

      INTEGER unirep,nobj,nobjx,nelft,nelf1(nobj),nfix
      LOGICAL deforb(nobjx),defcn(nobjx)
      DOUBLE PRECISION h(nobjx),g(nobjx),mass(nobjx),elem(6,nobj)
      DOUBLE PRECISION telem(nobj),cove(6,6,nobj),nore(6,6,nobj)
      CHARACTER*(*) elft(nfix),elf1(nfix,nobj)
      CHARACTER*(*) name(nobj),nameo(nobj),eltype(nobj),comele(nobj)

      INTEGER i,ln,lc,lno

      INTEGER lench
      EXTERNAL lench

      WRITE(unirep,122)
 122  FORMAT('Input of orbital elements:')

* No orbits defined yet
      DO 1 i=1,nobj
      deforb(i)=.false.
      defcn(i)=.false.
      h(i)=-1.D9
      g(i)=0.d0
      mass(i)=0.d0
 1    CONTINUE

* Input for particular files for each object
      DO 2 i=1,nobj
      CALL rdelem(unirep,nameo(i),1,elf1(1,i),nelf1(i),deforb(i),
     +            defcn(i),eltype(i),telem(i),elem(1,i),
     +            cove(1,1,i),nore(1,1,i),mass(i),h(i),g(i),comele(i))
 2    CONTINUE

* Input from common files
      CALL rdelem(unirep,nameo,nobj,elft,nelft,deforb,defcn,
     +            eltype,telem,elem,cove,nore,mass,h,g,comele)

      DO 3 i=1,nobj
      IF(deforb(i)) THEN
          ln=lench(name(i))
          IF(name(i).EQ.nameo(i)) THEN
              WRITE(unirep,123) name(i)(1:ln)
          ELSE
              lno=lench(nameo(i))
              WRITE(unirep,125) name(i)(1:ln),nameo(i)(1:lno)
          END IF
          lc=lench(comele(i))
          IF(lc.GT.0) WRITE(unirep,124) comele(i)(1:lc)
          CALL outele(unirep,elem(1,i),eltype(i),telem(i),' ',
     +                .true.,.false.)
      END IF
 3    CONTINUE
 123  FORMAT(5X,'Orbital elements for object ',A,':')
 125  FORMAT(5X,'Orbital elements for object ',A,' (with name ',A,'):')
 124  FORMAT(5X,'(origin: ',A,')')

      END

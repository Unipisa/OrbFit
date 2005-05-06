! Copyright (C) 1997-2000 by Mario Carpino (carpino@brera.mi.astro.it)
! Version: June 7, 2000
! ---------------------------------------------------------------------
!
!  *****************************************************************
!  *                                                               *
!  *                         O F I O R B                           *
!  *                                                               *
!  *                Auxiliary routine for ORBFIT:                  *
!  *                  input of orbital elements                    *
!  *                                                               *
!  *****************************************************************
!
! INPUT:    UNIREP    -  FORTRAN unit for report
!           ELFT      -  Input file names (for all objects)
!           NELFT     -  Number of input files (for all objects)
!           ELF1      -  Input file names (for each object)
!           NELF1     -  Number of input files (for each object)
!           NAME      -  Object names (true)
!           NAMEO     -  Names to be searched in orbital element files
!           NOBJ      -  Number of objects
!           NOBJX     -  Max number of objects
!           NFIX      -  First physical dimension of array ELF1
!
! OUTPUT:   ELEM      -  Orbital elements
!           ELEM_UNC  -  Orbital element uncertainty
!           DEFORB    -  Orbit definition flag
!           DEFCN     -  Tells whether covariance/normal matrices are defined
!           MASS      -  Mass
!           COMELE    -  Comment on orbital elements
!
      SUBROUTINE ofiorb(unirep,elft,nelft,elf1,nelf1,name,nameo,nobj,elem,elem_unc,deforb,defcn,mass,comele,nfix)

      USE orbit_elements
      IMPLICIT NONE

      INCLUDE 'parnob.h90'

      INTEGER unirep,nobj,nelft,nelf1(nobj),nfix
      TYPE(orbit_elem),DIMENSION(nobjx) :: elem
      TYPE(orb_uncert),DIMENSION(nobjx) :: elem_unc
      LOGICAL deforb(nobjx),defcn(nobjx)
      DOUBLE PRECISION mass(nobjx)
      CHARACTER*(*) elft(nfix),elf1(nfix,nobj)
      CHARACTER*(*) name(nobj),nameo(nobj),comele(nobj)

      DOUBLE PRECISION telem(nobj1x),elemv(6,nobj1x)
      DOUBLE PRECISION h(nobj1x),g(nobj1x)
      DOUBLE PRECISION cove(6,6,nobj1x),nore(6,6,nobj1x)
      CHARACTER*3 eltype(nobj1x)

      INTEGER i,ln,lc,lno

      INTEGER lench
      EXTERNAL lench

      WRITE(unirep,122)
  122 FORMAT('Input of orbital elements:')

! No orbits defined yet
      DO 1 i=1,nobj
      deforb(i)=.false.
      defcn(i)=.false.
      mass(i)=0.d0
    1 END DO

! Input for particular files for each object
      DO 2 i=1,nobj
      CALL rdelem(unirep,nameo(i),1,elf1(1,i),nelf1(i),deforb(i),       &
     &            defcn(i),eltype(i),telem(i),elemv(1,i),               &
     &            cove(1,1,i),nore(1,1,i),mass(i),h(i),g(i),comele(i))
    2 END DO

! Input from common files
      CALL rdelem(unirep,nameo,nobj,elft,nelft,deforb,defcn,eltype,telem,elemv,cove,nore,mass,h,g,comele)

      DO i=1,nobj
         elem(i)=undefined_orbit_elem
         IF(deforb(i)) THEN
            elem(i)%coo=eltype(i)
            elem(i)%ndim=6
            elem(i)%coord=elemv(1:6,i)
            elem(i)%t=telem(i)
            elem(i)%mag_set=(h(i) > 0.d0)
            elem(i)%h_mag=h(i)
            elem(i)%g_mag=g(i)
            ln=lench(name(i))
            IF(name(i).EQ.nameo(i)) THEN
               WRITE(unirep,123) name(i)(1:ln)
            ELSE
               lno=lench(nameo(i))
               WRITE(unirep,125) name(i)(1:ln),nameo(i)(1:lno)
            END IF
            lc=lench(comele(i))
            IF(lc.GT.0) WRITE(unirep,124) comele(i)(1:lc)
            CALL outele(unirep,elemv(1,i),eltype(i),telem(i),' ',.true.,.false.)
         END IF
         elem_unc(i)=undefined_orb_uncert
         IF(defcn(i)) THEN
            elem_unc(i)%c=nore(1:6,1:6,i)
            elem_unc(i)%g=cove(1:6,1:6,i)
            elem_unc(i)%succ=.true.
            elem_unc(i)%ndim=6
         END IF
      END DO
  123 FORMAT(5X,'Orbital elements for object ',A,':')
  125 FORMAT(5X,'Orbital elements for object ',A,' (with name ',A,'):')
  124 FORMAT(5X,'(origin: ',A,')')

      END

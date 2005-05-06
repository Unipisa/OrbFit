!
!  *****************************************************************
!  *                                                               *
!  *                         O F E P H E                           *
!  *                                                               *
!  *                Auxiliary routine for ORBFIT:                  *
!  *                          ephemerides                          *
!  *                                                               *
!  *****************************************************************
!
! INPUT:    UNIT      -  FORTRAN unit for output
!           NAME      -  Object names
!           DEFORB    -  Orbit definition flag
!           DEFCOV    -  Orbit covariance definition flag
!           ELEM      -  Orbital elements
!           ELEM_UNC  -  Orbital element uncertainty
!           MASS      -  Object masses
!           COMELE    -  Comment on orbital elements
!           NOBJ      -  Number of objects
!
      SUBROUTINE ofephe(unit,name,deforb,defcov,elem,elem_unc,mass,comele,nobj)
      USE orbit_elements
      IMPLICIT NONE

      INTEGER unit,nobj
      TYPE(orbit_elem),DIMENSION(nobj) :: elem
      TYPE(orb_uncert),DIMENSION(nobj) :: elem_unc
      DOUBLE PRECISION mass(nobj)
      LOGICAL deforb(nobj),defcov(nobj)
      CHARACTER*(*) name(nobj),comele(nobj)

! NEEDED common blocks:
      INCLUDE 'comeph.h90'

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
  100 FORMAT(/'Ephemeris for object ',A)
      lc=lench(comele(i))
      IF(lc.GT.0) WRITE(unit,101) comele(i)(1:lc)
  101 FORMAT(5X,'Origin of orbital elements: ',A)
      WRITE(unit,102)
  102 FORMAT(/)
!     CALL ephemc(unit,elem(i)%coo,elem(i)%t,elem(i)%coord,elem_unc(i)%g,defcov(i),teph1,teph2,dteph,mass(i),      &
!    &            elem(i)%h_mag,elem(i)%g_mag,idsta,ephtsc,ephfld)

      CALL ephemc(unit,elem(i),elem_unc(i),defcov(i),teph1,teph2,dteph,idsta,ephtsc,ephfld)

    1 END DO

      END

!
!  *****************************************************************
!  *                                                               *
!  *                         O F I N O D                           *
!  *                                                               *
!  *                Auxiliary routine for ORBFIT:                  *
!  *                 initial orbit determination                   *
!  *                                                               *
!  *****************************************************************
!
! INPUT:    UNIREP    -  FORTRAN unit for report
!           ELEOUT    -  Name of file for output of orbital elements
!           OPT       -  Initial orbit determination option
!           NAME      -  Object names
!           NAMOF     -  Names of observation files (w/o extension)
!           DEFORB    -  Orbit definition flag
!           DEFCN     -  Tells whether covariance/normal matrices
!                            are defined
!           DIR       -  Directory containing observation/residual file
!           NOBJ      -  Number of objects
!           OBS, OBSW -  OBSERVATIONS
!           N         -  Number of observations for each object
!           NT        -  Total number of observations
!           IP1       -  Pointer to first observation for each object
!
! IN/OUT:   ELEM      -  Orbital elements
!           COMELE    -  Comment on orbital elements
!
! Input variables DEFORB, DEFCN, SEL are modified by the routine
!
      SUBROUTINE ofinod(unirep,eleout,opt,name,namof,deforb,defcn,dir,nobj,obs,obsw,n,nt,ip1,elem,comele,error_model)
      USE astrometric_observations
      USE orbit_elements
      IMPLICIT NONE

      INTEGER unirep,opt,nobj,ip1(nobj)
! OBSERVATIONS
      INTEGER nt,n(nobj)
! new data types
      TYPE(ast_obs),DIMENSION(nt) :: obs
      TYPE(ast_wbsr),DIMENSION(nt) :: obsw
      CHARACTER(LEN=*) ::  error_model ! weighing model

      TYPE(orbit_elem),DIMENSION(nobj) :: elem

      LOGICAL deforb(nobj),defcn(nobj)
      CHARACTER*(*) eleout,name(nobj),namof(nobj),dir(nobj)
      CHARACTER*(*) comele(nobj)

      INTEGER i,ln,uniele,lc,ld
      DOUBLE PRECISION cove(6,6),nore(6,6)
      CHARACTER*150 rwofil
      LOGICAL opnd,doit,fail
      INCLUDE 'parcmc.h90'

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
      CALL io_det(unirep,rwofil,name(i),obs(ip1(i):ip1(i)+n(i)-1),       &
     &            obsw(ip1(i):ip1(i)+n(i)-1),n(i),error_model,2,         &
     &            elem(i)%coord,elem(i)%t,elem(i)%coo,comele(i),fail)
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
          cove=0.d0
          nore=0.d0
          CALL wromlr(uniele,name(i),elem(i)%coord,elem(i)%coo,elem(i)%t,cove, &
     &                .false.,nore,.false.,elem(i)%h_mag,elem(i)%g_mag,0.D0)
      END IF
  300 FORMAT(A,1X,A)

    1 END DO

      IF(opnd) CALL filclo(uniele,' ')

      END SUBROUTINE ofinod

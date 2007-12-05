! Copyright (C) 1997-1999 by Mario Carpino (carpino@brera.mi.astro.it)
! Version: February 11, 1999
! Version: November 4, 1999 (new inobs call, SRC)
! ---------------------------------------------------------------------
!
!  *****************************************************************
!  *                                                               *
!  *                         O F I O B S                           *
!  *                                                               *
!  *                Auxiliary routine for ORBFIT:                  *
!  *                    input of observations                      *
!  *                                                               *
!  *****************************************************************
!
! INPUT:    UNIREP    -  FORTRAN unit for report
!           NAME      -  Object file names
!           NAMOF     -  Names of observation files (w/o extension)
!           DIR       -  Directory containing observation/residual files
!           NOBJ      -  Number of objects
!           NOBX      -  Max number of observations
!           OPUPD     -  Option: update weight file
!
! OUTPUT:   N         -  Number of observations for each object
!           NT        -  Total number of observations
!           IP1       -  Pointer to first observation for each object
!           IP2       -  Pointer to first weight for each object
!           OBS,OBSW  -  Observations
!
      SUBROUTINE ofiobs(unirep,name,namof,dir,nobj,opupd,n,nt,     &
     &                  ip1,ip2,obs,obsw,error_model)
      USE astrometric_observations
      IMPLICIT NONE

      INCLUDE 'parobx.h90'

      INTEGER unirep,nobj,ip1(nobj),ip2(nobj)
      INTEGER obscod(nobx),sel(nobx),iobs(nobx),opupd
      CHARACTER*(*) name(nobj),namof(nobj),dir(nobj)

 ! OBSERVATIONS
      INTEGER nt,n(nobj)
! new data types
      TYPE(ast_obs),  DIMENSION(nobx) :: obs
      TYPE(ast_wbsr), DIMENSION(nobx) :: obsw
      CHARACTER(LEN=*),INTENT(IN) ::  error_model ! weighing model
      LOGICAL precob,ok

      INTEGER i,ln,ld
      LOGICAL change

      INTEGER lench
      EXTERNAL lench

      WRITE(unirep,120)
  120 FORMAT('Input of observations:')

      nt=0
      precob=.false.

      DO 1 i=1,nobj
      ip1(i)=nt+1
      ip2(i)=2*nt+1
! input data
      CALL input_obs(dir(i),namof(i),precob,error_model,ok,obs(ip1(i):),obsw(ip1(i):),n(i),unirep,change)
      IF(.NOT.ok) STOP '**** ofiobs: abnormal end ****'
      nt=nt+n(i)
      ld=lench(dir(i))
      ln=lench(name(i))
      WRITE(unirep,200) n(i),name(i)(1:ln)
  200 FORMAT(I9,' observations read for object ',A)
    1 END DO

      END SUBROUTINE ofiobs

! Copyright (C) 1999 by Mario Carpino (carpino@brera.mi.astro.it)
! Version: February 8, 1999
! ---------------------------------------------------------------------
!
!  *****************************************************************
!  *                                                               *
!  *                         O F I N I P                           *
!  *                                                               *
!  *       Initializations required for orbit propagation          *
!  *                                                               *
!  *****************************************************************
!
! INPUT:    RUN       -  Run name
!
      SUBROUTINE ofinip(run)
      USE output_control
      IMPLICIT NONE

      CHARACTER*(*) run

      INTEGER lf
      CHARACTER*80 file

      INTEGER lench
      EXTERNAL lench

      lf=lench(run)

! Output file for errors
      file=run(1:lf)//'.err'
      CALL filopn(ierrou,file,'UNKNOWN')
      numerr=0

! Output file for close approaches
      file=run(1:lf)//'.clo'
      CALL filopn(iuncla,file,'UNKNOWN')
      numcla=0

! Output file for propagator parameters
      file=run(1:lf)//'.pro'
      CALL filopn(ipirip,file,'UNKNOWN')

      END

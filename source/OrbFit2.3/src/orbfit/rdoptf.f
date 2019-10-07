* Copyright (C) 1997-1998 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: December 14, 1998
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         R D O P T F                           *
*  *                                                               *
*  *      Read options for least squares orbital fit(ORBFIT)       *
*  *                                                               *
*  *****************************************************************
*
      SUBROUTINE rdoptf
      IMPLICIT NONE

* Common blocks to be initialized:
      INCLUDE 'comlsf.h'

      LOGICAL found,fail1,fail

      fail=.false.

      delcr=1.d-5
      CALL rdnrea('lsfit.','conv_cntr',delcr,.false.,found,
     +            fail1,fail)

      IF(fail) STOP '**** rdoptf: abnormal end ****'

      iiclsf=36

      END

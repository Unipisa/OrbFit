* Copyright (C) 1997-1998 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: December 7, 1998
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         R D O P T I                           *
*  *                                                               *
*  *        Read options for orbit identification (ORBFIT)         *
*  *                                                               *
*  *****************************************************************
*
      SUBROUTINE rdopti
      IMPLICIT NONE

* Common blocks to be initialized:
      INCLUDE 'comidn.h'

      LOGICAL found,fail1,fail

      fail=.false.

      amfit=.true.
      CALL rdnlog('ident.','aM_fit',amfit,.false.,found,
     +            fail1,fail)
      delcr=1.d-5
      CALL rdnrea('ident.','conv_cntr',delcr,.false.,found,
     +            fail1,fail)

      IF(fail) STOP '**** rdopti: abnormal end ****'

      iicidn=36

      END

* Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: December 15, 1997
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         R D O P T E                           *
*  *                                                               *
*  *         Read options for ephemeris generation (ORBFIT)        *
*  *                                                               *
*  *****************************************************************
*
      SUBROUTINE rdopte
      IMPLICIT NONE

* Common blocks to be initialized:
      INCLUDE 'comeph.h'

      INTEGER mjd,mjde
      DOUBLE PRECISION sec,sece
      LOGICAL found,fail1,fail
      CHARACTER tmp*100,tsc*3

      fail=.false.

* List of objects
      CALL rdmint('ephem.','objects',kepobj,nepobj,3,'3',.false.,
     +            found,fail1,fail)
      IF(.NOT.found) THEN
          nepobj=3
          kepobj(1)=1
          kepobj(2)=2
          kepobj(3)=3
      END IF

* Limits and stepsize of ephemeris
      CALL rdntim('ephem.epoch.','start',tmp,mjd,sec,tsc,.true.,
     +            found,fail1,fail)
      IF(.NOT.fail1) THEN
          CALL cnvtim(mjd,sec,tsc,mjde,sece,'TDT')
          teph1=mjde+sece/86400.d0
      END IF
      CALL rdntim('ephem.epoch.','end',tmp,mjd,sec,tsc,.true.,
     +            found,fail1,fail)
      IF(.NOT.fail1) THEN
          CALL cnvtim(mjd,sec,tsc,mjde,sece,'TDT')
          teph2=mjde+sece/86400.d0
      END IF
      dteph=1
      CALL rdnrea('ephem.','step',dteph,.false.,found,fail1,fail)

* Observatory code
      idsta=500
      CALL rdnint('ephem.','obscode',idsta,.false.,found,fail1,fail)

* Timescale
      ephtsc='TDT'
      CALL rdncha('ephem.','timescale',ephtsc,.false.,found,fail1,fail)

* Output fields
      ephfld='cal,coord,delta,r,elong,phase,mag'
      CALL rdncha('ephem.','fields',ephfld,.false.,found,fail1,fail)

      IF(fail) STOP '**** rdopte: abnormal end ****'

      iiceph=36

      END

* Copyright (C) 1999 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: November 26, 1999
* ---------------------------------------------------------------------
* RMS classes of astrometric observations: additional information
*
* idcl        -  Class progressive number in index
* orstep      -  Decision step (RA/DEC) :
*                    0 = no information
*                    1 = single-station class
*                    2 = multi-station class
*                    3 = default
* tdtlim      -  Class limits (MJD, TDT)
*
      INTEGER idcl(2),orstep(2)
      COMMON/cmaow1/idcl,orstep
      DOUBLE PRECISION tdtlim(2,2)
      COMMON/cmaow2/tdtlim

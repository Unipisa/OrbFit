* Copyright (C) 1998 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: December 3, 1998
* ---------------------------------------------------------------------
* OPTIONS FOR ITERATION CONTROL IN DIFCOR
*
* batch       -  Logical to proceed noninteractively
* itmax       -  Max number of iterations
* itgmax      -  Max no. of iter. with non decreasing sigma0
* divrat      -  Divergency/non-convergency control
* iicdif      -  Initialization check
*
      INTEGER itmax,itgmax,iicdif
      DOUBLE PRECISION divrat
      LOGICAL batch
      COMMON/cmdif1/itmax,itgmax,iicdif
      COMMON/cmdif2/divrat
      COMMON/cmdif3/batch

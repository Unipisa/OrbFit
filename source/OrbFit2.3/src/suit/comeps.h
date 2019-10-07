* Copyright (C) 1999 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: November 11, 1999
* ---------------------------------------------------------------------
* ROUNDING-OFF PRECISION
*
* eps      -  Max value of delta for which 1+delta = 1
* l2eps    -  log2(eps), namely eps = 2**l2eps
* iiceps   -  Initialization check
*
      DOUBLE PRECISION eps
      COMMON/comep1/eps
      INTEGER l2eps,iiceps
      COMMON/comep2/l2eps,iiceps

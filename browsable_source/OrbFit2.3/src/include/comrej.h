* Copyright (C) 1998 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: December 4, 1998
* ---------------------------------------------------------------------
* OPTIONS FOR OUTLIER REJECTION/RECOVERY
*
* autrej      -  Automatic rejection flag
* rejopp      -  Permit the rejection of an entire opp.
* x2rej       -  Chi-square threshold for rejection
* x2rec       -  Chi-square threshold for recovery
* x2frac      -  Fraction of max chi-square at which start rejection
*                in the "large residual" regime
* delrej      -  Convergency control (correction norm)
* fomax       -  Max fraction of outliers (total)
* itmaxr      -  Max number of iterations
* iicrej      -  Initialization check
*
      DOUBLE PRECISION x2rej,x2rec,x2frac,delrej,fomax,x2mrej,x2mrec
      INTEGER itmaxr,iicrej
      LOGICAL autrej,rejopp
      COMMON/cmrej1/x2rej,x2rec,x2frac,delrej,fomax
      COMMON/cmrej2/itmaxr,iicrej
      COMMON/cmrej3/autrej,rejopp
      COMMON/magrej/x2mrej,x2mrec

* Copyright (C) 1999 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: November 12, 1999
* ---------------------------------------------------------------------
* IERS time series
*
* ieidir      -  IERS input directory
* ieidil      -  Length of IERS input directory name
* ieilis      -  IERS input file list
* ieilil      -  Length of IERS input file list name
* isamp       -  Undersampling factor
* tiers       -  Time of data points (MJD, TDT)
* xiers       -  Value of data points:
*                    xiers(1) = X pole (arcsec)
*                    xiers(2) = Y pole (arcsec)
*                    xiers(3) = TDT-UT1 (s)
*                    xiers(4) = Dpsi (arcsec)
*                    xiers(5) = Depsilon (arcsec)
* niers       -  Number of data points
* iutsmo      -  Smoothing used for UT1 (xiers(i,3)) :
*                   0 = none
*                   1 = UT1R
*                   2 = UT1S
* cciera      -  Apply consistency correction (flag)
* flcier      -  IERS consistency correction file
* cncor0(5)   -  Consistency correction (constant term)
* cncor1(5)   -  Consistency correction (linear term)
* cnep0       -  "Zero" epoch for consistency correction (Bess. year)
* nlpler      -  NL (interpolation length) for pcwlgi
* nvpler      -  NV (length of empty zone) for pcwlgi
* nspler      -  NS (length of superposition zone) for pcwlgi
* nsmopl      -  Order of smoothing for pcwlgi
* rmserx(5)   -  Max interpolation RMS error allowed for pcwlgi
* extra       -  Allow extrapolation
* blause      -  Use Bulletin A for extension in the future of real data
* blafil      -  Name of input file for Bulletin A
* tint1       -  Starting time for interpolation (end of extrapolation before data start)
* tint2       -  Ending time for interpolation (start of extrapolation after data end)
* dutd        -  Average time derivative of TDT-UT1
* utd1        -  Starting value for extrapolation of TDT-UT1 (tjme<tint1)
* utd2        -  Starting value for extrapolation of TDT-UT1 (tjme>tint2)
* iicier      -  Initialization check
*
      CHARACTER*100 ieidir,ieilis,flcier,blafil
      COMMON/cmier1/ieidir,ieilis,flcier,blafil
      INTEGER ieidil,ieilil,isamp,niers,iutsmo,iicier
      INTEGER nlpler,nvpler,nspler,nsmopl
      COMMON/cmier2/ieidil,ieilil,isamp,niers,iutsmo,nlpler,nvpler,
     +              nspler,nsmopl,iicier
      DOUBLE PRECISION tiers(niersx),xiers(niersx,5)
      DOUBLE PRECISION cncor0(5),cncor1(5),cnep0,rmserx(5),tint1,tint2
      DOUBLE PRECISION dutd,utd1,utd2
      COMMON/cmier3/tiers,xiers,cncor0,cncor1,cnep0,rmserx,tint1,tint2,
     +              dutd,utd1,utd2
      LOGICAL cciera,extra,blause
      COMMON/cmier4/cciera,extra,blause

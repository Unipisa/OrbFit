* Copyright (C) 1997 by A.Milani and M.Carpino
* Version: February 22, 1997
* Trigonometric constants
      double precision dpig,pig,radeg,radsec,degrad,secrad,radh,hrad
      parameter (dpig=6.28318530717958648d0)
      parameter (pig=dpig/2.0d0)
* Radians from degrees
      parameter (radeg=dpig/360.0d0)
* Radians from arcseconds
      parameter (radsec=radeg/3600.0d0)
* Degrees from radians
      parameter (degrad=360.0d0/dpig)
* Arcseconds from radians
      parameter (secrad=3600.0d0/radeg)
* Radians from hours
      parameter (radh=dpig/24.0d0)
* Hours from radians
      parameter (hrad=24.0d0/dpig)

* Copyright (C) 2000 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: November 27, 2000
* ---------------------------------------------------------------------
* Models of autocorrelation functions for all the observatories
* pto2f       -  pointer from observatory code to function
* nfo         -  number of functions for the model of each station
* nfunt       -  total number of functions
* npart       -  total number of parameters
* nparf       -  number of parameters for each function
* nparo       -  number of parameters for each observatory
* kfun        -  function integer code
* ptf2p       -  pointer from function to parameters
* par         -  parameter values
* aprmx       -  max a-priori RMS (arcsec)
* minapw      -  min a-priori weight (rad**(-2))
* maxdst      -  max "decision step" (1=specific obs; 2=mixed class)
* iiccor      -  initialization check
*
      INTEGER nfunt,npart,pto2f(obsc1:obsc2,2,2),nfo(obsc1:obsc2,2)
      INTEGER nparf(nparx),kfun(nparx),ptf2p(2,nparx)
      INTEGER nparo(obsc1:obsc2,2),maxdst,iiccor
      COMMON/cmcor1/nfunt,npart,nparf,kfun,pto2f,nfo,ptf2p,nparo,
     +              maxdst,iiccor
      DOUBLE PRECISION par(nparx),aprmx,minapw
      COMMON/cmcor2/par,aprmx,minapw

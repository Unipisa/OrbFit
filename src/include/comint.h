* controls for numerical integrator: Everhart
      integer llev
      double precision hev
      common/integr/hev,llev
*  controls for numerical integration: rkimp and multistep
      double precision hms,epms,eprk,deltos,error,eprk_c
      integer mms,isrk,lit1,lit2,imet,iusci,icha,icmet,lit1_r,lit2_r,
     +    lit1_c,lit2_c,isrk_c
      common/numer/hms,epms,eprk,deltos,error, eprk_c,mms,isrk,
     +         lit1,lit2,lit1_r,lit2_r,lit1_c,lit2_c,isrk_c,imet,
     +         iusci,icha,icmet

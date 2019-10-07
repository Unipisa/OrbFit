* Copyright (C) 1996 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: August 10, 1996
* ---------------------------------------------------------------------
* INPUT NAMELIST
*
* nne         -  Number of namelist items
* krecnm(i)   -  Record within input file (i=1,nne)
* kuord(i)    -  Order of use (i=1,nne)
* krtyp(i)    -  Keyword type:
*                      1 - integer
*                      2 - real
*                      3 - character string
*                      4 - logical
*                      5 - MJD
* kuorl       -  Last order of use
* iicnam      -  Initialization check
* keys(i)     -  Namelist keyword fields (i=1,nne)
* vals(i)     -  Namelist value fields (i=1,nne)
* namif(i)    -  Namelist input file (i=1,nne)
*
      INTEGER nne,krecnm,kuord,krtyp,kuorl,iicnam
      CHARACTER*(lckx) keys
      CHARACTER*(lcvx) vals
      CHARACTER*(lcfx) namif
      COMMON/cnam1/nne,krecnm(nnex),kuord(nnex),krtyp(nnex),
     +             kuorl,iicnam
      COMMON/cnam2/keys(nnex)
      COMMON/cnam3/vals(nnex)
      COMMON/cnam4/namif(nnex)

* Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: March 8, 1997
* ---------------------------------------------------------------------
* SIMPLIFIED FILE-HEADER NAMELIST
*
* nfne        -  Number of namelist items
* krcfnm(i)   -  Record within input file (i=1,nfne)
* kuorf(i)    -  Order of use (i=1,nfne)
* kuorlf      -  Last order of use
* iicfnm      -  Initialization check
* keysf(i)    -  Namelist keyword fields (i=1,nfne)
* valsf(i)    -  Namelist value fields (i=1,nfne)
* nmif        -  Namelist input file
* hnfuni      -  Namelist input unit
*
      INTEGER nfne,krcfnm,kuorf,kuorlf,hnfuni,iicfnm
      CHARACTER keysf*(lckx),valsf*(lcvx),nmif*(lcfx)
      COMMON/cmfnm1/nfne,krcfnm(nfnex),kuorf(nfnex),kuorlf,
     +              hnfuni,iicfnm
      COMMON/cmfnm2/keysf(nfnex)
      COMMON/cmfnm3/valsf(nfnex)
      COMMON/cmfnm4/nmif

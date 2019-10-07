* Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: March 8, 1997
* ---------------------------------------------------------------------
* Information on opened orbital element file
*
* orbunt      -  FORTRAN unit
* orbnr       -  Number of records read so far
* orbfn       -  File name
* dfrsty      -  default reference system type
* dfrsep      -  default reference system epoch
* rectyp      -  record type
* deltyp      -  default orbital element type
* depstr      -  default epoch (character string)
* dept0       -  default epoch (MJD, TDT)
* deft0       -  is dept0 defined?
* nxtend      -  end encountered at previous call
* iicorb      -  initialization check
*
      INTEGER orbunt,orbnr,iicorb
      CHARACTER*100 orbfn,depstr
      DOUBLE PRECISION dept0
      CHARACTER*4 dfrsty,rectyp,deltyp
      CHARACTER*10 dfrsep
      LOGICAL deft0,nxtend
      COMMON/cmorb1/orbunt,orbnr,iicorb
      COMMON/cmorb2/orbfn,depstr
      COMMON/cmorb3/dept0
      COMMON/cmorb4/dfrsty,rectyp,deltyp
      COMMON/cmorb5/dfrsep
      COMMON/cmorb6/deft0,nxtend

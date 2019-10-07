* Copyright (C) 1996 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: August 10, 1996
* ---------------------------------------------------------------------
* File names and assigned units
*
* filnam      -  File names
* allunt      -  Allocation indicator
* iicfil      -  Initialization check
*
      INTEGER iunf1,iunf2,iicfil
      PARAMETER (iunf1=50)
      PARAMETER (iunf2=99)
      CHARACTER*80 filnam(iunf1:iunf2)
      LOGICAL allunt(iunf1:iunf2)

      COMMON/cmfil1/filnam
      COMMON/cmfil2/allunt
      COMMON/cmfil3/iicfil

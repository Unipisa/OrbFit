* Copyright (C) 1997 by M.Carpino
* Version: February 22, 1997
*
* Binary (random access) ephemeris
*
* filbep   -  Name of the binary ephemeris file
* nbep     -  Number of bodies
* masbep   -  Masses
*
      CHARACTER*80 filbep
      INTEGER nbep
      DOUBLE PRECISION masbep
      COMMON/cmbep1/filbep
      COMMON/cmbep2/nbep
      COMMON/cmbep3/masbep(nbepx)

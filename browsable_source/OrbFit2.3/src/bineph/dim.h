* Max number of massive bodies (asteroids)
      INTEGER nbmx
      PARAMETER (nbmx=5)
* Max number of massless asteroids
      INTEGER nbax
      PARAMETER (nbax=1)
* Max number of planets from JPL ephemerides
      INTEGER nbjx
      PARAMETER (nbjx=10)
* Max number of asteroids in the state vector
      INTEGER nbtx
      PARAMETER (nbtx=nbax+nbmx)

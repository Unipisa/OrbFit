c selection of asteroids; mostly to avoid self-perturbations
c requires parbep.h to have dimensions
c note that asteroid are identified by number, since only
c numbered asteroids can have known mass
      INTEGER astid,iatrue
      COMMON/selast/astid(nbepx),iatrue

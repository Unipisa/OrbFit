c maximum number of planets (from JPL epehemrides)
      integer nplax,nmassx
      parameter (nplax=10)
c warning: masses.h depends upon parbep.h
c maximum number of perturbing bodies (including massive asteroids)
      parameter (nmassx=nplax+nbepx)
      double precision gm0,gmu,gm,gmse
c  masses of the planets, flags to find them in order
      integer npla,nmass,itarg,listpl(12)
      character*30 ordnam
      common/masse/gm0,gmu,gm(nmassx),gmse,npla,nmass,listpl,
     +itarg(nplax)
      common/planam/ordnam(nmassx)
c  distance below which close approach is stored
      double precision dmin
      common/ddis/dmin(nmassx)
c  asteroid masses: see combep.h

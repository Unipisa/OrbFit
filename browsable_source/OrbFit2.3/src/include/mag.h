c ===============================================================
c mag.h header file for MAGEST routine
c ===============================================================
c difference apparent minus absolute magnitude (to compute magnitude)
c dmagn is in original observations order, 
c dmagns is sorted as required by propagator 
c warning this header requires parox.h
c also slope parameter known for current object
c phav,phavs is phase
c dsunv,dsunvs distance from the sun
c disv,disvs distance from Earth
c adotmv,adots ddotmv,ddots proper motion
      double precision dmagn(nobx),dmagns(nobx),gmagc
      DOUBLE PRECISION dsuna(nobx),dsunas(nobx),disa(nobx),disas(nobx)
     +  ,phaa(nobx),phaas(nobx),adots(nobx),ddots(nobx),
     +   adotmv(nobx),ddotmv(nobx)
      common/mag/dmagn,dmagns,gmagc,phaa,phaas,dsuna,dsunas,disa,disas,
     +    adotmv,adots,ddotmv,ddots

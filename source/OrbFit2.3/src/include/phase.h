c ===============================================================
c phase.h header file for semi-public interface of PREOBX
c ===============================================================
c phase, distance to Earth, distance to Sun (to compute magnitude)
c elongation, galactic latitude, apparent motion
      double precision pha,dis,dsun,elo,gallat,adot,ddot
      common/phase/pha,dis,dsun,elo,gallat,adot,ddot
c galactic pole alpha=12h 51.3m, delta=27d 7.0min
      double precision algal,degal,gax,gay,gaz
      parameter (algal=3.36543113015807d0)
      parameter (degal=0.47327511549913d0)
      parameter (gax=-0.86787505968543d0)
      parameter (gay=-0.19757464383054d0)
      parameter (gaz=0.45580384036475d0)
c from http://www.seds.org/~spider/spider/ScholarX/coords.html#galactic
c The galactic north pole is at RA = 12:51.4, Dec = +27:07 (2000.0)
c the galactic center at RA = 17:45.6, Dec = -28:56 (2000.0). 
c The inclination of the galactic equator to Earth''s equator is thus 62.9 deg. 
c The intersection, or node line of the two equators is at 
c RA = 18:51.4, Dec = 0:00 (2000.0), and at l = 33 deg, b=0.

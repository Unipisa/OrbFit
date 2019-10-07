c =========MODULE moid_compute==============
c CONTAINS
c ROUTINES
c             nomoid
c             nodedi
c 
c  HEADERS
c moid_compute.o: 
c            parbep.h masses.h masses Earth and Moon
c            sunmass.h public
c
c =========================================
c NOMOID - compute moid using Gronchi's routine
c nodal distances and Minimum Orbital Intersection Distance
c with respect to the Earth
c for an orbit, given in equinoctal elements, at a given time
c =========================================
      SUBROUTINE nomoid(t0,eq0,moid,icon,dnp,dnm)
      IMPLICIT NONE
c ===========INPUT=====================
c elements of asteroid (equinoctal), epoch (MJD)
      DOUBLE PRECISION eq0(6),t0
c ===========OUTPUT=================================
c MOID, iteration count and succ. flag
c ascending node distance,  descending node distance
      DOUBLE PRECISION moid,dnp,dnm
c archaic argument
      INTEGER icon
c ==========END INTERFACE================================
c get gms (GM_sun) ans gmse (G(M_sun+M_ear))
      INCLUDE 'parbep.h'
      INCLUDE 'masses.h'
      INCLUDE 'sunmass.h'
c elements of Earth (equinoctal*, cartesian coordinates
      DOUBLE PRECISION eqp(6),xast(6),xea(6),enne
c cartesian coordinates asteroid and planet at minimum
      DOUBLE PRECISION cmin(6,20),cplmin(6,20)
c     SQUARED DISTANCE function 
      DOUBLE PRECISION d2(20)
c     number of relative minima found
      INTEGER nummin
c archaic parameter
      icon=0
c ======================================================
c get Earth elements
      CALL earth(t0,eqp)
c compute nodal distances
      CALL nodedi(eq0,eqp,dnp,dnm)
c transform to cartesian coordinates
      CALL coocha(eq0,'EQU',gms,xast,'CAR',enne)
      CALL coocha(eqp,'EQU',gmse,xea,'CAR',enne)
c compute moid      
      CALL compute_minima(xast,xea,3,cmin,cplmin,d2,nummin)
      moid=sqrt(d2(1))
      RETURN
      END
c ==========================================
c NODEDI
c nodal distances of two elliptic orbits
c ==========================================
      SUBROUTINE nodedi(eq,eqp,dnp,dnm)
      IMPLICIT NONE
c ===========INPUT=====================
c elements of asteroid, of Earth (equinoctal)
      DOUBLE PRECISION eq(6),eqp(6)
c ===========OUTPUT=================================
c output ascending node distance,  descending node distance
      DOUBLE PRECISION dnp,dnm
c ==========END INTERFACE================================
      DOUBLE PRECISION c(3),cp(3),x(6),xp(6),vlenz(3),vlenzp(3)
      INCLUDE 'parbep.h'
      INCLUDE 'masses.h'
      INCLUDE 'sunmass.h'
      DOUBLE PRECISION enne,ennep,vnod(3),vnl,ome,omep
      DOUBLE PRECISION ecc,eccp,cosf,cosfp,rp,rm,rpe,rme
      INTEGER i
      DOUBLE PRECISION prscal,vsize
c cartesian coordinates      
      CALL coocha(eq,'EQU',gms,x,'CAR',enne)
      CALL coocha(eqp,'EQU',gmse,xp,'CAR',ennep)
c  angular momentum      
      CALL prvec(x,x(4),c)
      CALL prvec(xp,xp(4),cp)
c ascending node
      CALL prvec(cp,c,vnod)
      vnl=vsize(vnod)
c  angular momentum unit vector, Lenz vector
      CALL  prvec(x(4),c,vlenz)
      DO i=1,3
        vlenz(i)=vlenz(i)/gms-x(i)/vsize(x)
      ENDDO
      ecc=vsize(vlenz)
      CALL  prvec(xp(4),cp,vlenzp)
      DO i=1,3
        vlenzp(i)=vlenzp(i)/gmse-xp(i)/vsize(xp)
      ENDDO
      eccp=vsize(vlenzp)
c true anomaly at mutual node= - arg. of perihelion
      cosf=prscal(vnod,vlenz)/(vnl*ecc)
      ome=acos(cosf)
      cosfp=prscal(vnod,vlenzp)/(vnl*eccp)
      omep=acos(cosfp)
c nodal points and distances
      rp=eq(1)*(1.d0-eq(2)**2-eq(3)**2)/(1.d0+ecc*cosf)
      rm=eq(1)*(1.d0-eq(2)**2-eq(3)**2)/(1.d0-ecc*cosf)  
      rpe=eqp(1)*(1.d0-eqp(2)**2-eqp(3)**2)/(1.d0+eccp*cosfp)
      rme=eqp(1)*(1.d0-eqp(2)**2-eqp(3)**2)/(1.d0-eccp*cosfp)            
      dnp=rp-rpe
      dnm=rm-rme
c      WRITE(*,*)rp,rpe,dnp,rm,rme,dnm,cosf,cosfp,ome,omep
      RETURN
      END


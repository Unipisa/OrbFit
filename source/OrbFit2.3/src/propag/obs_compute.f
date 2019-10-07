c=================MODULE obs_compute=============
c CONTAINS
c ROUTINES
c                  alfdel     optical observations
c                    ossdif       "
c                  appmag     apparent magnitude
c                  alfdel2    optical obs. only for attribute
c                    ossdif2       "
c                  aber1      aberration
c                  rrdot      radar observations
c                    deltau        "
c                    deldop1       "
c                    deldop2       "
c                  
c  HEADERS
c obs_compute.o:   trig.h : public
c                  model.h : ossdif, ossdif2 (flags for  vlight.h 
c                  jplhdr.h : rrdot receives data from jpl_ephem.f
c                  radius.h : rrdot, receives data from obs_wei.f
c                  sunmass.h : rrdot deltau deldop1 
c                  phase.h : hidden output (phase, proper motion,...)
c                   
c =====================================================================
c ALFDEL
c =====================================================================
c Computation of alpha and delta and their derivatives
c =====================================================================
c
c Input
c    tk: epoch time
c    tauj: observation time
c    ioj: station code
c    east: equinoctal orbital elements vector (a,h,k,p,q,lambda) 
c                    at time tk
c    ider: flag for derivatives options:
c        <1 no derivatives
c        <2 only first deriv. of alpha,delta w.r. to (east,tk,tauj)
c        >= 2 also second partial derivatives (aproximation with 2-body case)
c    twobo: logical flag for 2-body approximation; if .true., 2-body
c        approximation (but Earth position comes from JPL ephem); if .false.,
c        all orbit propagations are full n-body
c Output
c   alj,dej: alpha,delta  computed at time tauj
c   dade,ddde:  matrices of first derivatives (if required)
c   ddade,dddde: matrices of second deriv. (if required, only in 2-body approx)
c
c ==============INTERFACE============================================
      SUBROUTINE alfdel (east,tk,tauj,iocj,alj,dej,dade,ddde,ider,twobo,
     +          ddade,dddde)
      implicit none
c flag to control two-body approximation: if .false., full n-body
      logical twobo
c times: epoch time for asteroid elements, observation time (MJD)
      double precision tk,tauj
c asteroid equinoctal elements 
      double precision east(6)
c cartesian coordinates of the Earth
      double precision xea(6)
c observations: alpha (right ascension) delta (eclination), in RAD
      double precision alj,dej
c observatory code
      integer iocj
c flag to control computation of derivatives
      integer ider
c partial derivatives of alpha, delta, w.r. to asteroid coordinates
      double precision dade(6),ddde(6)
c second derivatives of alpha, delta, w.r. to asteroid coordinates, elements
      double precision  ddade(6,6),dddde(6,6)
c =============END INTERFACE=========================================
c asteroid cartesian coordinates
      double precision xast(6)
c first, second partial derivatives of alpha, delta, w.r. to ast. coordinates
      double precision dadx(3),dddx(3),ddadx(3,3),ddddx(3,3)
c first, second derivatives of cartesian coordinates with respect to elements
      double precision dxde(6,6),ddxde(3,6,6)
c temporary variables for matrix multiplication
      double precision tmp,tmpd
c loop variables i,j=1,6; k,m,3
      integer i,j,k,m
c double precision functions
      double precision prscal
c elongation,distance to Earth, distance to Sun (to compute magnitude)
      include 'phase.h'
c *************************************************************************
c****************
c   static memory not required
c****************
c Orbit propagation:
      if(twobo)then
c 2 body version
         call propa2 (tk,east,tauj,xast,xea,ider,dxde,ddxde)
      else
c full n-body numerical integration
         call propag (tk,east,tauj,xast,xea,ider,dxde,ddxde)
      endif
c Computation of observations
      call ossdif(xast,xea,tauj,iocj,alj,dej,ider,dadx,
     +              dddx,ddadx,ddddx)
      if(ider.lt.1)return
c Derivatives of $\alpha$ and $\delta$ w.r. to time
c currently not in use
c     dade(8)=prscal(dadx,d(4))
c     ddde(8)=prscal(dddx,d(4))
c     dade(7)=-dade(8)
c     ddde(7)=-ddde(8)
c derivatives with respect to equinoctal elements
      do 11 j=1,6
        dade(j)=prscal(dadx,dxde(1,j))
        ddde(j)=prscal(dddx,dxde(1,j))
 11   continue
c Chain rule for second derivatives (2-body case): we use eq. (6.2)
      if(ider.lt.2)return
      do 20 i=1,6
        do 23 j=1,6
          tmp=0.d0
          tmpd=0.d0
          do 21 k=1,3
            tmp=tmp+dadx(k)*ddxde(k,i,j)
            tmpd=tmpd+dddx(k)*ddxde(k,i,j)
            do 22 m=1,3
              tmp=tmp+dxde(m,i)*ddadx(m,k)*dxde(k,j)
              tmpd=tmpd+dxde(m,i)*ddddx(m,k)*dxde(k,j)
 22         continue
 21       continue
          dddde(i,j)=tmpd
          ddade(i,j)=tmp
 23     continue
 20   continue
      return
      end
c =====================================================================
c OSSDIF vers. 1.3 (equatorial coordinates for observations)
c A. Milani, June 14, 1997
c =====================================================================
c Corrections to observations
c =====================================================================
c Input
c xast asteroid cartesian coordinates at time tauj
c xea Earth cartesian coordinates at time tauj
c     both in the ecliptic system
c tauj observation time
c ioc station code
c ider flag for derivatives options:
c       =0 no derivatives
c       <2 only first deriv. of $\alpha,\delta$ w.r. 
c to positions
c       great equal 2 also second partial derivatives 
c (approximation with 2-body case)
c
c Output
c alj,dej,alpha,delta computed at time tauj (in the equatorial system)
c (if required)
c dadx,dddx matrices of first derivatives w. r. to ecliptic positions 
c (if required)
c ddadx,ddddx matrices of second deriv., 2-b approximation
c =====================================================================
      SUBROUTINE ossdif(xast,xea,tauj,ioc,alj,dej,ider,
     +     dadx,dddx,ddadx,ddddx)
      implicit none
c =====================================================================
      include 'trig.h'
      include 'model.h'
c =====================================================================
c observations alpha, delta (equatorial, radians) and time
      double precision alj,dej,tauj
c cartesian coordinates, ecliptic, of asteroid, of earth, vector difference
c WARNING: even if the velocities are not always propagated, they are available
c at the time of the observations
      double precision xast(6),xea(6),d(6)
c partials of equatorial alpha, delta w.r. to cartesian ecliptic coord. 
      double precision dadx(3),dddx(3),ddadx(3,3),ddddx(3,3)
c space for cal.
      double precision atmp(3,3)
c topocentric position of the observatory, obs. code
      double precision xo(3),vo(3)
      integer idst,ioc
c elongation,distance to Earth, distance to Sun (to compute magnitude)
      include 'phase.h'
c phase and elongation cosine
      double precision cospha,coselo
c control on no. derivatives, bits in mantissa
      integer ider,n
c rounding off, auxiliary var.
      double precision errm,dz
c real functions
      double precision roff,vsize,prscal
c integer loop index
      integer i,ii,ij,j
c rotation matrices, rotated vectors
      double precision rot(3,3),rotinv(3,3),deq(6),tmp(3)
c aux. var for computation of second derivatives
      double precision x,y,z
      double precision den,x2,y2,z2,x4,y4
      double precision x2y2,x2z2,y2z2
c =====================================================================
      integer lflag
*********************************
* static memory allocation only for:
      save lflag,errm, rot,rotinv
********************************* 
      data lflag/0/
      if(lflag.eq.0)then
         errm=roff(n)
* Change of reference system EQUM00 ---> ECLM00
c Rotation matrix from the reference system in which the
c observations are given (mean equatorial J2000) to the reference 
c system in which the orbit is computed
         call rotpn(rot,'EQUM','J2000',0.d0,'ECLM','J2000',0.d0)
         do 1 ii=1,3
           do 2 ij=1,3
             rotinv(ii,ij)=rot(ij,ii)
 2         continue
 1       continue
         lflag=1
      endif
c =====================================================================
c Difference vector
      call vdiff (xast,xea,d)
      call vdiff(xast(4),xea(4),d(4))
c =====================================================================
c Displacement of the station with respect to the center of the Earth
      if(istat.gt.0)then
         idst=ioc
         IF(idst.ne.500)THEN
            call pvobs(tauj,idst,xo,vo)
            call vdiff(d,xo,d)
            call vdiff(d(4),vo,d(4))
         ENDIF
      endif
c =====================================================================
c Aberration (only time delay)
      if(iaber.gt.0)then
         call aber1(d,xast(4),d)
      endif
c =====================================================================
c Computation of solar distance, earth distance, phase, elongation
      dsun=vsize(xast)
      dis=vsize(d)
      cospha=prscal(d,xast)/(dis*dsun)
      pha=acos(cospha)
      coselo=-prscal(d,xea)/(dis*vsize(xea))
      elo=acos(coselo)
c =====================================================================
c rotation to the equatorial reference system
      call prodmv(deq,rotinv,d)
      call prodmv(deq(4),rotinv,d(4))
c trick to change as little as possible from vers. 1.2 to 1.3
      do 3 j=1,6
        d(j)=deq(j)
 3    continue
c =====================================================================
c galactic latitude
      gallat=pig/2d0-acos((d(1)*gax+d(2)*gay+d(3)*gaz)/dis)
c =====================================================================
c Computation of observation: right ascension (radians)
      dz=d(1)**2+d(2)**2
      if (dz.le.errm) then
          alj=0.d0
      else
          alj=atan2(d(2),d(1))
          if (alj.lt.0.d0) then
              alj=alj+dpig
          endif
      endif
c Computation of observation: declination (radians)
      dej=asin(d(3)/dis)
c =====================================================================
c Computation of first derivatives of $\alpha$ and $\delta$ w.r. to positions 
c (if required): we derived eq. (2.20)
      dadx(1)=-d(2)/dz
      dadx(2)=d(1)/dz
      dadx(3)=0.d0
      dddx(1)=-d(3)*(d(1)/(sqrt(dz)*dis**2))
      dddx(2)=-d(3)*(d(2)/(sqrt(dz)*dis**2))
      dddx(3)=sqrt(dz)/dis**2
c =====================================================================
c Apparent motion: 
      adot=prscal(dadx,d(4))
      ddot=prscal(dddx,d(4))
      if(ider.eq.0)return
c =====================================================================
c rotation to the equatorial reference system
      call prodmv(tmp,rot,dadx)
      do 4 i=1,3
        dadx(i)=tmp(i)
 4    continue
      call prodmv(tmp,rot,dddx)
      do 5 i=1,3
        dddx(i)=tmp(i)
 5    continue
      if(ider.lt.2)return
c =====================================================================
c Computation of second derivatives of $\alpha$ w.r. to positions
c (if required)
      ddadx(1,1)=2.d0*d(1)*d(2)/dz**2
      ddadx(1,2)=(d(2)**2-d(1)**2)/dz**2
      ddadx(2,1)=ddadx(1,2)
      ddadx(2,2)=-ddadx(1,1)
      ddadx(3,1)=0.d0
      ddadx(3,2)=0.d0
      ddadx(3,3)=0.d0
      ddadx(2,3)=0.d0
      ddadx(1,3)=0.d0
c Computation of second derivatives of $\delta$ w.r. to positions
c (if required)
      den=1.d0/(dis**2*dz*sqrt(dz))
      x=d(1)
      y=d(2)
      z=d(3)
*
      x2=x*x
      y2=y*y
      z2=z*z
      x4=x2*x2
      y4=y2*y2
*
      x2y2=x2*y2
      x2z2=x2*z2
      y2z2=y2*z2
*
      ddddx(1,1)=z*(2.d0*x4+x2y2-y2z2-y4)*den
      ddddx(2,2)=z*(2.d0*y4+x2y2-x2z2-x4)*den
      ddddx(1,2)=x*y*z*(z2+3.d0*x2+3.d0*y2)*den
      ddddx(2,1)=ddddx(1,2)
      ddddx(3,3)=-2.d0*z*dz**2*den
      ddddx(1,3)=x*dz*(z2-x2-y2)*den
      ddddx(3,1)=ddddx(1,3)
      ddddx(2,3)=y*dz*(z2-x2-y2)*den
      ddddx(3,2)=ddddx(2,3)
c =====================================================================
c rotation to the equatorial reference system
c d2al/dy2=R * d2al/dx2 * Rt
      call prodmm(atmp,rot,ddadx)
      call prodmm(ddadx,atmp,rotinv)
c d2delta/dy2=R * d2delta/dx2 * Rt
      call prodmm(atmp,rot,ddddx)
      call prodmm(ddddx,atmp,rotinv)
      return
      end
c =====================================================================
*
*  ***************************************************************
*  *                                                             *
*  *                         A P P M A G                         *
*  *                                                             *
*  *     Calcolo della magnitudine apparente di un asteroide     *
*  *                                                             *
*  ***************************************************************
*
*
* INPUT:    H         -  Magnitudine assoluta
*           G         -  Parametro di slope
*           DS        -  Distanza dal Sole (AU)
*           DT        -  Distanza dalla Terra (AU)
*           BETA      -  Angolo di fase solare (rad): angolo tra il
*                        Sole e l'osservatore (Terra), visto dall'asteroide
*
      DOUBLE PRECISION FUNCTION appmag(h,g,ds,dt,beta)
      IMPLICIT NONE

      DOUBLE PRECISION h,g,ds,dt,beta
      DOUBLE PRECISION a1,a2,b1,b2
      DOUBLE PRECISION tb2,phi1,phi2
      SAVE a1,a2,b1,b2

* Costanti per il calcolo della magnitudine
      DATA a1,a2,b1,b2/3.33d0,1.87d0,0.63d0,1.22d0/

      tb2=TAN(beta/2.d0)
      phi1=EXP(-a1*tb2**b1)
      phi2=EXP(-a2*tb2**b2)
      appmag=5.d0*LOG10(ds*dt)+h-2.5d0*LOG10((1.d0-g)*phi1+g*phi2)
      RETURN
      END
c =====================================================================
c ALFDEL2
c =====================================================================
c Computation of alpha, delta, alphadot, deltadot  and their derivatives
c =====================================================================
c
c Input
c    tk: epoch time
c    tauj: observation time
c    ioj: station code
c    east: equinoctal orbital elements vector (a,h,k,p,q,lambda) 
c                    at time tk
c    ider: flag for derivatives options:
c        0 no derivatives
c        1 only first deriv. of alpha,delta w.r. to (east,tk,tauj)
c    twobo: logical flag for 2-body approximation; if .true., 2-body
c        approximation (but Earth position comes from JPL ephem); if .false.,
c        all orbit propagations are full n-body
c Output
c   alj,dej: alpha,delta  computed at time tauj
c   adot,ddot their time derivatives
c   dade,ddde:  matrices of first derivatives (if required)
c   ddade,dddde: matrices of second deriv. (if required, only in 2-body approx)
c
c ==============INTERFACE============================================
      SUBROUTINE alfdel2 (east,tk,tauj,iocj,obs,dobde,ider,twobo)
      implicit none
c ==============INPUT==========================
c flag to control two-body approximation: if .false., full n-body
      logical twobo
c times: epoch time for asteroid elements, observation time (MJD)
      double precision tk,tauj
c asteroid equinoctal elements 
      double precision east(6)
c cartesian coordinates of the Earth
      double precision xea(6)
c observatory code
      integer iocj
c flag to control computation of derivatives
      integer ider
c ============OUTPUT==============================
c observations: alpha (right ascension) delta (declination), in RAD
c alphadot, deltadot in RAD/day
      double precision obs(4)
c partial derivatives of obs, w.r. to asteroid coordinates
      double precision dobde(4,6)
c =============END INTERFACE=========================================
c asteroid cartesian coordinates
      double precision xast(6)
c first partial derivatives of alpha, delta, adot, ddotw.r. to ast. coordinates
      double precision dobdx(4,6)
c first, second derivatives of cartesian coordinates with respect to elements
      double precision dxde(6,6),ddxde(3,6,6)
c *************************************************************************
c****************
c   static memory not required
c****************
      IF(ider.gt.1.or.ider.lt.0)THEN
         WRITE(*,*)'alfdel2; ider=',ider,' not understood'
         STOP
      ENDIF
c Orbit propagation:
      if(twobo)then
c 2 body version
         call propa2 (tk,east,tauj,xast,xea,ider,dxde,ddxde)
      else
c full n-body numerical integration
         call propag (tk,east,tauj,xast,xea,ider,dxde,ddxde)
      endif
c Computation of observations
      call ossdif2(xast,xea,tauj,iocj,obs,ider,dobdx)
      if(ider.lt.1)return
c derivatives with respect to equinoctal elements
      CALL mulmat(dobdx,4,6,dxde,6,6,dobde)
      return
      end
c =====================================================================
c OSSDIF2
c =====================================================================
c Corrections to observations
c =====================================================================
c Input
c xast asteroid cartesian coordinates at time tauj
c xea Earth cartesian coordinates at time tauj
c     both in the ecliptic system
c tauj observation time
c ioc station code
c ider flag for derivatives options:
c       =0 no derivatives
c       <2 only first deriv. of $\alpha,\delta$ w.r. 
c to positions
c       great equal 2 also second partial derivatives 
c (approximation with 2-body case)
c
c Output
c obs=alpha,delta,adot,ddot computed at time tauj (in the equatorial system)
c (if required)
c dobdx matrix of first derivatives w. r. to ecliptic positions 
c (if required)
c =====================================================================
      SUBROUTINE ossdif2(xast,xea,tauj,ioc,obs,ider,dobdx)
      implicit none
c =====================================================================
      include 'trig.h'
      include 'model.h'
c =====================================================================
c observations alpha, delta (equatorial, radians), aodt,ddot and time
      double precision obs(4),tauj
c cartesian coordinates, ecliptic, of asteroid, of earth, vector difference
c WARNING: even if the velocities are not always propagated, they are available
c at the time of the observations
      double precision xast(6),xea(6),d(6)
c partials of equatorial alpha,delta,adot,ddot w.r.t. cartesian ecliptic coord. 
      double precision dobdx(4,6)
c topocentric position of the observatory, obs. code
      double precision xo(3),vo(3)
      integer idst,ioc
c elongation,distance to Earth, distance to Sun (to compute magnitude)
      include 'phase.h'
c phase and elongation cosine
      double precision cospha,coselo
c control on no. derivatives, bits in mantissa
      integer ider,n
c rounding off, auxiliary var.
      double precision errm,dz
c real functions
      double precision roff,vsize,prscal
c integer loop index
      integer i
c rotation matrices, rotated vectors
      double precision rot(3,3),rotinv(3,3),deq(6),tmp(3),ddd(3)
c scalar variables to preserve old code
      double precision alpha,delta
c partials of equatorial alpha, delta w.r. to cartesian ecliptic coord. 
      double precision dadx(3),dddx(3),ddadx(3,3),ddddx(3,3)
c aux. var for computation of second derivatives
      double precision x,y,z
      double precision den,x2,y2,z2,x4,y4
      double precision x2y2,x2z2,y2z2
c =====================================================================
      integer lflag
*********************************
* static memory allocation only for:
      save lflag,errm, rot,rotinv
********************************* 
      data lflag/0/
      if(lflag.eq.0)then
         errm=roff(n)
* Change of reference system EQUM00 ---> ECLM00
c Rotation matrix from the reference system in which the
c observations are given (mean equatorial J2000) to the reference 
c system in which the orbit is computed
         call rotpn(rot,'EQUM','J2000',0.d0,'ECLM','J2000',0.d0)
         call transp(rot,3,3,rotinv)
         lflag=1
      endif
c =====================================================================
c Difference vector
      call vdiff (xast,xea,d)
      call vdiff(xast(4),xea(4),d(4))
c =====================================================================
c Displacement of the station with respect to the center of the Earth
      if(istat.gt.0)then
         idst=ioc
         call pvobs(tauj,idst,xo,vo)
         call vdiff(d,xo,d)
         call vdiff(d(4),vo,d(4))
      endif
c =====================================================================
c Aberration (only time delay)
      if(iaber.gt.0)then
         call aber1(d,xast(4),d)
      endif
c =====================================================================
c Computation of solar distance, earth distance, phase, elongation
      dsun=vsize(xast)
      dis=vsize(d)
      cospha=prscal(d,xast)/(dis*dsun)
      pha=acos(cospha)
      coselo=-prscal(d,xea)/(dis*vsize(xea))
      elo=acos(coselo)
c =====================================================================
c rotation to the equatorial reference system
      call prodmv(deq,rotinv,d)
      call prodmv(deq(4),rotinv,d(4))
      call vcopy(6,deq,d)
c =====================================================================
c galactic latitude
      gallat=pig/2d0-acos((d(1)*gax+d(2)*gay+d(3)*gaz)/dis)
c =====================================================================
c Computation of observation: right ascension (radians)
      dz=d(1)**2+d(2)**2
      if (dz.le.errm) then
          alpha=0.d0
      else
          alpha=atan2(d(2),d(1))
          if (alpha.lt.0.d0) then
              alpha=alpha+dpig
          endif
      endif
c Computation of observation: declination (radians)
      delta=asin(d(3)/dis)
c =====================================================================
c Computation of first derivatives of $\alpha$ and $\delta$ w.r. to positions 
c (if required): we derived eq. (2.20)
      dadx(1)=-d(2)/dz
      dadx(2)=d(1)/dz
      dadx(3)=0.d0
      dddx(1)=-d(3)*(d(1)/(sqrt(dz)*dis**2))
      dddx(2)=-d(3)*(d(2)/(sqrt(dz)*dis**2))
      dddx(3)=sqrt(dz)/dis**2
c =====================================================================
c Apparent motion: 
      adot=prscal(dadx,d(4))
      ddot=prscal(dddx,d(4))
c store into obs vector
      obs(1)=alpha
      obs(2)=delta
      obs(3)=adot
      obs(4)=ddot
c check if observation partials are required
      if(ider.eq.0)RETURN
c =====================================================================
c partials of alpha, delta have already been computed
c partials of adot,ddot with respect to velocities are the same
c rotation to the equatorial reference system
      call prodmv(tmp,rot,dadx)
      do  i=1,3
        dobdx(1,i)=tmp(i)
        dobdx(1,i+3)=0.d0
        dobdx(3,i+3)=tmp(i)
      enddo
      call prodmv(tmp,rot,dddx)
      do i=1,3
        dobdx(2,i)=tmp(i)
        dobdx(2,i+3)=0.d0
        dobdx(4,i+3)=tmp(i)
      enddo
c =====================================================================
c partials of adot,ddot with respect to positions require the 
c second derivatives of alpha, delta with respect to the
c equatorial reference system
c =====================================================================
c Computation of second derivatives of $\alpha$ w.r. to positions
      ddadx(1,1)=2.d0*d(1)*d(2)/dz**2
      ddadx(1,2)=(d(2)**2-d(1)**2)/dz**2
      ddadx(2,1)=ddadx(1,2)
      ddadx(2,2)=-ddadx(1,1)
      ddadx(3,1)=0.d0
      ddadx(3,2)=0.d0
      ddadx(3,3)=0.d0
      ddadx(2,3)=0.d0
      ddadx(1,3)=0.d0
c Computation of second derivatives of $\delta$ w.r. to positions
      den=1.d0/(dis**2*dz*sqrt(dz))
      x=d(1)
      y=d(2)
      z=d(3)
*
      x2=x*x
      y2=y*y
      z2=z*z
      x4=x2*x2
      y4=y2*y2
*
      x2y2=x2*y2
      x2z2=x2*z2
      y2z2=y2*z2
*
      ddddx(1,1)=z*(2.d0*x4+x2y2-y2z2-y4)*den
      ddddx(2,2)=z*(2.d0*y4+x2y2-x2z2-x4)*den
      ddddx(1,2)=x*y*z*(z2+3.d0*x2+3.d0*y2)*den
      ddddx(2,1)=ddddx(1,2)
      ddddx(3,3)=-2.d0*z*dz**2*den
      ddddx(1,3)=x*dz*(z2-x2-y2)*den
      ddddx(3,1)=ddddx(1,3)
      ddddx(2,3)=y*dz*(z2-x2-y2)*den
      ddddx(3,2)=ddddx(2,3)
c =======================================================
c chain rule for derivatives of adot
      call prodmv(ddd,ddadx,d(4))
c =======================================================
c rotation to the equatorial reference system
      call prodmv(tmp,rot,ddd)
      do  i=1,3
        dobdx(3,i)=tmp(i)
      enddo
c =======================================================
c chain rule for derivatives of adot
      call prodmv(ddd,ddddx,d(4))
c =======================================================
c rotation to the equatorial reference system
      call prodmv(tmp,rot,ddd)
      do i=1,3
        dobdx(4,i)=tmp(i)
      enddo
c
      return
      end
* Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: February 24, 1997
*
*  ***************************************************************
*  *                                                             *
*  *                          A B E R 1                          *
*  *                                                             *
*  *          Correzione approssimata per aberrazione            *
*  *                  stellare e/o planetaria                    *
*  *                                                             *
*  ***************************************************************
*
*
* INPUT:    XREL(3)   -  Posizione relativa vera del corpo osservato (UA)
*           VREL(3)   -  Velocita` relativa (UA/d)
*
* OUTPUT:   XCOR(3)   -  Posizione relativa apparente (tenendo conto della
*                        aberrazione)
*
* NOTA: in generale per ottenere la correzione completa (comprendente le
*       cosiddette aberrazioni "stellare" + "planetaria") bisogna che VREL
*       sia la velocita` relativa del corpo osservato rispetto all'osserva-
*       tore:
*                 VREL  =   V(pianeta) - V(osservatore)
*
*       Se si vuole ottenere solo la correzione per l'aberrazione "stellare"
*       bisogna porre VREL =  - V(osservatore).
*
*
      SUBROUTINE aber1(xrel,vrel,xcor)
      implicit double precision (a-h,o-z)
      dimension xrel(3),vrel(3),xcor(3)
* Velocita` della luce in UA/d
      data c/173.14463272d0/
      ro=vsize(xrel)
c  effetto di ritardo
      dt=ro/c
      do 1 j=1,3
 1    xcor(j)=xrel(j)-dt*vrel(j)
      return
      end
c ==================================
c R R D O T
c
c radar observations
      SUBROUTINE rrdot (east,iobs,t0,tr,ioc,r,v,drde,dvde,ider,twobo)
      implicit none
c =============INPUT====================
c asteroid equinoctal elements 
      double precision east(6)
c observation type (used for surface bounce correction)
      integer iobs
c times: epoch time for asteroid elements, observation time (MJD)
      double precision t0,tr
c observatory code
      integer ioc
c flag to control computation of derivatives
      integer ider
c flag to control two-body approximation: if .false., full n-body
      logical twobo
c ============OUTPUT====================
c observations: range and range rate in AU, AU/day
      double precision r,v
c partial derivatives of range and range w.r. to asteroid coordinates
      double precision drde(6),dvde(6)
c =============END INTERFACE=========================================
c cartesian coordinates of the Earth, of the asteroid, id. at receive time
      double precision xea(6),xast(6),xastr(6)
c first partial derivatives of r, rdot, w.r. to ast. coordinates, vel.
      double precision drdx(6),dvdx(6)
c asteroid radius (for surface bounce correction)
      double precision rb
c station codes
      INTEGER iotr,iore
c station positions: geocentric, heliocentric
      double precision xre(3),yre(3),xtr(3),ytr(3)
      double precision rre(3),vre(3),rtr(3),vtr(3)
c difference vector, distance, velocity difference, size
      double precision rhorv(3),rhor,rhordv(3),rhord
      double precision rhotv(3),rhot,rhotdv(3),rhotd
      double precision vsize,prscal,prscag
c solar system barycentric velocities (for 1/c^2 Doppler level)
      double precision vressb(3),vastssb(3),vtrssb(3)
      double precision vressb2,vtrssb2
c down leg time, bounce time
      double precision taud,taudold,tb,tbold
c up leg time, transmit time
      double precision tauu,tt,ttold,tauuold
c delta tau correction function
      double precision deltau
c speed of light
      INCLUDE 'vlight.h'
c radius of asteroid
      include 'radius.h'
c iteration index, conrol and flag for the position/velocity subroutine
      INTEGER i,itmax,ifla
      PARAMETER (itmax=10)
c control on convergence set at 0.05 microseconds (Yeomans et al. AJ 1992)
      DOUBLE PRECISION ept
c     PARAMETER (ept=6.d-13)
c correction 9/1/2002: control on taud, not on time in MJD
      PARAMETER (ept=1.d-16)
      DOUBLE PRECISION tbdif(itmax),taudif(itmax)
c time scale
      double precision tretdb,temp,tdiffr,tdifft
c 2-body step
      DOUBLE PRECISION eqast(6),enne
      INCLUDE 'sunmass.h'
c for AU
      INCLUDE 'jplhdr.h'
c first derivatives of cartesian coordinates with respect to elements
      double precision dxde(6,6)
c second derivatives (dummy),planet coords (dummy)
      double precision ddxde(3,6,6),xdummy(6)
c auxiliary scalars for Doppler GR & tropospheric corrections
      double precision rgeo,rgeod,deldoptr,deldopre,rast,rastd
      double precision rsta,rstad,scal1,scal2,scal3,levelc
      double precision dtaud,dtauu
c =================================================
c two body approx must not be used
      IF(twobo)STOP'radar and two body approximations incompatible'
c deal with surface vs. mass center return:
      if(iobs-2000.ge.100)then
         rb=radius
         if(rb.le.0)then
           write(*,*)  '**** rrdot: internal error', radius
           stop
         endif
      else
c no correction is needed
         rb=0
      endif
c ======================================
c Displacement of the stations with respect to the center of the Earth
c find codes of two observatories
      iotr=ioc/10000
      iore=ioc-iotr*10000
      ifla=1
c compute position of receiver at time tr
      call pvobs2(tr,iore,xre,yre)
c =======================================
c get receive time in TDB
      call times(tr+2400000.5,temp,tdiffr)
      tretdb=tr+tdiffr/86400.d0
c Compute down leg time
c Orbit propagation at time tr
      call propag(t0,east,tretdb,xastr,xea,ider,dxde,ddxde)
c ... new initial conditions are xastr
c receive station position at time tr 
      CALL vsumg(3,xea,xre,rre)
      CALL vsumg(3,xea(4),yre,vre)
c initial guess is bounce time(position) equal to receive time(position)
      tb=tretdb
      taud=0.d0
      do i=1,6
         xast(i)=xastr(i)
      enddo
c Loop on down leg
      DO i=1,itmax
         CALL vdiff(xast,rre,rhorv)
         rhor=vsize(rhorv)      
         dtaud=deltau(xast,xre,rhorv,rre)
c correction
         taudold=taud
         taud=(rhor-rb)/vlight + dtaud
         tbold=tb
         tb=tretdb-taud
c convergence control
         tbdif(i)=tb-tbold
         taudif(i)=taud-taudold
c         WRITE(*,*)tb,tb-tbold,taud,taud-taudold
         IF(abs(taud-taudold).lt.ept) GOTO 9
c Orbit propagation at time tb
         CALL coocha(xastr,'CAR',gms,eqast,'EQU',enne)
         CALL prop2b(tretdb,eqast,tb,xast,gms,0,dxde,ddxde)
      ENDDO
c too many iterations
      WRITE(*,*)' slow conv. on down leg time ',(tbdif(i),i=1,itmax),
     +  (taudif(i),i=1,itmax)
c compute relative velocity between asteroid and receiver
 9    CONTINUE
      CALL vdiff(xast(4),vre,rhordv)
      rhord=prscal(rhorv,rhordv)/rhor
c solar velocity at the receive and bounce epochs to get the asteroid
c barycentric velocity (vastssb) and the receiver barycentric velocity 
c (vressb)
      ifla=2
c - receive
      CALL earcar(tretdb,xea,ifla)
      CALL vsumg(3,vre,xea(4),vressb)
      vressb2=vressb(1)*vressb(1)+vressb(2)*vressb(2)+
     +        vressb(3)*vressb(3)
c - bounce
      CALL earcar(tb,xea,ifla)
      CALL vsumg(3,xast(4),xea(4),vastssb)
c =======================================
c compute upleg time
c fist guess is up leg time = down leg time
      ifla=1
      tt=tb-taud
      tauu=taud
c Loop on upleg time
      DO i=1,itmax
c compute transmitter position at estimated transmit time tt
c call the Earth-rotation model in TDT
         call times(tt+2400000.5,temp,tdifft)
         CALL earcar(tt-(tdifft/86400.d0),xea,ifla)
         CALL pvobs2(tt,iotr,xtr,ytr)
         CALL vsumg(3,xea,xtr,rtr)
c upleg time
         CALL vdiff(xast,rtr,rhotv)
         rhot=vsize(rhotv)
         dtauu=deltau(xast,xtr,rhotv,rtr)
         tauuold=tauu
         tauu=(rhot-rb)/vlight + dtauu
         ttold=tt
         tt=tb-tauu
c convergence control
c        WRITE(*,*)tt,tt-ttold,tauu, tauu-tauuold
         IF(abs(tauu-tauuold).lt.ept) GOTO 19
      ENDDO
c too many iterations
      WRITE(*,*)' slow convergence up leg time ',tt-ttold, tauu-tauuold
 19   CONTINUE
c compute relative velocity between asteroid and transmitter
      CALL vsumg(3,xea(4),ytr,vtr)
      CALL vdiff(xast(4),vtr,rhotdv)
      rhotd=prscal(rhotv,rhotdv)/rhot
c ==========================================================
c compute distance
      r=0.5d0*(tauu+taud+(tdifft-tdiffr)/86400.d0)*vlight
c compute relative frequency shift (up to 1/c^3 level);
c solar velocity at the bounce epochs and the asteroid barycentric
c velocity (vastssb)
      ifla=2
      CALL earcar(tt,xea,ifla)
      CALL vsumg(3,vtr,xea(4),vtrssb)
      vtrssb2=vtrssb(1)*vtrssb(1)+vtrssb(2)*vtrssb(2)+
     +        vtrssb(3)*vtrssb(3)
      levelc=rhotd+rhord
      scal1=prscal(rhotv,vtrssb)/rhot/vlight
      scal2=prscal(rhorv,vastssb)/rhor/vlight
      scal3=0.5d0*(vtrssb2-vressb2)+
     +      gms*((1.d0/vsize(rtr))-(1.d0/vsize(rre)))
      v=0.5d0*(levelc+rhotd*scal1*(1.d0+scal1)
     +               -rhord*scal2*(1.d0-scal2)
     +               -(rhotd*rhord*(1.d0+scal1-scal2)
     +               -scal3*(1.d0-(levelc/vlight)))/vlight)
c - get GR and tropospheric corrections to Doppler
c -- GR stuff
      rast=vsize(xast)
      rastd=prscal(xast,xast(4))/rast
c a) upleg piece
      rsta=vsize(rtr)
      rstad=prscal(rtr,vtr)/rsta
      call deldop1(rast,rastd,rhot,rhotd,rsta,rstad,deldoptr)
c b) downleg piece
      rsta=vsize(rre)
      rstad=prscal(rre,vre)/rsta
      call deldop1(rast,rastd,rhor,rhord,rsta,rstad,deldopre)
      v=v-0.5d0*(deldoptr+deldopre)
c -- troposheric stuff
c a) at transmit passage -->
      rgeo=vsize(xtr)
      rgeod=prscal(xtr,ytr)/rgeo
      call deldop2(xast,xtr,ytr,rgeo,rgeod,rhotv,rhotdv,rhot,rhotd,
     . deldoptr)
c b) at receive passage <--
      rgeo=vsize(xre)
      rgeod=prscal(xre,yre)/rgeo
      call deldop2(xast,xre,yre,rgeo,rgeod,rhorv,rhordv,rhor,rhord,
     . deldopre)
      v=v-0.5d0*(deldoptr+deldopre)
c rem interplanetary environment effects (e^- plasma) neglected
c ==========================================================
c Derivatives
      IF(ider.eq.0)RETURN
c derivs of r,rdot wrt cartesian
      do i=1,3
c d(range)/d(r)
         drdx(i)=(rhotv(i)/rhot +rhorv(i)/rhor)/2.d0
c d(range)/d(v) = 0
         drdx(i+3)=0
c d(range-rate)/d(r)
         dvdx(i)=-((rhotd*rhotv(i)/rhot-rhotdv(i))/rhot + 
     +             (rhord*rhorv(i)/rhor-rhordv(i))/rhor)/2.d0
c d(range-rate)/d(v)= c * d(range)/d(r)
         dvdx(i+3)=drdx(i)
      enddo
c derivs of cartesian with respect to equinoctal elements
      do i=1,6
        drde(i)=prscag(6,drdx,dxde(1,i))
        dvde(i)=prscag(6,dvdx,dxde(1,i))
      enddo
      RETURN
      END

c ======================================================================
c DELTAU - "small" corrections to radar time of flight
c ======================================================================
c     xast - asteroid position, heliocentric
c     xsta - station position, relative to Earth center
c     rho - asteroid position, relative to station
c     r - station position, heliocentric
      DOUBLE PRECISION FUNCTION deltau(xast,xsta,rho,r)
      implicit none
      double precision xast(3),xsta(3),rho(3),r(3)
      include 'sunmass.h'
      include 'vlight.h'
      double precision vsize,prscag
      double precision rsta,e,p,q,cosz,cotz,deltau1,deltau2
c     double precision sinha,sin2ha,fghz,ampli,finte1,fun1,fun2
c     double precision aprim,bprim,deltau3,omeg1,alpha
c ================================
c Relativistic delay
      e=vsize(r)
      p=vsize(xast)
      q=vsize(rho)
      deltau1=2d0*gms*log(abs((e+p+q)/(e+p-q)))/vlight**3
c Earth ionospheric/tropospheric delay
c ref. EM Standish, A&A 233, 252 (1990)
      rsta=vsize(xsta)
      if(rsta.lt.1d-12)then
         write(*,*)'deltau: radar station at geocenter!'
         cosz=1d0-1d-8
      else
         cosz=prscag(3,xsta,rho)/rsta/q
      endif
      if(cosz.eq.1d0)cosz=cosz-1d-8
      cotz=cosz/sqrt(1d0-cosz**2)
      deltau2=(7d-9/86400d0)/(cosz+1.4d-3/(4.5d-2+cotz))
c Interplanetary medium propagation effect
c ref. EM Standish, A&A 233, 252 (1990) with constants of DE118
c rem. a more precise model might be needed here; e.g.
c      Muhleman & Anderson, ApJ 247, 1093 (1981) or newer
c      alpha=q/e
c      cosha=(r(1)*rho(1)+r(2)*rho(2)+r(3)*rho(3))/e/q
c      sin2ha=1.d0-cosha*cosha
c      sinha=dsqrt(sin2ha)
c X- or S-band;
cc !!! Information about the frequency is not passed here at the
cc     moment; it should be decided manually !!!
cc      fghz=2.38d0
c      fghz=8.51d0
c      ampli=(2.01094d-8)/fghz/fghz/e/86400.d0
c      aprim=1.237265d-6
c      bprim=9.524021d0
c      finte1=(datan((alpha+cosha)/sinha)-datan(cosha/sinha))/sinha
c      fun1=bprim*finte1
c      omeg1=1.d0+alpha*(2.d0*cosh+alpha)
c      fun2=aprim*(0.25d0*((alpha+cosha)*((1.d0/omeg1)+(1.5d0/sin2ha))
c     .           /omeg1-cosha*(1.d0+(1.5d0/sin2ha)))/sin2h
c     .           +0.375d0*finte1/sin2ha/sin2ha)/(e**4)
c      deltau3=ampli*(fun1+fun2)
c Add 'em up
cc      deltau=deltau1+deltau2+deltau3
      deltau=deltau1+deltau2
      return
      end

c ======================================================================
c DELDOP1 - "small" corrections to radar-rate measurements
c ======================================================================
      SUBROUTINE deldop1(p,pdot,q,qdot,e,edot,deldop)
      implicit none
      double precision p,pdot,q,qdot,e,edot,deldop
      include 'sunmass.h'
      include 'vlight.h'
      double precision brac1,brac2
c ================================
c relativistic range-rate correction
      brac1=q*(edot+pdot)-qdot*(e+p)
      brac2=((e+p)**2)-(q**2)
      deldop=4.d0*gms*brac1/brac2/(vlight**2)
      return
      end

c ======================================================================
c DELDOP2 - "small" corrections to radar-rate measurements
c ======================================================================
c     xast(6) - asteroid r & v heliocentric
c     xtr(3),ytr(3) - station r & v relative to Earth center
c     rsta,drsta - |xtr| & d|xtr|/dt
c     rhov(3),drhov(3) - asteroid r & v relative to station
c     rho,drho - |rhov| & d|rhov|/dt
      SUBROUTINE deldop2(xast,xsta,vsta,rsta,drsta,rhov,drhov,rho,drho,
     . deldop)
      implicit none
      double precision xast(6),xsta(3),vsta(3),rhov(3),drhov(3)
      double precision rsta,drsta,rho,drho,deldop
      include 'vlight.h'
      double precision cosz,sinz,cotz,dcoszdt,phiz
      double precision brac,brac1,brac2,scal1,scal2
c ================================
c rate of change of the Earth ionospheric/tropospheric ~ Doppler shift
      if(rsta.lt.1d-12)then
         write(*,*)'deltau: radar station at geocenter!'
         cosz=1d0-1d-8
      else
         cosz=(rhov(1)*xsta(1)+rhov(2)*xsta(2)+rhov(3)*xsta(3))/rsta/rho
      endif
      sinz=sqrt(1.d0-cosz**2)
      cotz=cosz/sinz
      brac=0.045d0+cotz
      brac1=1.d0-(0.0014d0/brac/brac/(sinz**3))
      brac2=cosz+(0.0014d0/brac)
      phiz=-vlight*(7.d-9/86400.d0)*brac1/(brac2**2)
      scal1=drhov(1)*xsta(1)+drhov(2)*xsta(2)+drhov(3)*xsta(3)
      scal2=rhov(1)*vsta(1)+rhov(2)*vsta(2)+rhov(3)*vsta(3)
      dcoszdt=((scal1+scal2)/rho/rsta)-cosz*((drho/rho)+(drsta/rsta))
      deldop=phiz*dcoszdt
      return
      end








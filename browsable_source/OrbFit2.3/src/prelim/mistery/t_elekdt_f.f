*
* Computation of the force function for an (N+1) + M body problem
* (Sun + N planets + M massless asteroids) in heliocentric
* coordinates (keplerian term + direct and indirect perturbations)
*
      subroutine force(x,v,tm,f)
      implicit double precision (a-h,o-z)
* Max no. of planets
      parameter (nplx=10)
* Max no. of asteroids
      parameter (nasx=10)
      parameter (ntotx=nplx+nasx)
      dimension recplm(*)
      dimension x(3,ntotx),v(3,ntotx),f(3,ntotx)
      dimension pm(nplx),gm1(ntotx),gm(nplx)
      dimension r(ntotx),r3(ntotx),dr(3)
      dimension pind(3,nplx),pintot(3)
      logical start
      save
      data start/.true./
      if(start)then
          write(*,*)' **** ERROR IN SUBROUTINE FORCE:'
          write(*,*)' **** CALLING FORCE WITHOUT INIZIALIZATION'
          stop
      end if
*
* Distance from the Sun and indirect perturbations
*
      do 3 k=1,ntot
      u=x(1,k)**2+x(2,k)**2+x(3,k)**2
      r(k)=dsqrt(u)
 3    r3(k)=r(k)**3
      do 33 k=1,npl
      do 33 i=1,3
 33   pind(i,k)=-gm(k)*x(i,k)/r3(k)
*
* Sum of indirect perturbations
*
      if(nast.eq.0)goto 6
      do 5 i=1,3
      u=0.d0
      do 4 k=1,npl
 4    u=u+pind(i,k)
 5    pintot(i)=u
 6    continue
*
* Inizializing total force
*
      do 7 k=1,ntot
      do 7 i=1,3
 7    f(i,k)=0.d0
*
* Mutual perturbations between planets
*
      do 10 k1=1,npl
      do 10 k2=1,npl
      if(k1.eq.k2)goto 10
      do 8 i=1,3
 8    dr(i)=x(i,k2)-x(i,k1)
      dr3=dr(1)**2+dr(2)**2+dr(3)**2
      dr3=dsqrt(dr3)**3
      do 9 i=1,3
 9    f(i,k1)=f(i,k1)+gm(k2)*dr(i)/dr3+pind(i,k2)
 10   continue
*
* Planetary perturbations on asteroids
*
      if(nast.eq.0)goto 14
      do 13 k1=npl+1,ntot
      do 12 k2=1,npl
      do 11 i=1,3
 11   dr(i)=x(i,k2)-x(i,k1)
      dr3=dr(1)**2+dr(2)**2+dr(3)**2
      dr3=dsqrt(dr3)**3
      do 12 i=1,3
 12   f(i,k1)=f(i,k1)+gm(k2)*dr(i)/dr3
      do 13 i=1,3
 13   f(i,k1)=f(i,k1)+pintot(i)
 14   continue
*
* Central term
*
      do 15 k=1,ntot
      do 15 i=1,3
 15   f(i,k)=-gm1(k)*x(i,k)/r3(k)+f(i,k)
      return
*
* Entry FORCIN: Initialization (definition of constants)
*
* N       =  no. of planets
* M       =  no. of asteroids
* CENTM   =  central mass (Sun)
* RECPLM  =  Reciprocal planetary masses
*
      entry forcin(n,m,centm,recplm)
      if(n.gt.nplx)stop' **** ERROR IN FORCIN: N > NPLX ****'
      if(m.gt.nasx)stop' **** ERROR IN FORCIN: M > NASX ****'
* Gauss gravitational constant
      gk=0.01720209895d0
      g=gk*gk
* No. of objects
      npl=n
      nast=m
      ntot=npl+nast
* Masses and related constants
      gm0=g*centm
      do 1 k=1,npl
      pm(k)=1.d0/recplm(k)
      gm(k)=g*pm(k)
 1    gm1(k)=g*(centm+pm(k))
      do 2 k=npl+1,ntot
 2    gm1(k)=gm0
      start=.false.
      return
      end

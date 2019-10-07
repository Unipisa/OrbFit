*
* Numerical integration of the Sun-Jupiter-Saturn 3-body
* problem to test the filter routines
*
      program numint
      implicit double precision (a-h,o-z)
      parameter (ntotx=20)
      parameter (nfx=1500)
      dimension x(3,ntotx),v(3,ntotx),elem(6,ntotx),gm(ntotx)
      dimension rm(ntotx)
      dimension datin(2),datout(2)
      dimension radin(3),radinp(3),radout(3)
      dimension nrin(3),nrinp(3),nrout(3)
      dimension h(nfx)
      character type*4
      logical out1,out2
      dpi=8.d0*atan(1.d0)
* Gauss' gravitational constant
      gk=0.01720209895d0
      g=gk*gk
* File opening
      open(1, file='jupsat.in',status='old')
      open(2, file='jupsat.out',status='unknown')
      open(3, file='jupsat.fil',status='unknown')
      open(10,file='filter.d24',status='old')
* Input of initial conditions:
* - number of planets and asteroids
      read(1,*)n,m
      nm=n+m
* - number of integration steps
      read(1,*)nstep
* - integration output interval
      read(1,*)dt
* - controls for RA15 routine
      read(1,*)xl,ll
      if(nm.gt.ntotx)stop' **** Error in numint: n+m > ntotx ****'
* - Sun mass (usually, =1)
      read(1,*)sunmas
* - planets' reciprocal masses
      do 7 k=1,n
      read(1,*)rm(k)
 7    gm(k)=g*(sunmas+1.d0/rm(k))
* - initial conditions of planets and asteroids (positions and
*   velocities)
      do 1 k=1,nm
 1    read(1,*)(x(i,k),i=1,3),(v(i,k),i=1,3)
* Setup of constants in subroutine force
      call forcin(n,m,sunmas,rm)
* Input of filter coefficients
      call rdfil(10,h,nfil,idec,nfx)
* Filter parameters
      ndat1=2
      ndat2=3
      nskip1=0
      nskip2=0
      hlen=(nfil-1)*dt/2.d0
      eps=dt*1.d-12
* Parameters for subroutine ra15
      nor=15
      neq=3*nm
* Here starts the integration loop
      do 10 kk=1,nstep
      dt1=dt
* Numerical integration
      call ra15 (x,v,dt1,xl,ll,neq,-2,nor)
      time=dt*kk
* Computation of keplerian elements
      do 3 k=1,nm
      call ccek1(elem(1,k),type,x(1,k),v(1,k),gm(k))
 3    continue
* Preparation of filter input.
* This part of the code is strictly dependent upon the particular
* problem. Usually you would probably pass to the filter all the
* orbital elements; here we are interested only in the semimajor
* axes, mean anomalies and the resonant argument of the 5:2
* mean motion resonance between Jupiter and Saturn.
      datin(1)=elem(1,1)
      datin(2)=elem(1,2)
      radin(1)=elem(4,1)
      radin(2)=elem(4,2)
      radin(3)=-2*radin(1)+5*radin(2)
* Count of the number of revolutions
      if(kk.eq.1)then
          do 4 j=1,ndat2
 4        nrin(j)=0
      else
          do 5 j=1,ndat2
          nrev=nint((radin(j)-radinp(j))/dpi)
 5        nrin(j)=nrinp(j)-nrev
      end if
      do 6 j=1,ndat2
      radinp(j)=radin(j)
 6    nrinp(j)=nrin(j)
* Output of sampled data
      years=time/365.25d0
      write(2,100)years,(datin(j),j=1,ndat1),
     +                  (radin(j)+dpi*nrin(j),j=1,ndat2)
* Filtering
      call filter(time,datin,ndat1,h,nfil,hlen,idec,nskip1,
     +            tout1,datout,out1)
      call filtan(time,radin,nrin,ndat2,h,nfil,hlen,idec,nskip2,
     +            tout2,radout,nrout,out2)
* Check on filter output
      if(out1.neqv.out2)stop' **** Filters are out of phase (1) ****'
      if(abs(tout1-tout2).gt.eps)
     +         stop' **** Filters are out of phase (2) ****'
* Output of filtered results
      if(out1)then
          years=tout1/365.25d0
          write(3,100)years,(datout(j),j=1,ndat1),
     +                      (radout(j)+dpi*nrout(j),j=1,ndat2)
      end if
* Here ends the integration loop
 10   continue
      stop
 100  format(f11.5,2f13.10,2f14.8,f14.10)
      end

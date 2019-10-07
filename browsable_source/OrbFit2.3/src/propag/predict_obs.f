c ===========MODULE predict_obs.f
c PUBLIC ROUTINES
c                 preobs
c                 preobc
c                 preobn 
c    Note: to be unified using optional arguments
c
c MODULE CONTAINS
c ROUTINES     
c                 linobs
c                 ellips
c                 elemov
c                 slinel
c                 graha  Note: these are required also by target_plane.f
c
c                 angupd(ang,vang,ng)
c                 outbc
c
c  HEADERS
c
c
c
c =====================================================================
c PREOBS-PREOBC
c =====================================================================
c PREOBS- version with angles only, no covariance
c =====================================================================
c  input:  coo   = coordinate type EQU, KEP, CAR, EQP
c          t0    = epoch time (MJD)
c          idsta = station code
c          t1    = prediction time (MJD)
c          east0 = orbital elements vector at time t0
c          h     = absolute magnitude
c          g     = opposition effect coefficient
c          iobs  = observation type
c  output: 
c          IF iobs=1000+x
c          x = right ascension (equatorial J2000), radians
c          y = declination (equatorial J2000), radians
c          hmagn = apparent magnitude, as predicted, from h and g given
c  WARNING: the magnitudes are often very poorly predictable
c
c          IF iobs=2000 + x
c          x = range in AU
c          y = range rate in AU/day
c                         (hmagn is dummy)
c          IF iobs=4000's
c          x = d alpha/dt (radians/day)
c          y = d delta/dt (radians/day)
c                         (hmagn is dummy)
c
c ============INTERFACE===================================================
      SUBROUTINE preobs(coo,t0,idsta,t1,east0,iobs,x,y,h,g,hmagn)
      implicit none
c elements and epoch times
      character*3 coo
      double precision east0(6),t0,t1
c observation type
      integer iobs
c magnitude input and output (only for iobs=1000's)
      double precision h,g,hmagn
c station code
      integer idsta
c observations
      double precision x,y
c ============END INTERFACE===============================================
c asteroid equinoctal elements and cartesian coordinates
      double precision eq(6),enne  
c partial derivatives of alpha, delta, w.r. to elements (not used)  
      double precision dxde(6),dyde(6)
c second derivatives of alpha, delta, w.r. to elements (not used)
      double precision ddxde(6,6),ddyde(6,6)
c vector observations
      double precision obs(4),dobde(4,6)
c functions
      double precision appmag
c control on derivatives
      integer ider
c elongation,distance to Earth, distance to Sun (to compute magnitude)
      include 'phase.h'
c ======== constant of gravitation ==============
      include 'sunmass.h'
c flag for 2-body approximation; must be .false. for full n-body computation
      logical twobo
      twobo=.false.
c****************
c   static memory not required
c****************
c =====================================================================
c coordinate change
      call coocha(east0,coo,gms,eq,'EQU',enne) 
c compute observation   
      ider=0
      IF(iobs/1000.eq.1)THEN   
c astrometric measurements of angles      
         call alfdel (eq,t0,t1,idsta,x,y,dxde,dyde,ider,twobo,
     +        ddxde,ddyde)
c compute apparent magnitude at time of observation
         hmagn=appmag(h,g,dsun,dis,pha)
      ELSEIF(iobs/1000.eq.2)THEN
c radar data
         call rrdot (eq,iobs,t0,t1,idsta,x,y,dxde,dyde,ider,twobo)
         hmagn=0.d0
      ELSEIf(iobs/1000.eq.4)THEN
c proper motion
         call alfdel2 (eq,t0,t1,idsta,obs,dobde,ider,twobo)
         x=obs(3)
         y=obs(4)
         hmagn=appmag(h,g,dsun,dis,pha)
      ELSE
         WRITE(*,*)'preobs: this we have not invented yet, iobs=',iobs
         STOP
      ENDIF
      return
      end
c =====================================================================
c PREOBC- version with covariance, linear theory
c =====================================================================
c  input:  coo   = coordinate type EQU, KEP, CAR, EQP
c          t0    = epoch time (MJD)
c          idsta = station code
c          t1    = prediction time (MJD)
c          east0 = orbital elements vector at time t0
c          h     = absolute magnitude
c          g     = opposition effect coefficient
c          iobs  = observation type
c  output: 
c          IF iobs=1000's
c          x = right ascension (equatorial J2000), radians
c          y = declination (equatorial J2000), radians
c          hmagn = apparent magnitude, as predicted, from h and g given
c  WARNING: the magnitudes are often very poorly predictable
c
c          IF iobs=2000's
c          x = range in AU
c          y = range rate in AU/day
c                         (hmagn is dummy)
c          IF iobs=4000's
c          x = d alpha/dt (radians/day)
c          y = d delta/dt (radians/day)
c                         (hmagn is dummy)
c
c  In the linear approximation, the ellipse of confidence has principal axes
c          along axes; the semiaxes lenghts are sig 
c
c ============INTERFACE===================================================
      SUBROUTINE preobc(coo,t0,idsta,t1,east0,h,g,gamm0,
     +    iobs,x,y,hmagn,gamad,sig,axes)
      implicit none
c ============= input ====================================================
c elements and epoch times, covraiance at t0
      character*3 coo
      double precision east0(6),t0,t1,gamm0(6,6)
c magnitude
      double precision h,g,hmagn
c station code
      integer idsta
c observation type
      integer iobs
c ============= output ===================================================
c observations
      double precision x,y
c covariance
      double precision gamad(2,2),axes(2,2),sig(2)
c ============END INTERFACE===============================================
c asteroid equinoctal elements and cartesian coordinates
      double precision eq(6),enne  
c partial derivatives of alpha, delta, w.r. to elements (by columns, by rows)
      double precision dxydet(6,2),dadde(2,6)
c second derivatives of alpha, delta, w.r. to elements (not used)
      double precision ddade(6,6),dddde(6,6)
c jacobian matrices of partial derivatives, eigenvalues, workspace
      double precision eigval(2),tmp(2,6),fv1(2),fv2(2)
c covariance, normal matrices
      DOUBLE PRECISION gameq(6,6)
c jacobian matrices of partial derivatives
      DOUBLE PRECISION deqde(6,6),deqdet(6,6),tmp6(6,6)
c vector observations
      double precision obs(4),dobde(4,6)
c control on derivatives
      integer ider
c error flag
      integer ierr
c loop indexes 
      integer i,j
c functions
      double precision appmag
c elongation,distance to Earth, distance to Sun (to compute magnitude)
      include 'phase.h'
c ======== constant of gravitation ==============
      include 'sunmass.h'
c flag for 2-body approximation; must be .false. for full n-body computation
      logical twobo
      twobo=.false.
c****************
c   static memory not required
c****************
c =====================================================================
c coordinate change
      call cooder(east0,coo,gms,eq,'EQU',enne,deqde)
c compute normal and covariance matrix for equinoctal elements
      CALL mulmat(deqde,6,6,gamm0,6,6,tmp6)
      CALL transp(deqde,6,6,deqdet)
      CALL mulmat(tmp6,6,6,deqdet,6,6,gameq)
c compute observation   
      ider=1
      IF(iobs/1000.eq.1)THEN   
c astrometric measurements of angles      
         call alfdel (eq,t0,t1,idsta,x,y,dxydet(1,1),dxydet(1,2),
     +        ider,twobo,ddade,dddde)
c compute apparent magnitude at time of observation
         hmagn=appmag(h,g,dsun,dis,pha)
      ELSEIF(iobs/1000.eq.2)THEN
c radar data
         call rrdot (eq,iobs,t0,t1,idsta,x,y,dxydet(1,1),dxydet(1,2),
     +        ider,twobo)
         hmagn=0.d0
      ELSEIf(iobs/1000.eq.4)THEN
c proper motion
         call alfdel2 (eq,t0,t1,idsta,obs,dobde,ider,twobo)
         x=obs(3)
         y=obs(4)
         DO i=1,2
           DO j=1,6
             dxydet(j,i)=dobde(i+2,j)
           ENDDO
         ENDDO
         hmagn=appmag(h,g,dsun,dis,pha)
      ELSE
         WRITE(*,*)'preobc2: this we have not invented yet, iobs=',iobs
         STOP
      ENDIF
c =====================================================================
c compute covariance of alpha, delta
      call transp(dxydet,6,2,dadde)
      call mulmat(dadde,2,6,gameq,6,6,tmp)
      call mulmat(tmp,2,6,dxydet,6,2,gamad)
c compute ellipse of confidence
c eigenvalues
      call rs(2,2,gamad,eigval,1,axes,fv1,fv2,ierr)
      do 1 i=1,2
        if(eigval(i).gt.0.d0)then
           sig(i)=sqrt(eigval(i))
        else
           write(*,*) 'non positive eigenvalue'
           sig(i)=0.d0
        endif
 1    continue

      return
      end
c =====================================================================
c PREOBN- version with confidence region, nonlinear theory
c =====================================================================
c  input:  coo   = coordinate type EQU, KEP, CAR
c          t0    = epoch time (MJD)
c          idsta = station code
c          t1    = prediction time (MJD)
c          east0 = orbital elements vector at time t0
c          h     = absolute magnitude
c          g     = opposition effect coefficient
c          gamm0 = covariance at time t0 for elements east0
c          c0    = normal matrix at time t0 (should be the inverse of gamm0)
c          sigma = level of the confidence boundary in RMS values
c          npo   = number of points in the boundary
c          ibv   = Type of depiction
c                       =0 for automatic selection
c                       =1 for confidence boundary
c                       =2 for line of maximum variation
c          inl   = handling of nonlinearity
c                       =0 for automatic selection
c                       =1 for linear ellipse
c                       =2 for 2-Body nonlinear propagation of covariance
c                       =3 for n-body nonlinear propagation of covariance
c          iob1  = observation type:
c                        1000's = RA/DEC, 
c                        2000's=R/RDOT, 
c
c  output: alpha = right ascension (equatorial J2000), radians
c          delta = declination (equatorial J2000), radians
c          hmagn = apparent magnitude, as predicted, from h and g given
c          gamad = covariance matrix of observations alpha, delta
c          sig   = sqrt(eigenvalues) of gamad
c          axes  = the eigenvectors of gamad are the columns of this matrix
c          npo1  = number of output dta points 9could be less than npo)
c          al(npo1),de(npo1) points on the confidence boundary  
c                   (difference with respect to best prediciton, radians)
c          elm(npo1) alternate elements for observation time
c
c  In the linear approximation, the ellipse of confidence has principal axes
c          along axes; the semiaxes lenghts are sig
c  In the nonlinear approximation, the boundary is a map of the
c          confidence ellipse in the elements space 
c
c  WARNING: the magnitudes are often very poorly predictable
c ============INTERFACE===================================================
      SUBROUTINE preobn(coo,t0,idsta,t1,east0,h,g,gamm0,
     +    c0,sigma,npo,ibv,inl,iob1,al,de,hmagv,elm,
     +    alpha,delta,hmagn,gamad,sig,axes,npo1)
      IMPLICIT NONE      
c ============= input ====================================================
c elements and epoch times, covariance and normal matrices at t0,
c sigmas for the boundary
      character*3 coo
      DOUBLE PRECISION east0(6),t0,t1,gamm0(6,6),c0(6,6),sigma
c number of points, flag for confidence bd/line of variation, nonlinearity
      INTEGER npo,ibv,inl,iob1
c magnitude
      DOUBLE PRECISION h,g,hmagn
c station code
      INTEGER idsta
c ============= output ===================================================
c points on the confidence boundary (difference w.r. to alpha,delta)
c WARNING! the output number of points is npo1.le.npo; 
c this beacuse hyperbolic points are discarded
      INCLUDE 'npoint.h'
      INCLUDE 'jplhdr.h'
      INTEGER npo1
      DOUBLE PRECISION al(npo),de(npo),hmagv(npo),allin,delin
c line of elements
      DOUBLE PRECISION elm(6,npo)
c best fit observations
      DOUBLE PRECISION alpha,delta
c covariance
      DOUBLE PRECISION gamad(2,2),axes(2,2),sig(2)
c ============END INTERFACE===============================================
c asteroid equinoctal elements, mean motion, covariance, normal matrices
      DOUBLE PRECISION eq(6),enne,gameq(6,6),ceq(6,6)  
c jacobian matrices of partial derivatives
      DOUBLE PRECISION dedeq(6,6),dedeqt(6,6),deqde(6,6),deqdet(6,6)
c return elements, workspace
      DOUBLE PRECISION east00(6),tmp(6,6)
c partial derivatives of alpha, delta, w.r. to elements (by columns)
      DOUBLE PRECISION daddet(6,2),dummy(6)
c second derivatives of alpha, delta, w.r. to elements (not used)
      DOUBLE PRECISION ddade(6,6),dddde(6,6)
c ===================================================================
c orthonormal basis, matrix defining the plane of the ellipse
      DOUBLE PRECISION v(6,6),ceicel(4,2)
c transformation matrix between the two planes
      DOUBLE PRECISION b(2,2)
c number of full revolutions around the sky
      INTEGER ng,nrev
c functions
      DOUBLE PRECISION appmag,prscag
c elongation,distance to Earth, distance to Sun (to compute magnitude)
      INCLUDE 'phase.h'
      DOUBLE PRECISION adot0,ddot0,pha0,dis0,dsun0,elov0,galla0
c ===================================================================
c constant of gravitation, trigonometric constants 
      INCLUDE 'sunmass.h'
      INCLUDE 'trig.h'
c temporaries, indexes
      DOUBLE PRECISION dal,ddl,maxsig,minsig
      INTEGER n
c flag for 2-body approximation; must be .false. for full n-body computation
      LOGICAL twobo
      twobo=.false.
c****************
c   static memory not required
c****************
c =====================================================================
      if(iob1/1000.eq.2.and.inl.eq.2)then
         WRITE(*,*)' preobn: mixing of radar and two-body '//
     +        'approximation not permitted'
         RETURN      
      endif
c coordinate change
      CALL cooder(east0,coo,gms,eq,'EQU',enne,deqde) 
c and inverse
      CALL cooder(eq,'EQU',gms,east00,coo,enne,dedeq)
c compute normal and covariance matrix for equinoctal elements
      CALL mulmat(deqde,6,6,gamm0,6,6,tmp)
      CALL transp(deqde,6,6,deqdet)
      CALL mulmat(tmp,6,6,deqdet,6,6,gameq)
      CALL transp(dedeq,6,6,dedeqt)
      CALL mulmat(dedeqt,6,6,c0,6,6,tmp)
      CALL mulmat(tmp,6,6,dedeq,6,6,ceq)
c =====================================================================
c compute observation; derivatives (of order 1) required            
c     write(*,*)'preobn calls alfdel', t0,t1
      if(iob1/1000.eq.1)then
         CALL alfdel (eq,t0,t1,idsta,alpha,delta,daddet(1,1),
     +        daddet(1,2),1,twobo,ddade,dddde)
c store true apparent motion, etc.
         adot0=adot
         ddot0=ddot
         pha0=pha
         dis0=dis
         dsun0=dsun
         elov0=elo
         galla0=gallat
c compute apparent magnitude at time of observation
         hmagn=appmag(h,g,dsun,dis,pha)
      elseif(iob1/1000.eq.2)then
         CALL rrdot (eq,iob1,t0,t1,idsta,alpha,delta,daddet(1,1),
     +        daddet(1,2),1,twobo)
         hmagn=0.d0
      ELSE
         WRITE(*,*)' preobn: this observation type not supported ',iob1
         RETURN           
      ENDIF

c =====================================================================
c compute ellipse of covariance of alpha,delta
      CALL ellips(daddet,gameq,sig,axes,gamad)
c If inl=0 then use automatic selection method
      IF(inl.eq.0)THEN
         maxsig=max(sig(1),sig(2))*degrad
         if(maxsig.le.1.0d0)then
            inl=1
c Is it safe to use two body if we may have close approach?
c         elseif(maxsig .le. 5.d0)then
c            inl=2
         else
            inl=3
         endif
      endif            
c If inl=0 then use automatic selection method
      IF(ibv.eq.0)THEN
         maxsig=max(sig(1),sig(2))
         minsig=min(sig(1),sig(2))
         if(maxsig/minsig.le.200.d0)then
            ibv=1
         else
            ibv=2
         endif
      endif    
      if(inl.eq.2)then
c 2-body aproximation for the central point
         CALL alfdel (eq,t0,t1,idsta,allin,delin,dummy,dummy,
     +        0,.true.,ddade,dddde)
      endif
c =====================================================================
c compute ellipse in the elements space 
      CALL slinel(daddet,gameq,ceq,ceicel,b,v)
c ===========================================================
c compute line of orbital elements
      CALL linobs(ibv,npo,eq,axes,sig,b,v,sigma,ceicel,elm,npo1)
c ===========================================================
      ng=0
      DO 7 n=1,npo1
c chose method to handle nonlinearity
        IF(inl.eq.1)THEN
c linear map from ellipse
           dal=prscag(6,elm(1,n),daddet(1,1))
           ddl=prscag(6,elm(1,n),daddet(1,2))
           al(n)=dal
           de(n)=ddl
           CALL vsumg(6,eq,elm(1,n),elm(1,n))
c apparent magnitude is the one of the nominal orbit
           hmagv(n)=hmagn
        ELSEIF(inl.eq.2)THEN
c 2-body propagation from ellipse
           CALL vsumg(6,eq,elm(1,n),elm(1,n))
           CALL alfdel (elm(1,n),t0,t1,idsta,al(n),de(n),
     +          dummy,dummy,0,.true.,ddade,dddde)
c compute apparent magnitude at time of observation
           hmagv(n)=appmag(h,g,dsun,dis,pha)
c difference is with respect to 2-body approx., used w.r. to true orbit
           al(n)=al(n)-allin
           de(n)=de(n)-delin
        ELSEIF(inl.eq.3)THEN
c full n-body propagation from ellipse 
           CALL vsumg(6,eq,elm(1,n),elm(1,n))
           CALL proele('EQU',t0,elm(1,n),t1,elm(1,n))
           if(iob1/1000.eq.1)then
              CALL alfdel (elm(1,n),t1,t1,idsta,al(n),de(n),
     +             dummy,dummy,0,twobo,ddade,dddde)
c other prediction data stored in common
              phav(n)=pha
              disv(n)=dis
              dsunv(n)=dsun
              elov(n)=elo
              gallav(n)=gallat
              adotv(n)=adot
              ddotv(n)=ddot
c compute apparent magnitude at time of observation
              hmagv(n)=appmag(h,g,dsun,dis,pha)
           elseif(iob1/1000.eq.2)then
              CALL rrdot (elm(1,n),iob1,t1,t1,idsta,al(n),de(n),
     +             dummy,dummy,0,twobo)
              hmagv(n)=0.d0
           ELSE
              stop'preobn: internal error'
           ENDIF
           al(n)=al(n)-alpha
           de(n)=de(n)-delta

        ELSE
           WRITE(*,*)' preobn: this we have not invented yet ', inl
           RETURN           
        ENDIF
c keep count of lost revolutions
        IF(n.eq.1)THEN
           IF(al(n).gt.pig)al(n)=al(n)-dpig
        ELSE
           CALL angupd(al(n),al(n-1),ng)
        ENDIF
c temporary output
        if(iob1/1000.eq.1)then
           write(*,*)n,', RA/DEC (deg)',al(n)*degrad,de(n)*degrad,ng
        elseif(iob1/1000.eq.2)then
           write(*,*)n,', R/RDOT (km,km/day)',al(n)*au,de(n)*au,ng
        endif
 7    continue
c =====================================================================
c ensure that LOV is consistent with nominal point
c first find midpoint of LOV, assume npo is even
      if(ibv.eq.2)then
         nrev=nint((al(npo/2)+al(npo/2+1))/2.d0/dpig)
         write(*,*)'debug: nrev:',nrev
         if(nrev.ne.0)then
            do n=1,npo1
               al(n)=al(n)-nrev*dpig
            enddo
         endif
      endif
c restore original apparent motion
      IF(iob1/1000.eq.1)THEN
         adot=adot0
         ddot=ddot0
         pha=pha0
         dis=dis0
         dsun=dsun0
         elo=elov0
         gallat=galla0
      ENDIF
      RETURN
      END
c ===========================================================
c common subroutines for preobn and fclan
c patch 1.6.1, A. Milani, May 2, 1998
c ===========================================================
c LINOBS defines line of changes in orbital elements to be used for
c confidence boundary/variations line 
c ===========================================================
      SUBROUTINE linobs(ibv,npo,eq,axes,sig,b,v,sigma,ceicel,elm,npo1)
      IMPLICIT NONE
c ====================INPUT==================================
      INTEGER ibv,npo
      DOUBLE PRECISION eq(6)
      DOUBLE PRECISION axes(2,2),sig(2),b(2,2),sigma
c matrix defining the plane of the ellipse,new orthonormal reference
      DOUBLE PRECISION ceicel(4,2),v(6,6)
c ===================OUTPUT==================================
      INTEGER npo1
      DOUBLE PRECISION elm(6,npo)
c ==================END INTERFACE============================
      INTEGER nn,n,i
      DOUBLE PRECISION s,x,y,vad(2),xv,yv,dn,dth,theta,xa,yd
      DOUBLE PRECISION eqnew(6)
      DOUBLE PRECISION alde(2),ecc
      INCLUDE 'trig.h'
c =====================================================================
c line of maximum variation: in the alpha-delta plane
      DO i=1,2
        vad(i)=axes(i,2)*sig(2)
      ENDDO
c in the elements space
      xv=(b(1,1)*vad(1)+b(1,2)*vad(2))
      yv=(b(2,1)*vad(1)+b(2,2)*vad(2))
c direction not used any more 
c     theta0=atan2(yv,xv)
c linear step for variation axis parametrisation
      dn=2.d0/float(npo-1)
c angular step for ellipse parametrisation
      dth=dpig/float(npo)
c ===========================================================
c main loop on the number of output points
      nn=0
      DO 7 n=1,npo
c ===========================================================
c choice between two output options
        IF(ibv.eq.2)THEN
c ===========================================================
c line of maximum variation in the elements space
           s=(n-1)*dn-1.d0
           x=sigma*s*xv
           y=sigma*s*yv
        ELSEIF(ibv.eq.1)THEN
c =====================================================================
c parametrisation of the ellipse in the subspace of elements, based upon the
c parametrisation of the ellipse in the alpha-delta plane
c WARNING: npo must be divisible by 2, otherwise one tip of the
c banana would be missed
           theta=(n-1)*dth
           xa=sig(1)*cos(theta)*sigma
           yd=sig(2)*sin(theta)*sigma
           CALL lincog(2,axes(1,1),xa,axes(1,2),yd,alde)
c transfer of parametrisation in the V1,V2 plane
           x=(b(1,1)*alde(1)+b(1,2)*alde(2))
           y=(b(2,1)*alde(1)+b(2,2)*alde(2))
        ELSE
           write(*,*)' linobs: this should not happen,ibv=',ibv
        ENDIF
c compute displacement on the confidence ellipsoid corresponding to x,y
        nn=nn+1    
        CALL elemov(x,y,v,ceicel,elm(1,nn))
c add to the original center of the ellipsoid of confidence
        CALL vsumg(6,eq,elm(1,nn),eqnew)
        ecc=sqrt(eqnew(2)**2+eqnew(3)**2)
        IF(ecc.ge.1.d0.or.eqnew(1).le.0.d0)THEN
           write(*,*)' Hyperbolic, ecc=',ecc,' a=',eqnew(1)
           nn=nn-1
        ELSEIF(ecc.ge.0.99d0)THEN
           write(*,*)' Almost Hyperbolic, ecc=',ecc,' a=',eqnew(1)
           nn=nn-1
        ENDIF
 7    continue
c final count of non hyperbolic orbits
      npo1=nn  
      RETURN
      END
c =====================================================================
c ELLIPS
c compute covariance ellipse of two observables
c =====================================================================
      SUBROUTINE ellips(daddet,gamm0,sig,axes,gamad)
      IMPLICIT NONE
c input covariance matrix
      DOUBLE PRECISION gamm0(6,6)
c input partial derivatives of alpha, delta, w.r. to elements (by columns)
      DOUBLE PRECISION daddet(6,2)
c output covariance
      DOUBLE PRECISION gamad(2,2),axes(2,2),sig(2)
c ==============END INTERFACE==========================================
c eigenvalues, workspace, transposed
      DOUBLE PRECISION eigval(2),tmp26(2,6),fv1(2),fv2(2),dadde(2,6)
c loop indexes
      INTEGER i
c error flag
      INTEGER ierr
c =====================================================================
      CALL transp(daddet,6,2,dadde)
      CALL mulmat(dadde,2,6,gamm0,6,6,tmp26)
      CALL mulmat(tmp26,2,6,daddet,6,2,gamad)
c =====================================================================
c compute ellipse of confidence
c eigenvalues
      CALL rs(2,2,gamad,eigval,1,axes,fv1,fv2,ierr)
      DO  i=1,2
        IF(eigval(i).gt.0.d0)THEN
           sig(i)=sqrt(eigval(i))
        ELSE
           write(*,*) 'non positive eigenvalue'
           sig(i)=0.d0
        ENDIF
      ENDDO
      RETURN
      END
c =====================================================================
c ELEMOV
c compute displacement on the confidence ellipsoid corresponding to x,y
c on the plane of the gradients of alpha-delta
c =====================================================================
      SUBROUTINE elemov(x,y,v,ceicel,del)
      IMPLICIT NONE
c inout/output
      DOUBLE PRECISION x,y,v(6,6),del(6)
      DOUBLE PRECISION ceicel(4,2)
c workspace
      DOUBLE PRECISION dee(4),deel(6)
c ===================
      CALL lincog(6,v(1,1),x,v(1,2),y,del)
      CALL lincog(4,ceicel(1,1),-x,ceicel(1,2),-y,dee)
      CALL mulmav(v(1,3),6,4,dee,4,deel)
      CALL vsumg(6,del,deel,del)
      RETURN
      END
c ===========================================================
c SLINEL
c semilinear boundary ellipse computation
c =========================================================== 
      SUBROUTINE slinel(dtpdet,gc,cc,ceicel,b,v)
      IMPLICIT NONE
c 6 by 2 matrix with columns= gradients
      DOUBLE PRECISION dtpdet(6,2)
c normal and covariance matrices
      DOUBLE PRECISION gc(6,6),cc(6,6)
c orthonormal basis
      DOUBLE PRECISION v(6,6),vt(6,6),gamv(6,6),cv(6,6),tmp(6,6)
c partial matrices
      DOUBLE PRECISION c4(4,4),cinv(4,4),c42(4,2),ceicel(4,2)
c line of maximum variation
      DOUBLE PRECISION a(2,2),b(2,2),deta
c loop indexes ii=1,2, ij,ijj=1,4
      INTEGER ii, ij, ijj
c for inversion with tcholevski: workspace, error flag
      DOUBLE PRECISION ws(4)
      INTEGER ierr
      DOUBLE PRECISION prscag
c =====================================================================
c adapted orthonormal basis, covariance and normal matrix in the new basis
      CALL graha(dtpdet,6,v)
      CALL transp(v,6,6,vt)
      CALL mulmat(vt,6,6,gc,6,6,tmp)
      CALL mulmat(tmp,6,6,v,6,6,gamv)
      CALL mulmat(vt,6,6,cc,6,6,tmp)
      CALL mulmat(tmp,6,6,v,6,6,cv)
c =====================================================================
c 4x4 and 4x2 submatrices of normal matrix
      do 15 ijj=1,4
        DO ij=1,4
          c4(ijj,ij)=cv(ijj+2,ij+2)
        ENDDO
        DO  ii=1,2
          c42(ijj,ii)=cv(ijj+2,ii)
        ENDDO
 15   continue
c ===========================================================
c Cholewski method for inversion
      CALL tchinv(c4,4,cinv,ws,ierr)
      IF(ierr.ne.0)THEN
         write(*,*)' decide what to do, ierr=',ierr
      ENDIF
c ===========================================================
c matrix to be used for out of plane component
      CALL mulmat(cinv,4,4,c42,4,2,ceicel)
c ===========================================================
c linear map from the elements space (with base V) and the alpha-delta plane 
      a(1,1)=prscag(6,dtpdet(1,1),v(1,1))
      a(1,2)=prscag(6,dtpdet(1,1),v(1,2))
      a(2,1)=prscag(6,dtpdet(1,2),v(1,1))
      a(2,2)=prscag(6,dtpdet(1,2),v(1,2))
      CALL inv22(a,b,deta)
      RETURN
      END
c **********************************************************
c  ANGUPD
c   given the principal value ang of an angle, and an angle vang 
c   with ng revolutions, the routine 
c   finds a new new ng in the assumption that
c   less than half a revolution has occurred; ang is modified to 
c   include the multiples of dpig
c   in case ang is not between -pig and pig, e.g. because it is 
c   between 0 and dpig, it is first reduced to principal value
      SUBROUTINE angupd(ang,vang,ng)
      IMPLICIT NONE
      DOUBLE PRECISION ang,vang,d
      INTEGER ng,ig
      INCLUDE 'trig.h'
c first reduce to principal value
      if(ang.gt.pig)ang=ang-dpig
      if(ang.lt.-pig)ang=ang+dpig
c count revolutions
      ig=nint((vang-ang)/dpig)
      ang=ang+dpig*ig
      d=ang-vang
c correct by one revolution if necessary
      if(d.gt.pig)then
         ng=ng-1
         ang=ang-dpig
      elseif(d.lt.-pig)then
         ng=ng+1
         ang=ang+dpig
      endif
      return
      end
c ====================================================================
c Graham- Schmidt procedure to generate an orthonormal basis v
c starting from 2  n-vectors a
c The new basis must be such that the first 2 vectors are a basis
c for the space spanned by the 2 columns of a
      SUBROUTINE graha(a,n,v)
      implicit none
      integer n,nx
      parameter (nx=10)
      double precision a(n,2),v(n,n)
      integer j,jok,jj
      double precision prscag,cc,cc1,cc2,epsi,vl
      double precision ws(nx)
      logical ize
c dimension check
      if(n.gt.nx)then
         write(*,*)'n =',n,' larger than nx=',nx,' in graha'
         stop
      endif 
c selection of the control for "zero" vectors
      cc1=sqrt(prscag(n,a(1,1),a(1,1)))
      cc2=sqrt(prscag(n,a(1,2),a(1,2)))
      epsi=1.d-12*min(cc1,cc2)
      if(epsi.eq.0.d0)then
         write(*,*)' a has rank zero'
c        stop
      endif
c start by orthonormalisation of the space spanned by the columns of a
c
c V1 is the versor of A1
      call versor(n,a(1,1),epsi,v(1,1),vl,ize)
      if(ize)then
         write(*,*)' first vector of a is too small'
c        stop
      endif 
c the following vectors are obtained
c by removing the components along the previous ones
      cc=-prscag(n,v(1,1),a(1,2))
      call lincog(n,a(1,2),1.d0,v(1,1),cc,v(1,2))
      call versor(n,v(1,2),epsi,v(1,2),vl,ize)
      if(ize)then
         write(*,*)' a has practically rank one'
c        stop
      endif
c we now use the vectors of the canonic basis to supplement the span of A1,A2
      jok=0
      do 1 j=1,n
c remove the components along span(A), that is along V1 and V2
        cc1=-v(j,1)
        cc2=-v(j,2)
        call lincog(n,v(1,1),cc1,v(1,2),cc2,ws)
        ws(j)=ws(j)+1.d0
        call versor(n,ws,epsi,v(1,3+jok),vl,ize)
        if(.not.ize)then
c now V(3+jok) is orthogonal to span(A); remove the components along
c the previous ones (unless it is the first)
           if(jok.gt.0)then
              do  jj=1,jok
                cc=-prscag(n,v(1,3+jok),v(1,2+jj))
                call lincog(n,v(1,3+jok),1.d0,v(1,2+jj),cc,v(1,3+jok))
              enddo
              call versor(n,v(1,3+jok),epsi,v(1,3+jok),vl,ize)
              if(ize)then
                 goto 1
              endif
           endif
c the new versor is a good one
           jok=jok+1
           if(jok.eq.n-2)then
              goto 2
           endif
        endif
 1    continue
 2    continue
      if(jok.lt.n-2)then
         write(*,*)' something went wrong, jok=',jok
      endif
      return
      end
      SUBROUTINE versor(n,a,epsi,b,vl,ize)
      implicit none
      integer n,i
      logical ize
      double precision a(n),b(n),epsi,prscag,vl
      vl=sqrt(prscag(n,a,a))
      if(vl.lt.epsi)then
         ize=.true.
      else
         ize=.false.
         do  i=1,n
           b(i)=a(i)/vl
         enddo
      endif
      return
      end
c
c =====================================================================
c OUTOBC
c =====================================================================
c  output of predicted observation, possibly with confidence ellipse
c   input: iun   = output unit
c          iobs  = observation type
c          ids = station code
c          t1 = time of observation (UTC)
c          alpha, delta, hmagn = observation
c          adot,ddot = proper motion
c          elo,dis = elongation, distance from Earth
c
c          icov  = 1 for observations only, 2 to add confidence ellipse
c          gamad,sig,axes = covariance matrix, sigmas along axes 
c                      (only for icov=2, otherwise dummy)
c =====================================================================
      SUBROUTINE outobc(iun,iobs,ids,t1,alpha,delta,hmagn,adot,ddot,
     +     elo,dis,icov,gamad,sig,axes)
      implicit none
      include 'trig.h'
c needs AU value in km
      INCLUDE 'jplhdr.h'
c output unit, station code, obs. type
      integer iun,ids,iobs
c observations
      double precision t1,alpha,delta,hmagn,adot,ddot,elo,dis
c covariance
      integer icov
      double precision gamad(2,2),axes(2,2),sig(2)
c ================end interface===============================
      double precision princ  
      integer i,j
c time variables
      integer ideg,iday,imonth,iyear,ihour,imin,isec,ln,truncat
      double precision hour,minu,sec
      CHARACTER*22 timstr
      CHARACTER*19 tmpstr
      CHARACTER*12 rastri,rdstri
      CHARACTER*1 signo
c convert time
      CALL mjddat(t1,iday,imonth,iyear,hour)
c convert hour to 12:12:12
      ihour=truncat(hour,1d-7)
      minu=(hour-ihour)*60.d0
      imin=truncat(minu,1d-5)
      sec=(minu-imin)*60.d0
      isec=truncat(sec,1d-3)
      WRITE(timstr,192) iyear,imonth,iday,ihour,imin,isec,sec-isec
 192  FORMAT(I4,'/',I2.2,'/',I2.2,1x,I2.2,':',I2.2,':',I2.2,f3.2)
c =================== select by observation type ===================
      IF(iobs/1000.eq.1)THEN
c %%%%%%%%%%%% ASTROMETRY %%%%%%%%%%%%%%%%
c convert RA
         alpha=princ(alpha)
         CALL sessag(alpha*degrad/15.d0,signo,ihour,imin,sec)
         IF(signo.eq.'-')STOP 'wrirms error: negative right ascension.'
c prepare RA string
         WRITE(tmpstr,FMT='(F6.3)') sec
         CALL rmsp(tmpstr,ln)
         IF(ln.lt.6)tmpstr='0'//tmpstr
         WRITE(rastri,130)ihour,imin,tmpstr
 130     FORMAT(I2.2,':',I2.2,':',a6)
c convert DEC
         CALL sessag(delta*degrad,signo,ideg,imin,sec)
c prepare DEC string
         WRITE(tmpstr,FMT='(F5.2)') sec
         CALL rmsp(tmpstr,ln)
         IF(ln.lt.5)tmpstr='0'//tmpstr
         WRITE(rdstri,170)signo,ideg,imin,tmpstr
 170     FORMAT(A1,I2.2,1x,I2.2,1x,a5)

         write(iun,101)timstr,t1,ids,
     +        rastri,alpha*degrad,
     +        rdstri,delta*degrad,
     +        secrad*adot/24.d0,secrad*ddot/24.d0,
     +        dis,elo*degrad,hmagn
         write(*,101)timstr,t1,ids,
     +        rastri,alpha*degrad,
     +        rdstri,delta*degrad,
     +        secrad*adot/24.d0,secrad*ddot/24.d0,
     +        dis,elo*degrad,hmagn
 101     format('Astrometric Observation Prediction'/
     +        'For ',a19,' (UTC); ',f12.5,'(MJD)'/
     +        'Observatory code= ',i3.3/
     +        'RA= ',a12,' (HH:MM:SS); ',f11.5,' (deg)'/
     +        'DEC= ',a12,' (deg min sec); ',f11.5,' (deg)'/
     +        'RA/DEC Apparent motion=',2(2x,f9.2),' (arcsec/hour)'/
     +        'Earth distance= ',f6.4,' (AU)'/
     +        'Solar elongation= ',f6.2,' (deg)'/
     +        'Apparent magnitude= ',f5.2)
         IF(icov.eq.1)RETURN
c rescaling in arcsec
         do  i=1,2
            sig(i)=sig(i)*secrad
            do  j=1,2
               gamad(i,j)=gamad(i,j)*secrad**2
            enddo
         enddo
         write(iun,201)(sig(j),(axes(i,j),i=1,2),j=1,2)
         write(*,201)(sig(j),(axes(i,j),i=1,2),j=1,2)
 201     format(
     +'Size and orientation of 1-sigma uncertainty ellipse'/
     +'Short axis : Size= ',1p,g12.6 ,' (arcsec); Direction= ',
     + 0p,2(1x,f8.5)/
     +'Long axis : Size= ',1p,g12.6 ,' (arcsec); Direction= ',
     + 0p,2(1x,f8.5))
      ELSEIF(iobs/1000.eq.2)THEN
c %%%%%%%%%%%% RADAR %%%%%%%%%%%%%%%%
         write(iun,102)t1,ids,alpha*au,delta*au
         write(*,102)t1,ids,alpha*au,delta*au
 102     format('time, MJD=',f13.6,'  station=',i4/
     +       ' range (KM)         = ',f16.5/
     +       ' range rate (KM/DAY)=  ',f15.5)
         IF(icov.eq.1)RETURN
c rescaling in km, km/day
         do  i=1,2
            sig(i)=sig(i)*au
            do  j=1,2
               gamad(i,j)=gamad(i,j)*au**2
            enddo
         enddo
         write(iun,202)(sig(j),(axes(i,j),i=1,2),j=1,2)
         write(*,202)(sig(j),(axes(i,j),i=1,2),j=1,2)
 202     format(' in the range (KM), range-rate (KM/DAY) plane'/
     +          ' sigma1 = ',1p,g14.7 ,' axis1= ',2(1x,g12.5)/
     +          ' sigma2 = ',1p,g14.7 ,' axis2= ',2(1x,g12.5))
      ELSEIF(iobs/1000.eq.4)THEN
c %%%%%%%%%%%% PROPER MOTION %%%%%%%%%%%%%%%%
         write(iun,109)t1,ids,alpha*secrad/24.d0,delta*secrad/24.d0
         write(*,109)t1,ids,alpha*secrad/24.d0,delta*secrad/24.d0
 109     format('time, MJD=',f13.6,'  station=',i4/
     +       ' RA motion (arcsec/hour)     = ',f9.2/
     +       ' DEC motion (arcsec/hour)    = ',f9.2)
         IF(icov.eq.1)RETURN
c rescaling in arcsec/hour
         do  i=1,2
            sig(i)=sig(i)*secrad/24.d0
            do  j=1,2
               gamad(i,j)=gamad(i,j)*(secrad/24.d0)**2
            enddo
         enddo
         write(iun,209)(sig(j),(axes(i,j),i=1,2),j=1,2)
         write(*,209)(sig(j),(axes(i,j),i=1,2),j=1,2)
 209     format('sigma1 (arcsec/hr)= ',1p,g14.7 ,' axis1= ',2(1x,g12.5)/
     +          'sigma2 (arcsec/hr)= ',1p,g14.7 ,' axis2= ',2(1x,g12.5))
      ELSE
         WRITE(*,*)'outobs: iobs=',iobs,' not understood'
      ENDIF
      return
      end

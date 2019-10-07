c ====================================================
c FOBPRE predict observations
c ====================================================
      SUBROUTINE fobpre(icov,ini0,cov0,iun20,iun8,ok,
     +   titnam,filnam,t0,eq0,h0,gmag,g0,c0,ids,iob1,t1,tut,aobs,dobs,
     +   t2,dt,astnam)
      IMPLICIT NONE
c =================INPUT=========================================
c requirements on covariance, handling of nonlinearity
      INTEGER icov,inl
c availability of initial conditions, covariance, all required
      LOGICAL ini0,cov0,ok
c output units
      INTEGER iun20,iun8
c epoch time, observation time, end of ephemerides time, step
      DOUBLE PRECISION t0,t1,t2,dt
c elements, absolute magnitude, covariance and normal matrix
      DOUBLE PRECISION eq0(6),h0,g0(6,6),c0(6,6)
c asteroid name etc.
      CHARACTER*80 titnam
      CHARACTER*60 filnam
c actual observations for comparison
      DOUBLE PRECISION aobs,dobs
c asteroid name
      CHARACTER*(*) astnam
c ===== predicted observations ===========
c UTC of observation
      DOUBLE PRECISION tut
c angles (best fit prediction)
      DOUBLE PRECISION alpha,delta
c station code, observation type
      INTEGER ids,iob1
c covariance of the observations
      DOUBLE PRECISION gamad(2,2),axes(2,2),sig(2)
c opposition effect, apparent magnitude of nominal solution
      DOUBLE PRECISION gmag,hmagn
c confidence boundary, line of max variation (alpha, delta, app. magnitude)
      INTEGER npo, npox, ibv, npo1, npop
      PARAMETER (npox=4000)
      DOUBLE PRECISION sigma, al(npox),de(npox),hmagv(npox)
c line of elements
      DOUBLE PRECISION elm(6,npox)
c multiple data for confidence boundary
      INCLUDE 'npoint.h'
c ====================================== 
c observation data (proper motion, elongation, distance)
      INCLUDE 'phase.h'
c menu empty items
      CHARACTER*20 menunam
      CHARACTER*70 s4,s5,s6,s7,s8,s9,s10
      INTEGER i
c epehemrides output
      CHARACTER*100 file,fields
      CHARACTER*3 scale
      DOUBLE PRECISION mass
      INTEGER ln,iuneph
c =====================================================================
c options
c ===================================================================
c ===================================================================
c chose handling of nonlinearity
 57   IF(icov.eq.3.or.icov.eq.4)THEN
         menunam='prednonl'
         CALL menu(inl,menunam,3,'How to handle nonlinearity?=',
     +         'linear map=',
     +         '2-body nonlinearity=',
     +         'full n-body nonlinearity=',
     +          s4,s5,s6,s7,s8,s9,s10)
         IF(inl.eq.0)GOTO 57
      ENDIF 
c =====================================================================
c check availability of initial conditions (also covariance for icov=2)
      CALL chereq(icov,ini0,cov0,t0,iun20,iun8,ok)
      IF(.not.ok)RETURN
c =====================================================================
c check availability of JPL ephemerides and ET-UT table
      CALL chetim(t1,t1,ok)
      IF(.not.ok)RETURN
c =====================================================================
c compute prediction; without and with covariance
c =====================================================================
      IF(icov.eq.1)THEN
c =====================================================================
c only alpha, delta, magnitude
         CALL preobs('EQU',t0,ids,t1,eq0,iob1,alpha,delta,
     +        h0,gmag,hmagn)
         CALL outobc(iun20,iob1,ids,tut,alpha,delta,hmagn,adot,ddot,
     +     elo,dis,icov,gamad,sig,axes)
      ELSEIF(icov.eq.2)THEN
c =====================================================================
c alpha, delta, magnitude, covariance and ellipse of confidence
         CALL preobc('EQU',t0,ids,t1,eq0,h0,gmag,g0,
     +        iob1,alpha,delta,hmagn,gamad,sig,axes)
         CALL outobc(iun20,iob1,ids,tut,alpha,delta,hmagn,adot,ddot,
     +     elo,dis,icov,gamad,sig,axes)
c ===================================================================
c generation of sky epehemrides
      ELSEIF(icov.eq.5)THEN
c check availability of JPL ephemerides and ET-UT table for entire time span
         CALL chetim(t1,t2,ok)
         IF(.not.ok)RETURN
         IF(nint(abs(t2-t1)/dt).gt.500)THEN
            write(*,*)' Too many ephemerides points:',
     +           nint(abs(t2-t1)/dt)
            write(*,*)'Select a time interval and span to ',
     +           'ensure that there are fewer than 500 points.'
         ELSE
c open ephemerides file in current directory
            file=astnam//'.eph'
            CALL rmsp(file,ln)
            CALL filopn(iuneph,file(1:ln),'unknown')
            fields='cal,mjd,coord,mag,elong,glat,r,delta,appmot,skyerr'
            scale='UTC'
            CALL ephemc(iuneph,'EQU',t0,eq0,g0,.true.,t1,t2,
     +           dt,mass,h0,gmag,ids,scale,fields)
            CALL filclo(iuneph,' ')
            WRITE(*,*)' Generated ephemeris in file: ',file(1:ln)
         ENDIF
c =====================================================================      
         ELSEIF(icov.eq.3.or.icov.eq.4)THEN
c =====================================================================
c alpha, delta, magnitude, covariance and confidence boundary;
c input specification of set of points
         CALL asscbd(iun20,npox,npo,sigma,ibv)          
c =====================================================================
c compute prediction, boundary
         CALL preobn('EQU',t0,ids,t1,eq0,h0,gmag,g0,
     +             c0,sigma,npo,ibv,inl,iob1,al,de,hmagv,elm,
     +             alpha,delta,hmagn,gamad,sig,axes,npo1)
         CALL outobc(iun20,iob1,ids,tut,alpha,delta,hmagn,adot,ddot,
     +     elo,dis,icov,gamad,sig,axes)
         IF(npo1.le.0)THEN
            WRITE(*,*)'fobpre: no elliptic orbits ',npo1
            RETURN
         ENDIF
         IF(ibv.eq.1)THEN
c confidence boundary; one point added to close line
            al(npo1+1)=al(1)
            de(npo1+1)=de(1)
            hmagv(npo1+1)=hmagv(1)
            DO i=1,6
              elm(i,npo1+1)=elm(i,1)
            ENDDO
            IF(inl.eq.3)disv(npo+1)=disv(1)
            npop=npo1+1
         ELSEIF(ibv.eq.2)THEN
c line of variations
            npop=npo1
         ENDIF
c if no observation is given, use the nominal marked with a cross
         IF(icov.eq.3)THEN
            aobs=alpha
            dobs=delta
         ENDIF     
c =====================================================================
c output observation, apparent motion, confidence boundary
         CALL outmul(titnam,filnam,tut,sigma,alpha,delta,
     +              al,de,hmagv,1,npop,1,icov-2,aobs,dobs,iob1)
      ENDIF
      RETURN
      END









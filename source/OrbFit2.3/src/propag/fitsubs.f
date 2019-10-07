c ============================================
c  FITSUBS
c ============================================
c routines used by orbit determination interactive programs, such as fitobs
c
c  CONTAINS:
c
c   whichcor   selection of parameters to be solved for
c   chereq     check availability of data required
c   cheobs     id for observations
c   chetim     id for time series such as JPL ephem.)
c   seleph     interactive selection of interval
c   asstim        "          "       of time
c   asscbd        "          "       of confidence boundary
c   orbsel        "          "       of orbit
c   start      initial guess for identification (by interpolation)
c ===================================================================
c WHICOR
c ===================================================================
c interrogation routine for inew, icor
      SUBROUTINE whicor(inter,icor,ncor,inew)
      implicit none
      integer iansw,inew,ncor,icor(6),inter,i
c Choose method (NOW FIXED AT PSEUDO-NEWTON)
      inew=2
c60   write(*,*)' 1=true Newton 2=pseudo Newton'
c     read(*,*)inew
c     if(inew.ne.1.and.inew.ne.2)then
c           write(*,*)'This we have not invented yet'
c           goto 60
c     endif
c interactive/automatic default version
      if(inter.eq.0)then
         do 59 i=1,6
           icor(i)=1
 59      enddo
         ncor=6
         return
      else
c interactive version
c
c Component of orbital element vector that need to be corrected
         ncor=0
         do 63 i=1,6
           write(*,*)'Element:  ',i,   ' 1=correct, 0=no'
           read(*,*) iansw
           if(iansw.ne.0)then
              ncor=ncor+1
              icor(i)=iansw
           else
              icor(i)=0
           endif
 63      enddo
      endif
      return
      end
c ====================================================================
c CHEREQ
c ====================================================================
c checking the availability of the required data, for propagation
c and preciction of observations
c ====================================================================
      SUBROUTINE chereq(icov,ini,cov,t,iun20,iun8,ok)
      implicit none
      integer icov
      logical ini,cov
      integer iun20,iun8
      double precision t
      logical ok
c ===================================================================
c logical check: initial conditions available
c ===================================================================
      if(ini)then
         IF(iun20.gt.0)write(iun20,222) t
         ok=.true.
      else
         write(*,*)' initial conditions not available'
         ok=.false.
      endif
c ===================================================================
c logical check: if covariance has to be propagated,
c it needs to be available for the initial time
c ===================================================================
      if(icov.gt.1)then
         if(.not.cov)then
            write(*,*)' initial covariance not available'
            ok=.false.
         elseif(iun8.gt.0)THEN
            write(iun8,*)
            write(iun8,222)t
         endif 
      endif
  222 format('  Initial epoch (MJD): ',f8.1)
      return
      end
c ====================================================================
c CHEOBS
c ====================================================================
c check availability of observations and initial condition
c ===================================================================
      SUBROUTINE cheobs(obs,ini,ok)
      logical obs,ini,ok
      if(ini)then
         if(obs)then
            ok=.true. 
         else
            write(*,*)'missing observations for this arc'
            ok=.false.
         endif
      else
         write(*,*)'missing initial conditions for this arc'
         ok=.false.
      endif
      return
      end
c ====================================================
c CHETIM check availability of required time-dependent data
c     input: t1,t2 dates MJD; assumed t1.le.t2
c     ok: .true. if available
c ====================================================
      SUBROUTINE chetim(t1,t2,ok)
      implicit none
      double precision t1,t2
      logical ok
c ======== time spans for JPL data etc. =========
      include 'timespan.h'
      ok=.true.
c ==================================================================
c check availability of JPL ephemerides
      if(t1.lt.tejpl1.or.t2.gt.tejpl2)then
         write(*,*)' JPL epehemerides not available for =',t1,t2
         write(*,*)' but only for interval ',tejpl1,tejpl2
         ok=.false.
      endif
c ===================================================================
c check availability of ET-UT table
      if(t1.lt.temut1.or.t2.gt.temut2)then
         if(temute) then
c           write(*,*)' ET-UT not available for ',t1,t2
c           write(*,*)' but only for interval ',temut1,temut2
c           write(*,*)' however, extrapolation will be used'
         else
            write(*,*)' ET-UT not available for ',t1,t2
            write(*,*)' but only for interval ',temut1,temut2
            ok=.false.
         endif
      endif
      return
      end
c ===================================================
c SELEPH
c select time interval, step
      SUBROUTINE seleph(tut1,tdt1,tut2,tdt2,dt,idsta)
      IMPLICIT NONE
c output
      DOUBLE PRECISION  tut1,tdt1,tut2,tdt2,dt     
      INTEGER idsta
c times in various formats used internally
      CHARACTER*3 scale
      INTEGER mjd,mjdtdt
      DOUBLE PRECISION sec,sectdt
      WRITE(*,*)' Initial time (MJD UTC)?'
      READ(*,*)tut1
      WRITE(*,*)' Final time (MJD UTC)?'
      READ(*,*)tut2
      WRITE(*,*)' Time interval (days)?'
      READ(*,*)dt
      WRITE(*,*)' Observatory code?'
      READ(*,*)idsta
      scale='UTC'
c =========== TIME CONVERSION ================
c starting time
      mjd=tut1
      sec=(tut1-mjd)*86400
      call cnvtim(mjd,sec,scale,mjdtdt,sectdt,'TDT')
      tdt1=sectdt/86400.d0+float(mjdtdt)
c stopping time
      mjd=tut2
      sec=(tut2-mjd)*86400
      call cnvtim(mjd,sec,scale,mjdtdt,sectdt,'TDT')
      tdt2=sectdt/86400.d0+float(mjdtdt)
      RETURN
      END
c ====================================================
c ASSTIM assign time for predictions
c operations controlled by flag icov
c  ICOV=1,2,3  ask the user to assign prediciton time t1
c  ICOV=4      ask the user to assign observation number im
c                      then t1=tau(im)
c  for ICOV.ne.4 also station code ids and obs. type iob1
c        have to be assigned by the user
c ====================================================
      SUBROUTINE asstim(icov,iobs,tau,tut,idsta,m,mall,im,
     +       iob1,t1,tut1,ids)
      IMPLICIT NONE
c ==============input=============
      INTEGER icov,m,mall,mp,idsta(mall),iobs(mall)
      DOUBLE PRECISION tau(mall),tut(mall)
c ==============output=============
      INTEGER im,iob1,ids
      DOUBLE PRECISION t1,tut1
c time conversion
      INTEGER mjd1,mjd2,intlo
      DOUBLE PRECISION sec1,sec2
c ===================================================================
      if(mall.lt.m)then
         write(*,*)'asstim: this should not happen, m,mall ',m,mall
         mp=0
      else
         mp=mall-m
      endif
c assign observation time
      if(icov.eq.4)then
c compare confidence boundary with observations
 184     write(*,*)' observed arcs: from, to, no. obs'
         write(*,182)tau(1),tau(m),m
 182     format('arc 1: ',2f8.1,i6)
         IF(mall.gt.m)THEN
            write(*,183)tau(m+1),tau(mall),mp
 183        format('arc 2: ',2f8.1,i6)
         ENDIF
         write(*,*)' observation number?   '
         read(*,*)im
         if(im.lt.0.or.im.gt.mall)then
            write(*,*)' observation no. im=',im,' not available'
            goto 184
         endif
         t1=tau(im)
         tut1=tut(im)
         ids=idsta(im)
         iob1=iobs(im)
      else
c dummy obs. number (not used, to avoid out of bounds)
         im=1
c assign arbitrary time
         write(*,*)' give time of prediction (MJD)   '
         read(*,*)t1
c universal time of the required observation 
         mjd1=intlo(t1)
         sec1=(t1-float(mjd1))*86400.d0
         CALL cnvtim(mjd1,sec1,'TDT',mjd2,sec2,'UTC')
         tut1=sec2/86400.d0+float(mjd2)
c         write(*,*)t1,tut1
c assign observation type
c185     write(*,*)' observation type 1=RA,DEC  9=proper motion?'
 185     write(*,*)' observation type 1=RA,DEC 2=radar 4=proper motion?'
         read(*,*)iob1
         IF(iob1.ne.1.and.iob1.ne.2.and.iob1.ne.4)GOTO 185
         iob1=iob1*1000            
c assign observatory code
         write(*,*)
         write(*,*)' observatory code (geocenter=500)?   '
         read(*,*) ids
c secret nationalistic feature
         if(ids.lt.0)then
            ids=599
         endif
      endif
      return
      end

c ===================================================================
c  ASSCBD
c ===================================================================
c alpha, delta, magnitude, covariance and confidence boundary;
c input specification of set of points
c ===================================================================
      SUBROUTINE asscbd(iun20,npox,npo,sigma,ibv)
      IMPLICIT NONE
      INTEGER npox,npo,ibv,iun20
      DOUBLE PRECISION sigma
      WRITE(*,*)' How many sigmas?   '
      READ(*,*) sigma
 1    WRITE(*,*)' 1=confidence boundary 2=line of max variation 0=auto'
      READ(*,*) ibv
      IF(ibv.ne.1.and.ibv.ne.2.and.ibv.ne.0)THEN
         WRITE(*,*)' option not understood ',ibv
         GOTO 1
      ENDIF
      WRITE(*,*)' how many points (even, please)?   '
      READ(*,*) npo
      IF(npo+2.gt.npox)THEN
         WRITE(*,*)' npo=',npo,' too large for npox=',npox
         npo=npox-2
         WRITE(*,*)' used npox-2'
      ENDIF
      WRITE(iun20,*)'no. points ',npo     
      RETURN
      END

c =====================================================
c ORB_SEL
      SUBROUTINE orb_sel(ini0,inip,initwo,isel)
      IMPLICIT NONE
c input: initial conditions available?
      LOGICAL ini0,inip,initwo
c output: selection 1=arc1 2=arc2 3=joint solution
      INTEGER isel
c characters for menu 
      CHARACTER*20 menunam
      CHARACTER*70 s3,s4,s5,s6,s7,s8,s9,s10
c cases without questions
      IF(ini0.and.(.not.inip.and..not.initwo))THEN
         WRITE(*,*)' ARC 1'
         isel=1
         RETURN
      ELSEIF(inip.and.(.not.ini0.and..not.initwo))THEN
         WRITE(*,*)' ARC 2'
         isel=2
         RETURN
      ELSEIF(.not.ini0.and..not.inip.and..not.initwo)THEN
         WRITE(*,*)' No initial conditions available'
         isel=0
         RETURN
      ENDIF
c cases in which there is a choice
      menunam='orbsel'
      CALL menu(isel,menunam,3,'which orbit?=',
     +      'arc 1=','arc 2=',
     +      'joint computed orbit=',
     +      s4,s5,s6,s7,s8,s9,s10)
      RETURN
      END
c ================================================
c START   Estimate of the mean motion and semimajor axis
c         allowing to combine two arcs
c Input   eq0,eqp elements of both arcs
c         t0,tp epoch times of elements of 1st and 2nd arc
c         gms Gauss constant
c         iorb select orbit 1,2
c Output  ng number of revolutions 
c         eng mean motion value for central time
c         am  semimajor axis value
c         plm mean longitude value
c ===============INTERFACE========================
      SUBROUTINE start(eq0,eqp,t0,tp,gms,iorb,ng,enm,am,plm)
      implicit none
c ===========INPUT==================
      double precision eq0(6),eqp(6),t0,tp,gms
      integer iorb
c ===========OUTPUT=================
      double precision enm,am,plm
      integer ng
c ==========END INTERFACE===========
      double precision en1,pl12,pl12p,en2,pl21,pl21p,enm0,enmp
      double precision tm,pll1,pll2
      integer ng0,ngp
c functions
      DOUBLE PRECISION primea
      integer ifloor
c ========INCLUDE HEADERS==========================
      include 'trig.h'
c =================================================
c Prediction of $\lambda$ at time tp, starting from elements at time t0
      en1=sqrt(gms/eq0(1)**3)
      pl12=eq0(6)+en1*(tp-t0)
      ng0=ifloor(pl12/dpig)
      pl12p=pl12-ng0*dpig
c The number of revolutions is computed to reduce the discrepancies in 
c lambda between the two arcs
      if(eqp(6)-pl12p.gt.pig)then
         ng0=ng0-1
      elseif(pl12p-eqp(6).gt.pig)then
         ng0=ng0+1
      endif
c Prediction of $\lambda$ at time t0, starting from elements at time tp
      en2=sqrt(gms/eqp(1)**3)
      pl21=eqp(6)+en2*(t0-tp)
      ngp=ifloor(pl21/dpig)
      pl21p=pl21-ngp*dpig
c The number of revolutions is computed to reduce the discrepancies in
c $\lambda$ \hfil\break
c between the two arcs
      if(eq0(6)-pl21p.gt.pig)then
         ngp=ngp-1
      elseif(pl21p-eq0(6).gt.pig)then
         ngp=ngp+1
      endif
c Computation of mean mean motion using the 2-body model, starting
c firstly from arc 1 and then from arc 2
      enm0=(eqp(6)-eq0(6)+ng0*dpig)/(tp-t0)
      enmp=(eq0(6)-eqp(6)+ngp*dpig)/(t0-tp)
c Choose  the best orbit to start from
      write(*, 177)ng0,ngp
 177  format(' number of rev. ',i4, i4) 
      if(iorb.eq.1)then
         enm=enm0
         ng=ng0
      elseif(iorb.eq.2)then
          enm=enmp
          ng=-ngp
       else
          write(*,*)' start: iorb=',iorb
          stop
       endif
      am=(gms/enm**2)**(1.d0/3.d0)
c 1997 change:use estimated mean motion to compute mean
c longitude at time in the middle
      tm=(t0+tp)/2.d0
      pll1=eq0(6)+enm*(tm-t0)
      pll2=eqp(6)+enm*(tm-tp)
      plm=primea(pll1,pll2)
      return
      end
c =====================================================================
      INTEGER FUNCTION ifloor(a)
      double precision a
      ifloor=a
      if(a.lt.0)ifloor=ifloor-1
      return
      end


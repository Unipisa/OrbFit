c ============ephem_prop===================
c PUBLIC ROUTINES:
c                  fstpro
c                  proele
c                  proelc
c                  fsteph
c                  ephemc
c INTERNAL 
c ROUTINES
c                  srtept
c                  outco
c
c  HEADERS
c ephem_prop.o: phase.h
c
c                 trig.h, sunmass.h  public
c
c ====================================================
c FSTPRO state propagation for FITOBS
c ====================================================
      SUBROUTINE fstpro(batch,icov,ini0,cov0,iun20,iun8,ok,
     +         t0,eq0,h0,g0,c0,tr,eq1,h1,g1,c1)
      IMPLICIT NONE
c =================INPUT=========================================
c requirements on covariance
      INTEGER icov
c batch mode
      LOGICAL batch
c availability of initial conditions, covariance, all required
      LOGICAL ini0,cov0,ok
c output units
      INTEGER iun20,iun8
c epoch time, target time
      DOUBLE PRECISION t0,tr
c elements, magnitude, covariance and normal matrix
      DOUBLE PRECISION eq0(6),h0,g0(6,6),c0(6,6)
c ================OUTPUT=================================
c elements, magnitude, covariance and normal matrix at epoch tr
      DOUBLE PRECISION eq1(6),h1,g1(6,6),c1(6,6)
c ================END INTERFACE==========================
c trigonometric constants
      INCLUDE 'trig.h'
c keplerian elements, mena motion, opposition effect
      DOUBLE PRECISION ekr(6),enne
c     DOUBLE PRECISION gmag
c loop indexes
      INTEGER ii,j
c ======== constant of gravitation ==============
      INCLUDE 'sunmass.h'
c =====================================================================
c check availability of required data
      IF(batch)THEN
         CALL chereq(icov,ini0,cov0,t0,iun20,-1,ok)
      ELSE
         CALL chereq(icov,ini0,cov0,t0,iun20,iun8,ok)
      ENDIF
      IF(.not.ok)RETURN
c =====================================================================
c check availability of JPL ephemrides
         CALL chetim(t0,tr,ok)
         IF(.not.ok) THEN
            WRITE(*,*)' JPL ephemrides not available for tr=',tr
            ok=.false.
            RETURN
         ENDIF
c =====================================================================
c propagation to time tr
      IF(icov.eq.1)THEN
c state vector only
         CALL proele('EQU',t0,eq0,tr,eq1)
         cov0=.false.
      ELSEIF(icov.eq.2)THEN
         CALL proelc('EQU',t0,eq0,g0,c0,tr,eq1,g1,c1)
      ENDIF 
      h1=h0
c =====================================================================
c output
c =====================================================================
      IF(.not.batch.and.iun20.gt.0)THEN
         WRITE(iun20,223) tr
         WRITE(*,223) tr
 223     FORMAT(' elements at time ',f8.1,' (MJD):')
         WRITE(*,105) eq1
         WRITE(iun20,105) eq1
         CALL coocha(eq1,'EQU',gms,ekr,'KEP',enne)
         DO  ii=3,6
            ekr(ii)=ekr(ii)*degrad
         ENDDO
         WRITE(*,105)ekr
         WRITE(iun20,105)(ekr(j),j=6,1,-1)
         WRITE(*,*)
         WRITE(iun20,*)
 105     FORMAT(6f13.7)
         IF(icov.eq.2.and.iun8.gt.0)THEN
            WRITE(iun8,*) 'COVARIANCE MATRIX FOR NEW EPOCH'
            WRITE(iun8,223) tr
            CALL outco(iun8,g1,c1)
         ENDIF
      ENDIF
      RETURN
      END
c =====================================================================
c PROELE-PROELC 
c =====================================================================
c interface to N+1-body problem propagator (subroutine propag)
c
c This version uses JPL ephemerides as source for the planetary positions
c
c  WARNING: the input and output elements east0, east1 
c           are in the ecliptic (mean of J2000.0) system
c
c  WARNING: input and output use the same coo. To transform
c           also the coordinate type, combine with coocha/cooder
c =====================================================================
c PROELE- version with elements only, no covariance
c =====================================================================
c
c  input: coo coordinate type EQU, KEP, CAR, EQP
c         t0 epoch time (MJD)
c         t1 prediction time (MJD)
c         east0  orbital elements vector  at time t0
c  output:
c         east1  orbital elements for the asteroid at time t1
c ============INTERFACE===================================================
      SUBROUTINE proele(coo,t0,east0,t1,east1)
      implicit none
c equinoctal elements and epoch times
      character*3 coo
      double precision east0(6),east1(6),t0,t1
c ============END INTERFACE===============================================
c cartesian coordinates, mean motion, derivatives(not used)
      double precision xastr(6),xear(6),enne,dxde(6,6),ddxde(3,6,6)
c equinoctal elements
      double precision eq(6)
c ======== constant of gravitation ==============
      include 'sunmass.h'
c****************
c   static memory not required
c****************
c =====================================================================
c change coordinates
      call  coocha(east0,coo,gms,eq,'EQU',enne)
c =====================================================================
c call private propagation routine, derivatives not required
      call propag(t0,eq,t1,xastr,xear,0,dxde,ddxde)
c compute new elements
      call coocha(xastr,'CAR',gms,east1,coo,enne)
      return
      end
c =====================================================================
c PROELC- version with elements and covariance propagation
c =====================================================================
c
c  input: coo coordinate type EQU, KEP, CAR, EQP
c         t0 epoch time (MJD)
c         t1 prediction time (MJD)
c         east0 orbital elements vector at time t0
c         gamma0 covariance matrix for the elements east0 at epoch t0
c         c0 normal matrix for the same
c  output:
c         east1 orbital elements for the asteroid at time t1
c         gamma1 covariance matrix for the elements east1 at epoch t1
c         c1 normal matrix for the same 
c ============INTERFACE===================================================
      SUBROUTINE proelc(coo,t0,east0,gamma0,c0,t1,east1,gamma1,c1)
      implicit none
c elements and epoch times
      character*3 coo
      double precision east0(6),east1(6),t0,t1
c normal matrices
      DOUBLE PRECISION c0(6,6),c1(6,6)
c covariance matrices
      double precision gamma0(6,6),gamma1(6,6)
c ============END INTERFACE===============================================
c equinoctal elements
      double precision eq(6)
c cartesian coordinates, mean motion, derivatives
      double precision xastr(6),xear(6),enne,dxde(6,6),ddxde(3,6,6)
c derivatives of elements w.r. to cartesian, ele w.r. elements, workspaces
      double precision dedx(6,6),dede(6,6),dedet(6,6),tmp(6,6),der0(6,6)
      DOUBLE PRECISION dedein(6,6),det,dedint(6,6)
      INTEGER ising,invop
c ======== constant of gravitation ==============
      include 'sunmass.h'
c****************
c   static memory not required
c****************
c =====================================================================
c     change coordinates; der0=d(eq)/d(east0)
      call  cooder(east0,coo,gms,eq,'EQU',enne,der0)
c =====================================================================
c call private propagation routine, requiring derivatives
      call propag(t0,eq,t1,xastr,xear,1,dxde,ddxde)
c compute new elements, wtih partial derivatives
      call cooder(xastr,'CAR',gms,east1,coo,enne,dedx)
c chain rule to obtain d(east1)/d(east0)
      call mulmat(dedx,6,6,dxde,6,6,tmp)
      call mulmat(tmp,6,6,der0,6,6,dede)
c covariance matrix is propagated by similarity transformation
      call mulmat(dede,6,6,gamma0,6,6,tmp)
      call transp(dede,6,6,dedet)
      call mulmat(tmp,6,6,dedet,6,6,gamma1)
c normal matrix is propagated by similarity transformation
      call mcopy(6,6,dede,dedein)
      invop=1
      CALL matin(dedein,det,6,0,6,ising,invop)
      IF(ising.ne.0)THEN
         WRITE(*,*)'PROELC: this should not happen'
      ENDIF
      CALL mulmat(c0,6,6,dedein,6,6,tmp)
      CALL transp(dedein,6,6,dedint)
      CALL mulmat(dedint,6,6,tmp,6,6,c1)
      return
      end
c Copyright 1998 Orbfit Consortium
c Version December 18, 1998 Steven Chesley (chesley@dm.unipi.it)
c ====================================================
c FSTEPH Ephemerides Generation for FITOBS
c This is a hacked version of FSTPRO
c First integrate from t0 back to tr, saving at interval step.
c Then write this data in chronological order
c Finally integrate from to forward to tf, writing after each step
c inputs:
c     name - asteroid name (for file and ephem listing)
c     dir - location for files
c     defele - logical for existence of elements
c     ok - success flag
c     coord - coordinate type for input elements(EQU,CAR,KEP)
c     t0 - current epoch
c     eq0 - input elements
c     tr - initial time for ephemerides
c     tf - final time
c     step - stepsize
c     numsav - an upper bound for the number of steps necessary before t0
c     ephefl - logical flag for output of ephemrides file
c     eltype - coordinate type for ephemerides (output) elements (EQU,CAR,KEP)
c     moidfl - logical flag for output of file of MOID's and nodal distances
c output: NONE
c     
c REMARK: IF moidfl=.false. and ephefl=.false. you will still 
c     get a close approach file.
c ====================================================
      SUBROUTINE fsteph(name,dir,defele,ok,coord,t0,eq0,
     +     tr,tf,step,numsav,ephefl,eltype,moidfl)
      IMPLICIT NONE
c =================INPUT=========================================
c name, place, output element type
      CHARACTER*(*) name,dir
c availability of initial conditions, covariance, all required
      LOGICAL defele,ok,moidfl,ephefl
c necessary size of temporary storage array
      INTEGER numsav
c epoch times, elements
      DOUBLE PRECISION t0,tr,tf,step,eq0(6)
      CHARACTER*(3) coord,eltype
c ================OUTPUT=================================
c none
c ================END INTERFACE==========================
c loop indexes
      INTEGER i,n
c converted elements
      DOUBLE PRECISION elem0(6),enne,eqtmp(6)
c temporary storage
      INTEGER numsavx
      PARAMETER (numsavx=10000)
      DOUBLE PRECISION tsav(numsavx),eqsav(6,numsavx),t1,elem1(6),t2
      INTEGER jst
c output
      INTEGER unit
      INTEGER ln,lnm,lnnam
      CHARACTER*(60) file,filem
c ======== time spans for JPL data etc. =========
      INCLUDE 'sunmass.h'
c moid
      double precision moid,dnp,dnm
      double precision msav(numsavx),ndsav(numsavx,2)
      integer munit,iconv,icsav(numsavx)
c =====================================================================
c check dimensions
      IF(numsav.gt.numsavx)THEN
         WRITE(*,*)'fsteph: numsav=',numsav,' is > numsavx=', numsavx
         RETURN
      ENDIF
c check availability of required data
      ok=.true.
      IF (.not.defele)THEN
         WRITE(*,*) 'Sorry: You must provide an orbit first.'
         ok = .false.
      ENDIF
      IF (tr .ge. tf)THEN
         WRITE(*,*) 'Sorry: Initial time must be before final time.'
         ok = .false.
      ENDIF
      IF(.not.ok)RETURN
c =====================================================================
c check availability of JPL ephemrides
      call chetim(tr,tf,ok)
      if(.not.ok)return
c========= OPEN FILES =================================
      if(ephefl)then
c open ephemerides ouput file
         call filnam(dir,name,'ele',file,ln)
         call filopn(unit,file(1:ln),'unknown')
         call wro1lh(unit,'ECLM','J2000',eltype)
      endif
      if(moidfl)then
c open moid file
         call filnam(dir,name,'moid',filem,lnm)
         call filopn(munit,filem(1:lnm),'unknown')
      endif
      call rmsp(name,lnnam)
c========= CONVERT INPUT ELEMENTS TO REQUESTED TYPE ===============
      call coocha(eq0,coord,gms,elem0,eltype,enne)
c========= PROPAGATE, SAVE, AND THEN WRITE IN TIME ORDER ===============
c     step back in time first (if necessary)
      CALL set_restart(.true.)
      if(tr .le. t0)then
         t1=t0
         do i=1,6
            elem1(i)=elem0(i)
         enddo
c        but first propagate back to tf if t0 isn't within interval
         if(tf .lt. t0)then
            t1=tf
            call proele(eltype,t0,elem0,t1,elem1)
         endif
         if(moidfl)then
            call coocha(elem1,eltype,gms,eqtmp,'EQU',enne)
            call nomoid(t1,eqtmp,moid,iconv,dnp,dnm)
         endif
         n=1
         tsav(n)=t1
         msav(n)=moid
         icsav(n)=iconv
         ndsav(n,1)=dnp
         ndsav(n,2)=dnm
         do i=1,6
            eqsav(i,n)=elem1(i)
         enddo
c        start loop
         DO jst=1,10000
            t2=t1-jst*step
            IF(t2.lt.tr) GOTO 5
c           WRITE(*,*)t2
c         do t1=t1,tr+step,-step
            call proele(eltype,t0,elem0,t2,elem1)
            CALL set_restart(.false.)
            if(moidfl)then
               call coocha(elem1,eltype,gms,eqtmp,'EQU',enne)
               call nomoid(t1-step,eqtmp,moid,iconv,dnp,dnm)
            endif
            n=n+1
            tsav(n)=t2
            msav(n)=moid
            ndsav(n,1)=dnp
            ndsav(n,2)=dnm
            icsav(n)=iconv
            do i=1,6
               eqsav(i,n)=elem1(i)
            enddo
         enddo
c        Get last point
 5       CONTINUE
         if(tsav(n).gt.tr)then
            t2=tr
c           WRITE(*,*)t2
            call proele(eltype,t0,elem0,t2,elem1)
            if(moidfl)then
               call coocha(elem1,eltype,gms,eqtmp,'EQU',enne)
               call nomoid(tr,eqtmp,moid,iconv,dnp,dnm)
            endif
            n=n+1
            tsav(n)=tr
            msav(n)=moid
            ndsav(n,1)=dnp
            ndsav(n,2)=dnm
            icsav(n)=iconv
            do i=1,6
               eqsav(i,n)=elem1(i)
            enddo
         endif
         CALL set_restart(.true.)
c        write in forward time order
c MAGNITUDE IS WRONG!!!!
         do i=n,1,-1
            if(ephefl) call wro1lr(unit,name(1:lnnam),eqsav(1,i),
     +           eltype,tsav(i),-1.d4,1.0d0)
            if(moidfl) write(munit,199)tsav(i),msav(i),
     +           ndsav(i,1),ndsav(i,2)
 199        format(f13.4,1x,f8.5,1x,f8.5,1x,f8.5)
         enddo
      endif
c========= PROPAGATE AND WRITE SIMULTANEOUSLY ========================
      if(tf .ge. t0)then
         t1=t0
         do i=1,6
            elem1(i)=elem0(i)
         enddo
c        now step to the beginning (if necessary)
         if(tr .gt. t0)then
            t1=tr
            call proele(eltype,t0,elem0,t1,elem1)
            if(ephefl) call wro1lr(unit,name(1:lnnam),elem1,eltype,
     +           t1,-1.d4,1.0d0)
            if(moidfl)then
               call coocha(elem1,eltype,gms,eqtmp,'EQU',enne)
               call nomoid(t1,eqtmp,moid,iconv,dnp,dnm)
               write(munit,199)t1,moid,dnp,dnm
            endif
         endif
c         do t1=t1,tf-step,step
         DO jst=1,100000
            t2=t1+jst*step
            IF(t2.gt.tf) GOTO 6
c           WRITE(*,*)t2
            call proele(eltype,t0,elem0,t2,elem1)
            CALL set_restart(.false.)
            if(ephefl) call wro1lr(unit,name(1:lnnam),elem1,eltype,
     +           t2,-1.d4,1.0d0)
            if(moidfl)then
               call coocha(elem1,eltype,gms,eqtmp,'EQU',enne)
               call nomoid(t1+step,eqtmp,moid,iconv,dnp,dnm)
               write(munit,199)t1+step,moid,dnp,dnm
            endif
         enddo
 6       CONTINUE
c        Get last point
         if(t2.lt.tf)then
            t2=tf
c           WRITE(*,*)t2
            call proele(eltype,t0,elem0,tf,elem1)
            if(ephefl) call wro1lr(unit,name(1:lnnam),elem1,eltype,
     +           tf,-1.d4,1.0d0)
            if(moidfl)then
               call coocha(elem1,eltype,gms,eqtmp,'EQU',enne)
               call nomoid(tf,eqtmp,moid,iconv,dnp,dnm)
               write(munit,199)tf,moid,dnp,dnm
            endif
         endif
         CALL set_restart(.true.)
      endif
c =====================================================================
      if(ephefl) call filclo(unit,' ')
      if(moidfl) call filclo(munit,' ')
      if(ephefl) write(*,*)'Ephemeris file generated:',file(1:ln)
      if(moidfl) write(*,*)'Moid file generated:',filem(1:lnm)

      RETURN
      END
* Copyright (C) 1997-2000 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: December 11, 2000
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         E P H E M C                           *
*  *                                                               *
*  *                 Computation of ephemerides                    *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    UNIT      -  Output FORTRAN unit
*           ELTYPE    -  Type of orbital elements (EQU/KEP/CAR)
*           T0        -  Epoch of orbital elements (MJD, TDT)
*           ELEM      -  Orbital elements (ECLM J2000)
*           COVE      -  Covariance matrix of orbital elements
*           DEFCOV    -  Tells whether the covariance matrix is defined
*           T1        -  Starting time for ephemeris (MJD, TDT)
*           T2        -  Ending time for ephemeris (MJD, TDT)
*           DT        -  Ephemeris stepsize (d)
*           MASS      -  Mass (solar masses)
*           HMAG      -  H absolute magnitude (if <-100, missing)
*           GMAG      -  G slope parameter
*           IDSTA     -  Station identifier
*           SCALE     -  Timescale for output
*           FIELDS    -  Output fields (separated by commas)
*
* SUPPORTED OUTPUT FIELDS:
*     cal       calendar date
*     mjd       Modified Julian Day
*     coord     coordinates (RA and DEC, or ecliptic long. and lat.)
*     mag       magnitude
*     delta     distance from the Earth
*     r         distance from the Sun
*     elong     Sun elongation angle
*     phase     Sun phase angle
*     glat      galactic latitude
*     appmot    apparent motion
*     skyerr    sky plane error
*
      SUBROUTINE ephemc(unit,eltype,t0,elem,cove,defcov,t1,t2,dt,mass,
     +                  hmag,gmag,idsta,scale,fields)
      IMPLICIT NONE

      INTEGER unit,idsta
      CHARACTER*(*) eltype,scale,fields
      DOUBLE PRECISION t0,t1,t2,elem(6),dt,mass,hmag,gmag,cove(6,6)
      LOGICAL defcov

      DOUBLE PRECISION epst
      PARAMETER (epst=1.d-6)

* Max number of output fields
      INTEGER nfx
      PARAMETER (nfx=30)
* Max number of ephemeris points
      INTEGER nephx
      PARAMETER (nephx=600)
* Max length of output records
      INTEGER lrx
      PARAMETER (lrx=200)

      INTEGER nf,lh,ider,i,lf,lr,day,month,year,ia,ma,id,md,ls,neph,k,ip
      INTEGER srtord(nephx),lrv(nephx),obstyp,iepfor,la,lad,lfo,lf1,lf2
      INTEGER inb1,inb2
      PARAMETER (obstyp=1000)
      DOUBLE PRECISION tdt,alpha,delta,mag,hour,sa,sd,difft,signdt,cvf
      DOUBLE PRECISION gamad(2,2),sig(2),axes(2,2),err1,err2,pa
      DOUBLE PRECISION teph(nephx),velsiz
      CHARACTER*1 siga,sigd,anguni
      CHARACTER*3 cmonth(12)
      CHARACTER*20 field(nfx),frameo,amuni,amunid,amfor
      CHARACTER*(lrx) head1,head2,head3,head4,outrec,recv(nephx),blank
      CHARACTER cval*80
      LOGICAL outmot,outerr,outmag,oldrst,found,fail1,fail,usexp
* Time conversion (added by Steve Chesley)
      INTEGER mjdt,mjdout
      DOUBLE PRECISION sect,secout,tout

      INTEGER lench,intlo
      EXTERNAL lench,intlo

      INCLUDE 'trig.h'
      INCLUDE 'phase.h'

      DATA cmonth/'Jan','Feb','Mar','Apr','May','Jun',
     +            'Jul','Aug','Sep','Oct','Nov','Dec'/

      inb1=3
      inb2=2

* Parameter which should become options
      frameo='EQUATORIAL'

      signdt=1
      IF(dt.LT.0) signdt=-1

* List of ephemeris epochs
      neph=0
      tdt=t1
 2    CONTINUE
      difft=signdt*(tdt-t2)
      IF(difft.LE.epst) THEN
          neph=neph+1
          IF(neph.GT.nephx) STOP '**** ephemc: neph > nephx ****'
          teph(neph)=tdt
          tdt=t1+neph*dt
          GOTO 2
      END IF

* Sorting of ephemeris epochs
      CALL srtept(teph,neph,t0,srtord)

* List of output fields
      CALL spflds(fields,field,nf,nfx)

* COMPOSITION OF HEADER LINES
      head1=' '
      head2=' '
      head3=' '
      head4=' '
      blank=' '
      lh=0
      outmot=.false.
      outerr=.false.
      ider=0

      DO 5 i=1,nf
      IF(field(i).EQ.'cal') THEN
          head1(lh+1:)=blank
          head2(lh+1:)='    Date      Hour '
          ls=lench(scale)
          WRITE(head3(lh+1:),300) scale(1:ls)
          head4(lh+1:)=' =========== ======'
          lh=lh+19
      ELSEIF(field(i).EQ.'mjd') THEN
          head1(lh+1:)=blank
          head2(lh+1:)='     MJD     '
          ls=lench(scale)
          WRITE(head3(lh+1:),301) scale(1:ls)
          head4(lh+1:)=' ============'
          lh=lh+13
      ELSEIF(field(i).EQ.'coord') THEN
          IF(frameo.EQ.'EQUATORIAL') THEN
              head1(lh+1:)='     Equatorial coordinates  '
              head2(lh+1:)='       RA            DEC     '
          ELSEIF(frameo.EQ.'ECLIPTICAL') THEN
              head1(lh+1:)='      Ecliptic coordinates   '
              head2(lh+1:)='    Longitude      Latitude  '
          ELSE
              lf=lench(frameo)
              WRITE(*,331) frameo(1:lf)
              STOP '**** ephemc: Abnormal end ****'
          END IF
          head3(lh+1:)='    h  m  s        d  ''  "   '
          head4(lh+1:)='  =============  ============'
          lh=lh+29
      ELSEIF(field(i).EQ.'delta') THEN
          head1(lh+1:)=blank
          head2(lh+1:)='  Delta '
          head3(lh+1:)='   (AU) '
          head4(lh+1:)=' ======='
          lh=lh+8
      ELSEIF(field(i).EQ.'r') THEN
          head1(lh+1:)=blank
          head2(lh+1:)='    R   '
          head3(lh+1:)='   (AU) '
          head4(lh+1:)=' ======='
          lh=lh+8
      ELSEIF(field(i).EQ.'elong') THEN
          head1(lh+1:)=blank
          head2(lh+1:)=' Elong'
          head3(lh+1:)=' (deg)'
          head4(lh+1:)=' ====='
          lh=lh+6
      ELSEIF(field(i).EQ.'phase') THEN
          head1(lh+1:)=blank
          head2(lh+1:)=' Phase'
          head3(lh+1:)=' (deg)'
          head4(lh+1:)=' ====='
          lh=lh+6
      ELSEIF(field(i).EQ.'mag') THEN
          outmag=(hmag.GT.-100.d0)
          IF(outmag) THEN
              head1(lh+1:)=blank
              head2(lh+1:)='  Mag '
              head3(lh+1:)='      '
              head4(lh+1:)=' ====='
              lh=lh+6
          END IF
      ELSEIF(field(i).EQ.'glat') THEN
          head1(lh+1:)=blank
          head2(lh+1:)=' Glat '
          head3(lh+1:)=' (deg)'
          head4(lh+1:)=' ====='
          lh=lh+6
      ELSEIF(field(i).EQ.'skyerr') THEN
          IF(defcov) THEN
              outerr=.true.
              ider=1
              head1(lh+1:)=blank
              head2(lh+1:)='       Sky plane error    '
              head3(lh+1:)='     Err1      Err2    PA '
              head4(lh+1:)='  ========  ======== ====='
              lh=lh+26
          END IF
      ELSEIF(field(i).EQ.'appmot') THEN
          fail=.false.
* Default format for apparent motion:
* IEPFOR = 1 -> (Vx, Vy)
* IEPFOR = 2 -> (V, PosAng)
          iepfor=1
          CALL sv2int('ephem.appmot.','format',cval,iepfor,.false.,
     +                found,fail1,fail)
          amunid='"/min'
          amuni=amunid
          CALL rdncha('ephem.appmot.','units',amuni,.false.,
     +                found,fail1,fail)
          IF(fail) STOP '**** ephemc: abnormal end ****'
* CVF = conversion factor from rad/d to selected unit
          CALL angvcf(amuni,cvf,fail)
          IF(fail) THEN
              la=lench(amuni)
              lad=lench(amunid)
              WRITE(*,320) amuni(1:la),amunid(1:lad)
              amuni=amunid
              CALL angvcf(amuni,cvf,fail)
              IF(fail) STOP '**** ephemc: internal error (01) ****'
          END IF
          la=lench(amuni)
* Choice of the output format: normally F format is preferred
          amfor='(F10.4)'
* LFO = length of the output string
          lfo=10
          usexp=.false.
* Check whether F format can supply the required dynamic range
          IF(1.D0*cvf.GE.500.D0) usexp=.true.
          IF(3.D-3*cvf.LE.0.1D0) usexp=.true.
* Otherwise use exponential format
          IF(usexp) THEN
              amfor='(1P,E12.4)'
              lfo=12
          END IF
          lf1=lfo
          IF(iepfor.EQ.1) THEN
              lf2=lfo
          ELSE
              lf2=6
          END IF
          CALL filstr('App. motion',cval,lf1+lf2,inb1,0)
          head1(lh+1:)=cval
          IF(iepfor.EQ.1) THEN
              CALL filstr('RA',cval,lf1,inb1,0)
          ELSE
              CALL filstr('Vel',cval,lf1,inb1,0)
          END IF
          head2(lh+1:)=cval
          CALL filstr(amuni,cval,lf1,inb1,0)
          head3(lh+1:)=cval
          head4(lh+1:)='  =================='
          lh=lh+lf1
          IF(iepfor.EQ.1) THEN
              CALL filstr('DEC',cval,lf2,inb1,0)
              head2(lh+1:)=cval
              CALL filstr(amuni,cval,lf2,inb1,0)
              head3(lh+1:)=cval
              head4(lh+1:)='  =================='
          ELSE
              CALL filstr('PA',cval,lf2,inb2,0)
              head2(lh+1:)=cval
              CALL filstr('deg',cval,lf2,inb2,0)
              head3(lh+1:)=cval
              head4(lh+1:)=' ====='
          END IF
          lh=lh+lf2
          outmot=.true.
      ELSE
          lf=lench(field(i))
          WRITE(*,330) field(i)(1:lf)
          STOP '**** ephemc: abnormal end ****'
      END IF
      IF(lh.GT.lrx) STOP '**** ephemc: lh > lrx ****'
 5    CONTINUE
 300  FORMAT('             (',A,') ')
 301  FORMAT('    (',A,')    ')
 330  FORMAT('Sorry, I don''t know how to produce field "',A,'"')
 331  FORMAT('Sorry, I don''t know "',A,'" reference system')
 320  FORMAT('WARNING(ephemc): I do not know how to use "',A,
     +    '" as units for apparent motion;'/
     +    17X,'using instead default units "',A,'"')
      WRITE(unit,100) head1(1:lh)
      WRITE(unit,100) head2(1:lh)
      WRITE(unit,100) head3(1:lh)
      WRITE(unit,100) head4(1:lh)
 100  FORMAT(A)

      CALL set_restart(.true.)

* Start loop on ephemeris epochs
      DO 1 k=1,neph
      ip=srtord(k)
      tdt=teph(ip)
* Numerical integration
      IF(outerr) THEN
          CALL preobc(eltype,t0,idsta,tdt,elem,hmag,gmag,
     +                cove,obstyp,alpha,delta,mag,gamad,sig,axes)
      ELSE
          CALL preobs(eltype,t0,idsta,tdt,elem,obstyp,alpha,delta,
     +                hmag,gmag,mag)
      END IF
      CALL set_restart(.false.)

* COMPOSITION OF OUTPUT RECORD
* Time scale conversion
      mjdt=intlo(tdt)
      sect=(tdt-mjdt)*86400.d0
      CALL cnvtim(mjdt,sect,'TDT',mjdout,secout,scale)
      tout=secout/86400.d0+mjdout
*
      outrec=' '
      lr=0
      DO 6 i=1,nf
* Calendar date
      IF(field(i).EQ.'cal') THEN
          CALL mjddat(tout,day,month,year,hour)
          IF(month.LT.1 .OR. month.GT.12)
     +        STOP '**** ephemc: internal error (02) ****'
          WRITE(outrec(lr+1:),201) day,cmonth(month),year,hour
          lr=lr+19
* Modified Julian Day
      ELSEIF(field(i).EQ.'mjd') THEN
          WRITE(outrec(lr+1:),202) tout
          lr=lr+13
* Astrometric coordinates
      ELSEIF(field(i).EQ.'coord') THEN
          CALL sessag(alpha*hrad,siga,ia,ma,sa)
          IF(siga.NE.'+') STOP '**** ephemc: internal error (03) ****'
          CALL sessag(delta*degrad,sigd,id,md,sd)
          WRITE(outrec(lr+1:),203) ia,ma,sa,sigd,id,md,sd
          lr=lr+29
* Distance from the Earth
      ELSEIF(field(i).EQ.'delta') THEN
          WRITE(outrec(lr+1:),204) dis
          lr=lr+8
* Distance from the Sun
      ELSEIF(field(i).EQ.'r') THEN
          WRITE(outrec(lr+1:),204) dsun
          lr=lr+8
* Solar elongation
      ELSEIF(field(i).EQ.'elong') THEN
          WRITE(outrec(lr+1:),205) elo*degrad
          lr=lr+6
* Solar phase angle
      ELSEIF(field(i).EQ.'phase') THEN
          WRITE(outrec(lr+1:),205) pha*degrad
          lr=lr+6
* Magnitude
      ELSEIF(field(i).EQ.'mag') THEN
          IF(outmag) THEN
              WRITE(outrec(lr+1:),205) mag
              lr=lr+6
          END IF
* Sky plane error
      ELSEIF(field(i).EQ.'skyerr') THEN
          IF(outerr) THEN
              err1=sig(1)*degrad
              err2=sig(2)*degrad
* correction A. Milani 19/3/2000: remember right ascension increases
* from right to left
*             pa=ATAN2(axes(2,1),axes(1,1))
              pa=ATAN2(axes(2,1),-axes(1,1))
              anguni='d'
              IF(MAX(err1,err2).LT.1.d0) THEN
                  err1=err1*60
                  err2=err2*60
                  anguni=''''
                  IF(MAX(err1,err2).LT.1.d0) THEN
                      err1=err1*60
                      err2=err2*60
                      anguni='"'
                  END IF
              END IF
              WRITE(outrec(lr+1:),208) err2,anguni,err1,anguni,pa*degrad
              lr=lr+26
          END IF
* Galactic latitude
       ELSEIF(field(i).EQ.'glat') THEN
          WRITE(outrec(lr+1:),205) gallat*degrad
          lr=lr+6
* Apparent motion
      ELSEIF(field(i).EQ.'appmot') THEN
          IF(iepfor.EQ.1) THEN
              WRITE(outrec(lr+1:),amfor) adot*cvf
              lr=lr+lf1
              WRITE(outrec(lr+1:),amfor) ddot*cvf
              lr=lr+lf2
          ELSE
              velsiz=SQRT(adot**2+ddot**2)
              WRITE(outrec(lr+1:),amfor) velsiz*cvf
              lr=lr+lf1
              pa=ATAN2(adot,ddot)
              IF(pa.LT.0.D0) pa=pa+dpig
              WRITE(outrec(lr+1:),205) pa*degrad
              lr=lr+6
          END IF
      ELSE
          lf=lench(field(i))
          WRITE(*,340) field(i)(1:lf)
          STOP '**** ephemc: internal error (04) ****'
      END IF
 6    CONTINUE
 201  FORMAT(I3,1X,A3,I5,F7.3)
 202  FORMAT(F13.6)
 203  FORMAT(2X,2I3,F7.3,2X,A1,I2,I3,F6.2)
 204  FORMAT(F8.4)
 205  FORMAT(F6.1)
 207  FORMAT(2F8.4)
 208  FORMAT(2(F9.3,A1),F6.1)
 340  FORMAT(' ERROR: illegal output field "',A,'"')
      IF(lr.GT.lrx) STOP '**** ephemc: lr > lrx ****'
      recv(ip)=outrec(1:lr)
      lrv(ip)=lr
 1    CONTINUE

      CALL set_restart(.true.)


* Output of ephemeris records in the required order
      DO 3 ip=1,neph
      lr=lrv(ip)
      WRITE(unit,100) recv(ip)(1:lr)
 3    CONTINUE

      END
* Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: December 12, 1997
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         S R T E P T                           *
*  *                                                               *
*  *                 Sorting of ephemeris epochs                   *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    T         -  Ephemeris epochs
*           N         -  Number of ephemeris epochs
*           T0        -  Epoch of initial conditions
*
* OUTPUT:   IPT       -  Sorting pointer
*
      SUBROUTINE srtept(t,n,t0,ipt)
      IMPLICIT NONE

      INTEGER n,ipt(n)
      DOUBLE PRECISION t(n),t0

      INTEGER i,i1,i2,k1,k2,nit
      LOGICAL ge1,ge2,rev,change

* Initial guess for IPT (assuming ephemeris times to be sorted
* in ascending order)
      i1=n+1
      DO 1 i=1,n
      IF(t(i).GE.t0) THEN
          i1=i
          GOTO 2
      END IF
 1    CONTINUE
 2    CONTINUE

* I1 is now the number of the first t(i)>=t0
      i2=0
      DO 3 i=i1,n
      i2=i2+1
      ipt(i2)=i
 3    CONTINUE
      DO 4 i=i1-1,1,-1
      i2=i2+1
      ipt(i2)=i
 4    CONTINUE
      IF(i2.NE.n) STOP '**** srtept: internal error (01) ****'

* Sorting
      nit=0
 5    CONTINUE
      nit=nit+1
      IF(nit.GT.n+3) STOP '**** srtept: internal error (02) ****'
      change=.false.
      DO 6 k1=1,n-1
      k2=k1+1
      i1=ipt(k1)
      i2=ipt(k2)
      ge1=(t(i1).GE.t0)
      ge2=(t(i2).GE.t0)
      IF(ge1) THEN
          IF(ge2) THEN
              rev=(t(i1).GT.t(i2))
          ELSE
              rev=.false.
          END IF
      ELSE
          IF(ge2) THEN
              rev=.true.
          ELSE
              rev=(t(i1).LT.t(i2))
          END IF
      END IF
      IF(rev) THEN
          i=ipt(k1)
          ipt(k1)=ipt(k2)
          ipt(k2)=i
          change=.true.
      END IF
 6    CONTINUE
      IF(change) GOTO 5

      END
c =====================================================================
c OUTCO
c =====================================================================
c  output of covariance, computation of eigenvalues
c   input: iun   = output unit
c          gamma = covariance matrix
c          c     = normal matrix
      SUBROUTINE outco(iun,gamma,c)
      implicit none
c output unit, error flag
      integer iun,ierr
c loop indexes
      integer i,j
c covariance, normal matrix
      double precision gamma(6,6),c(6,6)
c eigenvalues, eigenvectors
      double precision eigvec(6,6),eigval(6),fv1(6),fv2(6)
c output covariance
      write(iun,*)
      write(iun,*) 'COVARIANCE MATRIX'
      do j=1,6
        write(iun,109) (gamma(i,j),i=1,6)
 109    format(6e24.16)
      enddo
c eigenvalues
      call rs(6,6,gamma,eigval,1,eigvec,fv1,fv2,ierr)
      write(iun,*)
      write(iun,*) 'EIGENVALUES '
      write(iun,109) (eigval(i),i=1,6)
      write(iun,*)
      write(iun,*) 'EIGENVECTORS'
      do 3 j=1,6
c by columns (check)
        write(iun,109) (eigvec(i,j),i=1,6)
 3    continue
      write(iun,*)
      write(iun,*) 'NORMAL MATRIX'
      do j=1,6
        write(iun,109) (c(i,j),i=1,6)
      enddo
      return 
      end

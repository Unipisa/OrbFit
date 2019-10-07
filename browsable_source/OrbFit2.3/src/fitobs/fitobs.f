c Copyright (C) 1997,1998 ORBFIT consortium
c Version 1.6.0: A. Milani and Z. Knezevic, March 23, 1998
c patch 1.6.1: A. Milani, May 2, 1998
c Version 1.7.0: A. Milani, July 8, 1998
c patch 1.7.2: modified by Z. Knezevic, November 18, 1998
c Version 1.8.0: S. Chesley and A. Milani, December 15, 1998
c Version 1.9.0: A. Milani and S. Chesley, February 1999
c Version 1.9.2: A. Milani, S. Chesley and D. Vokrouhlichy
c Version 2.0.4: A. Milani, March 2000
C patch 2.1.3: modified by Z. Knezevic July 5, 2000
c version 2.2.6: A. Milani August 2001
c ===================================================================
c PROGRAM FITOBS
c ==================================================================
c Two arc differential corrector, state propagator, 
c     observation predictor
c ===================================================================
      PROGRAM fitobs
      IMPLICIT NONE
c ===== observational data ===========================
c observation numbers: maximum and actual 
      INCLUDE 'parobx.h'
      INTEGER m,mp,mall
c observations both arcs: alpha, delta, time (ET and UT), station code, type
      DOUBLE PRECISION aln(nobx),den(nobx),tau(nobx),tut(nobx)
      INTEGER idsta(nobx),iobs(nobx)
c asteroid identifier, apparent magnitude
      CHARACTER*9 objid(nobx)
      CHARACTER*6 smag(nobx)
c magnitude estimation data
      INCLUDE 'mag.h'
c observation rms, magnitude rms,weight
      DOUBLE PRECISION rmsd(nobx),rmsa(nobx),rmsmag(nobx)
c selection flags; number of observations not discarded
      INTEGER sel(nobx),iob0,iobp,iobm
c weight file names
      CHARACTER*60 rwofi0,rwofip,rwotwo,rwofil
c successful input flags, for observations
      LOGICAL obs0,obsp,obstwo
c radar data present?
      INCLUDE 'radar.h'
c =====state variables: first arc, second arc, joint, current========
c equinoctal orbital elements, epoch times
      DOUBLE PRECISION eq0(6),eqp(6),eq(6)
      DOUBLE PRECISION t0,tp,tm
c flag to select prelim orbit metod
      INTEGER imeth
c asteroid names (18 CHARACTERs)
      CHARACTER*18 astna0,astnap
      CHARACTER*25 astnaj
c  successful input flags, for elements
      LOGICAL ini0,inip,inide,initwo,iniboth
c ===== differential corrections output ==============
c normal, covariance matrices, existence flags for covariance
      DOUBLE PRECISION c0(6,6),cp(6,6),c(6,6)
      DOUBLE PRECISION g0(6,6),gp(6,6),g(6,6)
      LOGICAL cov0,covp,covtwo
c norms of correction, of residuals
      DOUBLE PRECISION delno0,delnop,delnor
      DOUBLE PRECISION csino0,csinop,csinor
c temporary output
      DOUBLE PRECISION eqc(6),gmag
c success flag (for differential corrections)
      LOGICAL succ
c ============ propagation ============================
c target time
      DOUBLE PRECISION tr
c generate ephemerides 
      DOUBLE PRECISION tf,step,interv
      CHARACTER*(8) ans
      INTEGER numsav,ityp
      CHARACTER*3 reqtyp
c control for MTP analysis
      INTEGER iclan
c batch mode for fclan
      LOGICAL batchcl
      DOUBLE PRECISION tlim
c ======== proposed identifications =======
c equinoctal orbital elements (guess from identif), epoch time
      DOUBLE PRECISION eqide(6),tide,tid
c values for joint orbits 
      DOUBLE PRECISION enm
c  number of revolutions for identifications
      INTEGER ng
c identification success flag, proposed elements
      LOGICAL ff
      DOUBLE PRECISION eqf(6)
c ===== predicted observations ===========
c angles, time (TDT, UTC)
      DOUBLE PRECISION  tut1,tut2,t1,t2,dt     
c magnitudes: absolute for all orbits, opp. effect, apparent
      DOUBLE PRECISION h0,hp,hm,hide,gma0,gmap,gmam
c station code, obs. type
      INTEGER ids,iob1
c confidence boundary, line of max variation
      INTEGER  im
c apparent motion
      INCLUDE 'phase.h'
c multiple data for confidence boundary
      INCLUDE 'npoint.h'
c ======== multiple solutions ==========================
      DOUBLE PRECISION sigma
      INCLUDE 'parmul.h'
      INTEGER imult,ifff,iff,imi,nmult,iarm
c first, last, reference solution, interval
      INTEGER imip,imim,imi0,m1,m2
c multiple output arrays
      DOUBLE PRECISION eqm(6,mulx),csinom(mulx),delnom(mulx)
      DOUBLE PRECISION gm(6,6,mulx),cm(6,6,mulx),tc,hmu(mulx)
c ===========close approach analysis======================
c propagation times, close appr. analysys time
      DOUBLE PRECISION trmult,tcmult,tmcla
c input multiple solution catalog
      CHARACTER*160 catname
c minimum v_infty with respect to Earth (circular approx.)
      DOUBLE PRECISION vel_inf
c max value of sigma change allowed for newton step
      DOUBLE PRECISION siglim
c ======== output elements =====================
      CHARACTER*80 elefi0,elefip,eletwo
      INTEGER iunel0,iunelp,iunelt
c ======== output moid =====================
      DOUBLE PRECISION moid, dnp, dnm
      INTEGER iconv
c ========= calendar to julian =================
      DOUBLE PRECISION jd,sec
      INTEGER ihr,imin,iy,imo,iday
c ========= input control ======================
c available data
      LOGICAL ok
c file names depending upon run identifier
      CHARACTER*80 run
      CHARACTER*80 titnam
      CHARACTER*60 filnam
      INTEGER le,lnam
      CHARACTER*6 progna
c logical units
      INTEGER iun20,iun8
c function for time of difcor
      double precision meanti
c ======== controls and flags ===============
c main menus, test derivatives, choice of arc, 
      INTEGER ifun,ifobs,iele,iarc,igue,ider2,iprop,icov,ipred,icop
     + ,iprob
c iope=1 experimental; iope=0 distribution version
      INTEGER iope
      PARAMETER (iope=0)
c characters for menu 
      CHARACTER*20 menunam
      CHARACTER*70 s3,s4,s5,s6,s7,s8,s9,s10
c short circuit 
      LOGICAL init,init2
c asteroids with mass
      LOGICAL found
      INTEGER nfound,lench
c ======== loop indexes =====================
c  ii=subset of 1,6 etc,j=1,6,i=1,m
      INTEGER j,ii,i
c ====== INCLUDE files: trigonometric constants ==================
      INCLUDE 'trig.h'
c ======== constant of gravitation ==============
      INCLUDE 'sunmass.h'
c ===========units for err,pro,clo files ========
      INCLUDE 'proout.h'
      INCLUDE 'verbosity.h'
c =====================================================================
c Run name
      WRITE(*,*) 'Run name ='
      READ(*,100) run
 100  FORMAT(a)
      IF(run.eq.'')stop 'No run specified.'
c input options
      progna='fitobs'
      CALL finopt(progna,run,astna0,astnap,iun20,iun8)
c check for asteroid masses
      IF(lench(astnap).ne.0)THEN
         CALL selpert2(astna0,astnap,nfound)
      ELSE
         CALL selpert(astna0,found)
      ENDIF
c =====================================================================
c initialisations
c =====================================================================
c verbosity levels for an interactive program
      verb_pro=10
      verb_clo=10
      verb_dif=10
      verb_mul=10
      verb_rej=10
c setting of logical flags: nothing is available at the beginning
      obs0=.false.
      obsp=.false.
      obstwo=.false.
      ini0=.false.
      inip=.false.
      inide=.false.
      initwo=.false.
      cov0=.false.
      covp=.false.
      covtwo=.false.
c multiple solutions not yet computed
      imip=0
      imim=0
c setting to zero of counters for data not yet available
      m=0
      mp=0
      mall=0
      t0=0.
      tp=0.
      tm=0.
      tide=0.
c magnitude data: default values
      gma0=0.15d0
      gmap=0.15d0
      gmam=0.15d0
      h0=-1.d9
      hp=-1.d9
      hm=-1.d9
c unit numbers for orbital elements;zero if not opened
      iunel0=0
      iunelp=0
      iunelt=0
c ================SHORT CIRCUIT ==============================
      init=.true.
      init2=.false.
c ================MAIN MENU===================================
c Choice of function
 50   continue
      IF(init2)THEN
         ifun=2
         GOTO 60
      ENDIF
      IF(init)THEN
         ifun=1
         GOTO 60
      ENDIF
      menunam='mainmenu'
      CALL menu(ifun,menunam,10,'What would you like?=',
     +   'input of observational data=',
     +   'acquire orbital elements=',
     +   'differential corrections=',
     +   'first guess for identification=',
     +   'state propagation=',
     +   'predictions of observations=',
     +   'multiple solutions=',
     +   'close approach analysis=',
     +   'status=',
     +   'date conversion=')
 60   IF(ifun.eq.0)THEN
c ==========TERMINATE CLEANLY=========================
         CALL filclo(iun20,' ')
         CALL filclo(iun8,' ')
c close close approach file
         IF(numcla.gt.0)THEN
            CALL filclo(iuncla,' ')
         ELSE
            CALL filclo(iuncla,'DELETE')
         ENDIF
c close propagator parameters file
         CALL filclo(ipirip,'  ')
c close error file
         IF(numerr.gt.0)THEN
            CALL filclo(ierrou,' ')
         ELSE
            CALL filclo(ierrou,'DELETE')
         ENDIF
         STOP
      ELSEIF(ifun.eq.1)THEN
c ================MENU 1: INPUT OBS============================
         IF(init)THEN
            ifobs=3
            GOTO 61
         ENDIF
         WRITE(*,*)' INPUT OF OBSERVATIONAL DATA'
 51      menunam='inputobs'
         CALL menu(ifobs,menunam,5,' which data to input?=',
     +      'first arc=','second arc=','both=',
     +      'undiscard outliers, arc 1=',
     +      'undiscard outliers, arc 2=',
     +      s6,s7,s8,s9,s10)
 61      IF(ifobs.eq.0) GOTO 50

c =====================================================================
c input data, according to request
         IF(ifobs.eq.1.or.ifobs.eq.3)THEN
c ===================================================================== 
c Arc 1 input
            CALL tee(iun20,' INPUT OF OBSERVATIONAL DATA, ARC 1=')
            CALL finobs(progna,1,astna0,objid,obs0,m,iobs,tau,aln,den,
     +         tut,idsta,sel,rmsa,rmsd,rmsmag,smag,nobx,rwofi0,iun20)
c compose elements  file full name
            IF(obs0.and.iunel0.eq.0)THEN
               elefi0=astna0//'.fel'
               CALL rmsp(elefi0,le)
               CALL filopn(iunel0,elefi0(1:le),'unknown')
c output header 
               CALL wromlh (iunel0,'ECLM','J2000')
            ENDIF
         ENDIF
         IF(ifobs.eq.2.or.ifobs.eq.3)THEN
c =====================================================================
c Arc 2 input
            CALL tee(iun20,' INPUT OF OBSERVATIONAL DATA, ARC 2=')
            CALL finobs(progna,2,astnap,objid(m+1),obsp,mp,iobs(m+1),
     +        tau(m+1),aln(m+1),den(m+1),tut(m+1),idsta(m+1),sel(m+1),
     +        rmsa(m+1),rmsd(m+1),rmsmag(m+1),smag(m+1),nobx-m,rwofip,
     +        iun20)
c compose elements  file full name
            IF(obsp.and.iunelp.eq.0)THEN
               elefip=astnap//'.fel'
               CALL rmsp(elefip,le)
               CALL filopn(iunelp,elefip(1:le),'unknown')
c output header 
               CALL wromlh (iunelp,'ECLM','J2000')
            ENDIF
         ENDIF
c =====================================================================
c reintroduce outliers
c Arc 1
         IF(ifobs.eq.4)THEN
            DO i= 1,m
              sel(i)=1
            ENDDO
c Arc 2
         ELSEIF(ifobs.eq.5)THEN
            DO i= 1,mp
              sel(i+m)=1
            ENDDO
         ENDIF
c =====================================================================
c  end input observational data
c 
c summary of input status
         obstwo=obs0.and.obsp
         mall=m+mp
         IF(mp.eq.0)tau(m+1)=0.d0
         iob0=0
         iobp=0
         iobm=0
c find if there are radar data
         CALL radar_ob(iobs,mall)
c weights and elements files for identification
         IF(obstwo.and.iunelt.eq.0)THEN
            CALL titast(3,astna0,astnap,titnam,rwofil,le)
            CALL rmsp(rwofil,le)
            rwotwo=rwofil(1:le)//'.rwo'
c compose elements file full name
            eletwo=rwofil(1:le)//'.fel'
            CALL rmsp(eletwo,le)
            CALL filopn(iunelt,eletwo(1:le),'unknown')
c output header
            CALL wromlh (iunelt,'ECLM','J2000')
         ENDIF
c =====================================================================
      ELSEIF(ifun.eq.2)THEN
c =====================================================================
c input orbital elements
         IF(init2)THEN
            iele=3
            GOTO 62
         ENDIF
         WRITE(*,*)' INPUT OF ORBITAL ELEMENTS'
c ================MENU 2: INPUT ELEMENTS======================
 52      menunam='inputele'
         CALL menu(iele,menunam,8,
     +      ' Which orbital elements to input/compute?=',
     +      ' input arc 1=',' input arc 2=',
     +      ' input both arcs=',
     +      ' compute arc 1 by Gauss/Vaisala method=',
     +      ' compute arc 2 by Gauss/Vaisala method=',
     +      ' compute both arcs by Gauss/Vaisala=',
     +      ' give to arc 2 ele of arc 1=',
     +      ' input first guess from identif=',
     +      s9,s10)
 62      IF(iele.eq.0)GOTO 50
         IF(iele.eq.1.or.iele.eq.3)THEN
c ====================================================================
c initial conditions for arc 1
            CALL tee(iun20,' INPUT OF ORBITAL ELEMENTS, ARC 1=')
            CALL finele(progna,1,obs0,astna0,t0,eq0,h0,gma0,ini0,
     +           cov0,c0,g0,iun20)
            IF(.not.ini0)THEN
               IF(.not.init2)GOTO 52
            ENDIF
         ENDIF
         IF(iele.eq.2.or.iele.eq.3)THEN
c ===================================================================
c initial conditions for arc 2
            CALL tee(iun20,' INPUT OF ORBITAL ELEMENTS, ARC 2=')
            CALL finele(progna,2,obsp,astnap,tp,eqp,hp,gmap,inip,
     +              covp,cp,gp,iun20)
            IF(.not.inip)THEN
               IF(.not.init2)GOTO 52
            ENDIF
         ENDIF
         IF(iele.eq.4.or.iele.eq.6)THEN
c ===================================================================
c use Gauss/Vaisala method for preliminary orbit, arc 1
            IF(.not.obs0)THEN
               WRITE(*,*)'missing observations for arc 1'
               GOTO 52
            ENDIF
            menunam='prelimet'
            CALL menu(imeth,menunam,3,' Which method to use?=',
     +            ' Automatic=',
     +            ' Gauss=',
     +            ' Vaisala=',
     +            s4,s5,s6,s7,s8,s9,s10)
            IF(imeth.eq.1)THEN
              CALL tee(iun20,' AUTO SELECT METHOD, ARC 1=')
            ELSEIF(imeth.eq.2)THEN
              CALL tee(iun20,' GAUSS METHOD, ARC 1=')
            ELSEIF(imeth.eq.3)THEN
              CALL tee(iun20,' VAISALA METHOD, ARC 1=')
            ENDIF
            CALL fgauss(iun20,iunel0,astna0,ini0,cov0,h0,gma0,
     +           rwofi0,tau,tut,aln,den,iobs,idsta,
     +           rmsa,rmsd,objid,smag,rmsmag,sel,m,imeth,eq0,t0)
         ENDIF
         IF(iele.eq.5.or.iele.eq.6)THEN
c =====================================================================
c use Gauss/Vaisala method for preliminary orbit, arc 2
            IF(.not.obsp)THEN
               WRITE(*,*)'missing observations for arc 2'
               GOTO 52
            ENDIF
            menunam='inputmeth'
            CALL menu(imeth,menunam,3,' Which method to use?=',
     +            ' Automatic=',
     +            ' Gauss=',
     +            ' Vaisala=',
     +            s4,s5,s6,s7,s8,s9,s10) 
            IF(imeth.eq.1)THEN
              CALL tee(iun20,' AUTO SELECT METHOD, ARC 2=')
            ELSEIF(imeth.eq.2)THEN
              CALL tee(iun20,' GAUSS METHOD, ARC 2=')
            ELSEIF(imeth.eq.3)THEN
              CALL tee(iun20,' VAISALA METHOD, ARC 2=')
            ENDIF
            CALL fgauss(iun20,iunelp,astnap,inip,covp,hp,gmap,rwofip,
     +        tau(m+1),tut(m+1),aln(m+1),den(m+1),iobs(m+1),idsta(m+1),
     +           rmsa(m+1),rmsd(m+1),objid(m+1),smag(m+1),rmsmag(m+1),
     +           sel(m+1),mp,imeth,eqp,tp)
         ENDIF
         IF(iele.eq.7)THEN
c =====================================================================
c copy elements of arc 1 into elements of arc 2
            IF(ini0)THEN
               tp=t0
               CALL vcopy(6,eq0,eqp)
               hp=h0
               gmap=gma0
               inip=.true.
c covariance is not copied
               covp=.false.
            ELSE
               WRITE(*,*)' initial conditions for arc 1 not available'
            ENDIF
         ENDIF
         IF(iele.eq.8)THEN
c =====================================================================
c  First guess from identif
            IF(iope.eq.0)THEN
               WRITE(*,*)' THIS FUNCTION IS NOT READY'
               GOTO 50
            ENDIF
            CALL tee(iun20,' INPUT PROPOSED IDENTIFICATION=')
c           CALL iniide(astna0,astnap,eqide,tide,inide,iun20)
            IF(.not.inide)THEN
               GOTO 52
            ENDIF
c  magnitude is guessed at average (not very good...)
            hide=(h0+hp)/2.d0
         ENDIF
c initialisation of Yarkovski after acquiring elements 
c (the physical model of the first asteroid is asumed)
         IF(ini0) CALL yarkinit(astna0,'EQU',eq0)
c =====================================================================
      ELSEIF(ifun.eq.3)THEN
         CALL tee(iun20,' DIFFERENTIAL CORRECTIONS=')
c =================MENU 3: DIFFERENTIAL CORRECTIONS====================
 53      menunam='diffcorr'
         CALL menu(iarc,menunam,4,' which arc to correct?=',
     +      'first arc=','second arc=','both=','identification=',
     +      s5,s6,s7,s8,s9,s10)
         IF(iarc.eq.0)GOTO 50
c =====================================================================
         IF(iarc.eq.1.or.iarc.eq.3)THEN
c =====================================================================
c Arc 1 : differential corrections
            CALL tee(iun20,' ARC 1=')
c set slope parameter for current object
            gmagc=gma0            
            CALL fdifco(1,obs0,ini0,ok,cov0,t0,eq0,m,objid,iobs,tau,
     +          tut,idsta,aln,den,rmsa,rmsd,rmsmag,smag,h0,sel,iob0,
     +          rwofi0,iun20,iun8,eqc,g0,c0,csino0,delno0,succ)
c availability of observations, initial condition, JPL and ET-UT data
            IF(.not.ok) GOTO 53
            IF(succ)THEN
c output new elements
               CALL wromlr (iunel0,astna0,eq0,'EQU',t0,g0,.true.,
     +                  c0,.true.,h0,gma0,0.d0)
               CALL nomoid(t0,eq0,moid,iconv,dnp,dnm)
               write(*,199)moid,iconv,dnp,dnm
 199           format('orb.dist.      dist.n+  dist.n-'/ 
     +              f8.5,1x,i4,1x,f8.5,1x,f8.5)
               write(*,*)
               write(iunel0,198)moid,iconv,dnp,dnm
 198           format('!MOID ',f8.5,1x,i4/'!NODES ',f8.5,1x,f8.5)
            ENDIF
         ENDIF
c =====================================================================
         IF(iarc.eq.2.or.iarc.eq.3)THEN
c =====================================================================
c Arc 2 : differential corrections
            CALL tee(iun20,' ARC 2=')
c set slope parameter for current object
            gmagc=gmap
            CALL fdifco(2,obsp,inip,ok,covp,tp,eqp,mp,objid(m+1),
     +        iobs(m+1),tau(m+1),tut(m+1),idsta(m+1),aln(m+1),den(m+1),
     +        rmsa(m+1),rmsd(m+1),rmsmag(m+1),smag(m+1),hp,
     +        sel(m+1),iobp,rwofip,iun20,iun8,eqc,gp,cp,
     +        csinop,delnop,succ)
c availability of observations, initial condition, JPL and ET-UT data
            IF(.not.ok) GOTO 53
            IF(succ)THEN
c output new elements
               CALL wromlr (iunelp,astnap,eqp,'EQU',tp,gp,.true.,
     +                  cp,.true.,hp,gmap,0.d0)
               CALL nomoid(tp,eqp,moid,iconv,dnp,dnm)
               write(*,199)moid,iconv,dnp,dnm
               write(*,*)
               write(iunelp,198)moid,iconv,dnp,dnm
            ENDIF
         ENDIF
c =====================================================================
         IF(iarc.eq.4)THEN
c =====================================================================
c Both arcs to be identified: differential corrections
            CALL tee(iun20,' BOTH ARCS=')
c set slope parameter for current object
            gmagc=gmam
            CALL fdifco(4,obstwo,initwo,ok,covtwo,tm,eq,
     +         mall,objid,iobs,tau,tut,idsta,aln,den,rmsa,rmsd,
     +         rmsmag,smag,hm,sel,iobm,
     +          rwotwo,iun20,iun8,eqc,g,c,csinor,delnor,succ)
c availability of observations, initial condition, JPL and ET-UT data
            IF(.not.ok) GOTO 53
            IF(succ)THEN
c output new elements: problem of the name (used arc 1, but...)
               CALL wromlr (iunelt,astna0,eq,'EQU',tm,g,.true.,
     +                  c,.true.,hm,gmam,0.d0)
               CALL nomoid(tm,eq,moid,iconv,dnp,dnm)
               write(*,199)moid,iconv,dnp,dnm
               write(*,*)
               write(iunelt,198)moid,iconv,dnp,dnm
            ENDIF
c =====================================================================
         ENDIF
c =====================================================================
      ELSEIF(ifun.eq.4)THEN
c =====================================================================
c check availability of JPL ephemerides and ET-UT table
            CALL chetim(tau(1),tau(mall),ok)
            IF(.not.ok) GOTO 50
c ok, go on with arc
         CALL tee(iun20,' FIRST GUESS FOR IDENTIFICATION=')
c ==================MENU 4: FIRST GUESS=======================
 54      menunam='firstgue'
         CALL menu(igue,menunam,6,' which initial guess?=',
     +      'use averages, fit longitudes=',
     +      'use elements of arc 1=',
     +      'use elements of arc 2=',
     +      'copy from ident file=',
     +      'recompute id5=',
     +      'recompute id6=',
     +      s7,s8,s9,s10)
         IF(igue.eq.0)GOTO 50
c =====================================================================
         IF(igue.eq.1)THEN
c =====================================================================
c first guess of a unique solution for both arcs, with averages
c   and fit to longitudes
c  =====================================================================
c check availability of observations and initial condition
            iniboth=ini0.and.inip
            CALL cheobs(obstwo,iniboth,ok)
            IF(.not.ok) GOTO 54
c can be done
            CALL tee(iun20,' GUESS FROM LONGITUDE FIT=')
c =====================================================================
c epoch time in the middle, unless they are too close
            IF(abs(t0-tp).lt.1.d0)THEN
               WRITE(*,*)' initial times ',t0,tp,' too close'
               GOTO 54
            ENDIF
            tm=(t0+tp)/2.d0
            WRITE(iun20,123) tm
            WRITE(*,123) tm
 123        FORMAT(1x,'tm =',f8.1)
            WRITE(iun20,*)
c =====================================================================
c magnitude is mean (not knowing any better)
            hm=(h0+hp)/2.d0
            gmam=(gma0+gmap)/2.d0
c =====================================================================
c use as first guess linear interpolation for h,k,p,q
            DO  ii=2,5
              eq(ii)=(eq0(ii)+eqp(ii))/2.d0
            ENDDO
c Estimate of the mean motion (hence semimajor axis)
c allowing to combine two arcs
            CALL start(eq0,eqp,t0,tp,gms,1,ng,enm,eq(1),eq(6))
            initwo=.true.
            CALL wriequ(iun20,astna0,tm,eq)
            WRITE(*,*)' INITIAL COND. AVAILABLE, NOW USE DIFF. CORR.'
c =====================================================================
         ELSEIF(igue.eq.2)THEN
c =====================================================================
            CALL cheobs(obstwo,ini0,ok)
            IF(.not.ok) GOTO 54
c can be done
            CALL tee(iun20,' USE ELEMENTS OF ARC 1=')
            tm=t0
            CALL vcopy(6,eq0,eq)
            initwo=.true.
            WRITE(*,*)' INITIAL COND. AVAILABLE, NOW USE DIFF. CORR.'
c =====================================================================
c magnitude is first arc
            hm=h0
            gmam=gma0
c =====================================================================
         ELSEIF(igue.eq.3)THEN
c =====================================================================
            CALL cheobs(obstwo,inip,ok)
            IF(.not.ok) GOTO 54
c can be done
            CALL tee(iun20,' USE ELEMENTS OF ARC 2=')
            tm=tp
            CALL vcopy(6,eqp,eq)
            initwo=.true.
            WRITE(*,*)' INITIAL COND. AVAILABLE, NOW USE DIFF. CORR.'
c =====================================================================
c magnitude is second arc
            hm=hp
            gmam=gmap
c =====================================================================
         ELSEIF(igue.eq.4)THEN
c =====================================================================
c use first guess as computed from identif for an arc including
c both sets of observations
c =====================================================================
            IF(iope.eq.0)THEN
               WRITE(*,*)' THIS FUNCTION IS NOT READY'
               GOTO 54
            ENDIF
c check availability of observations and initial condition
            CALL cheobs(obstwo,inide,ok)
            IF(.not.ok) GOTO 54
c =====================================================================
c check availability of JPL ephemerides and ET-UT table
            CALL chetim(tide,tide,ok)
            IF(.not.ok) GOTO 54
c can be done
            CALL tee(iun20,' GUESS FROM IDENTIF=')
c =====================================================================
            tm=tide
            CALL vcopy(6,eqide,eq)
            initwo=.true.
            WRITE(*,*)' INITIAL COND. AVAILABLE, NOW USE DIFF. CORR.'
c =====================================================================
c magnitude is mean (not knowing any better)
            hm=(h0+hp)/2.d0
            gmam=(gma0+gmap)/2.d0
c =====================================================================
         ELSEIF(igue.ge.5)THEN
c =====================================================================
c recompute first guess with the identif algorithm
c =====================================================================
            CALL fident(igue,t0,tp,cov0,covp,eq0,eqp,
     +            g0,c0,gp,cp,tau(m),tau(mall),iun20,eqf,tid,ff)
            IF(ff)THEN
               CALL vcopy(6,eqf,eq)
               tm=tid
               initwo=.true.
               WRITE(*,*)' INITIAL COND. AVAILABLE, NOW USE DIFF. CORR.'
            ELSE
               GOTO 54
            ENDIF
c =====================================================================
c magnitude is mean (not knowing any better)
            hm=(h0+hp)/2.d0
            gmam=(gma0+gmap)/2.d0
c =====================================================================
         ENDIF
c =====================================================================
      ELSEIF(ifun.eq.5)THEN
c =====================================================================
c Computation of asteroid elements at a required time 
         CALL tee(iun20,' PROPAGATION OF ELEMENTS=')
c ================= MENU 5: PROPAGATION =======================
 55      menunam='propagat'
         CALL menu(iprop,menunam,8,'what to propagate?=',
     +      'propagate arc 1, variable time=',
     +      'propagate arc 2, variable time=',
     +      'propagate identified orbit, variable time=',
     +      'propagate arc 1, weighted center of observations=',
     +      'propagate arc 2, weighted center of observations=',
     +      'generate time history of orbital elements, arc 1=',
     +      'generate time history of orbital elements, arc 2=',
     +   'generate time history of orbital elements, identified orbit=',
     +      s9,s10)
         IF(iprop.eq.0) GOTO 50
c =====================================================================
c state vector only? also covariance? (not meaningful for ephemerides)
c =====================================================================
         IF(iprop .le. 5)THEN
            menunam='null'
            CALL menu(icov,menunam,2,'What is required?=',
     +           'orbital elements only=',
     +           'also covariance matrix=',
     +           s3,s4,s5,s6,s7,s8,s9,s10)
            IF(icov.eq.0) GOTO 55
         ELSE
            icov=1
         ENDIF
c =====================================================================
c selection of target epoch
c =====================================================================
         IF(iprop.eq.4)THEN
            tr= meanti(tau,rmsa,rmsd,m)
         ELSEIF(iprop.eq.5)THEN
            tr= meanti(tau(m+1),rmsa(m+1),rmsd(m+1),mp)
         ELSEIF(iprop.ge.6.and. iprop.le.8)THEN
            WRITE(*,*)' Current time is : ',t0,'(MJD).'
            WRITE(*,*)' begin ephemerides from epoch (MJD)?   '
            READ(*,*)tr
            WRITE(*,*)' end ephemerides from epoch (MJD)?   '
            READ(*,*)tf
            WRITE(*,*)' time step in days?'
            READ(*,*)step
            WRITE(*,*) 'Is data correct? (y/n)'
            READ(*,*)ans
            IF(ans(1:1).eq.'n' .or. ans(1:1).eq.'N') GOTO 55 
c           determine number of steps before t0
            IF(tf .lt. t0)THEN
               interv=tf-tr
            ELSEIF(tr .gt. t0)THEN
               interv=0
            ELSE
               interv=t0-tr
            ENDIF
            numsav=interv/step+10
c determine type of elements output
            menunam='null'
            call menu(ityp,menunam,3,'What type of elements?=',
     +           'Keplerian=',
     +           'Equinoctial=',
     +           'Cartesian=',
     +           s4,s5,s6,s7,s8,s9,s10)
            if(ityp.eq.0)then
               goto 55
            elseif(ityp.eq.1)then
               reqtyp='KEP'
            elseif(ityp.eq.2)then
               reqtyp='EQU'
            elseif(ityp.eq.3)then
               reqtyp='CAR'
            endif
c warning: funny result if mp=0; check obsp?
         ELSE
            WRITE(*,*)' propagate to epoch (MJD)?   '
            READ(*,*)tr
         ENDIF
c =====================================================================
c propagation
c =====================================================================
         IF(iprop.eq.1.or.iprop.eq.4)THEN
            CALL tee(iun20,' ARC 1=')
            CALL fstpro(.false.,icov,ini0,cov0,iun20,iun8,ok,
     +         t0,eq0,h0,g0,c0,tr,eq0,h0,g0,c0)
            IF(ok)THEN
               t0=tr
               IF(obs0)THEN
                  CALL wromlr (iunel0,astna0,eq0,'EQU',t0,g0,cov0,
     +                  c0,cov0,h0,gma0,0.d0)
                  CALL nomoid(t0,eq0,moid,iconv,dnp,dnm)
                  write(*,199)moid,iconv,dnp,dnm
                  write(iunel0,198)moid,iconv,dnp,dnm
               ENDIF
            ENDIF
         ELSEIF(iprop.eq.2.or.iprop.eq.5)THEN
            CALL tee(iun20,' ARC 2=')
            CALL fstpro(.false.,icov,inip,covp,iun20,iun8,ok,
     +         tp,eqp,hp,gp,cp,tr,eqp,hp,gp,cp)
            IF(ok)THEN
               tp=tr
               IF(obsp)THEN
                  CALL wromlr (iunelp,astnap,eqp,'EQU',tp,gp,covp,
     +                  cp,covp,hp,gmap,0.d0)
                  CALL nomoid(tp,eqp,moid,iconv,dnp,dnm)
                  write(*,199)moid,iconv,dnp,dnm
                  write(iunelt,198)moid,iconv,dnp,dnm
               ENDIF
            ENDIF
         ELSEIF(iprop.eq.3)THEN
            CALL tee(iun20,' BOTH ARCS=')
            CALL fstpro(.false.,icov,initwo,covtwo,iun20,iun8,ok,
     +         tm,eq,hm,g,c,tr,eq,hm,g,c)
            IF(ok)tm=tr
            IF(ok)THEN
               tm=tr
               IF(obstwo)THEN
c problem with name for identification; used arc 1, but...
                  CALL wromlr (iunelt,astna0,eq,'EQU',tm,g,covtwo,
     +                 c,covtwo,hm,gmam,0.d0)
                  CALL nomoid(tm,eq,moid,iconv,dnp,dnm)
                  write(*,199)moid,iconv,dnp,dnm
                  write(iunelt,198)moid,iconv,dnp,dnm
               ENDIF
            ENDIF
c =====================================================================
c ephemerides (orbital elements) generation
c =====================================================================
         ELSEIF(iprop.eq.6)THEN
            CALL fsteph(astna0,'.',ini0,ok,'EQU',t0,eq0,
     +           tr,tf,step,numsav,.true.,reqtyp,.true.)
         ELSEIF(iprop.eq.7)THEN
            CALL fsteph(astnap,'.',inip,ok,'EQU',tp,eqp,
     +           tr,tf,step,numsav,.true.,reqtyp,.true.)
         ELSEIF(iprop.eq.8)THEN
            astnaj=astna0//'joint'
            CALL rmsp(astnaj,le)
            CALL fsteph(astnaj(1:le),'.',initwo,ok,'EQU',tm,eq,
     +           tr,tf,step,numsav,.true.,reqtyp,.true.)
         ENDIF
c if some data are not available, this cannot be done
         IF(.not.ok)THEN
            WRITE(iun20,*)'    DATA NOT AVAILABLE'
            GOTO 55
         ENDIF
c =====================================================================
      ELSEIF(ifun.eq.6)THEN
c =====================================================================
c prediction of observations
         CALL tee(iun20,' PREDICTION OF OBSERVATIONS=')
c MENU 6: PREDICTIONS 
 56      CALL orb_sel(ini0,inip,initwo,ipred)
c        menunam='prediobs'
         IF(ipred.eq.0) GOTO 50
c =====================================================================
c observations vector only? also covariance?
c ===================================================================
         menunam='predicbd'
         CALL menu(iprob,menunam,5,'What is required?=',
     +      'observations (alpha, delta) only=',
     +      'also covariance matrix=',
     +      'confidence boundary=',
     +      'compare CB with observations=',
     +      'ephemerides (on the sky)=',
     +       s6,s7,s8,s9,s10)
         ok=.true.
         IF(iprob.eq.0)GOTO 56
c ===================================================================
c assign observation time
         IF(iprob.le.4)THEN 
            CALL asstim(iprob,iobs,tau,tut,idsta,m,mall,im,
     +       iob1,t1,tut1,ids)
         ELSE
            CALL seleph(tut1,t1,tut2,t2,dt,ids)
            im=1
c poi passare tdt1,tdt2 a fobopre per epehemc
         ENDIF
c setup title string for graphics output
         CALL titast(ipred,astna0,astnap,titnam,filnam,lnam)
c ===================================================================
c predict
c =====================================================================
         IF(ipred.eq.1)THEN
            CALL tee(iun20,' ARC 1=')
            CALL fobpre(iprob,ini0,cov0,iun20,iun8,ok,
     +           titnam,filnam,t0,eq0,h0,gma0,g0,c0,ids,iob1,t1,
     +           tut1,aln(im),den(im),t2,dt,astna0)
         ELSEIF(ipred.eq.2)THEN
            CALL tee(iun20,' ARC 2=')
            CALL fobpre(iprob,inip,covp,iun20,iun8,ok,
     +           titnam,filnam,tp,eqp,hp,gmap,gp,cp,ids,iob1,t1,
     +           tut1,aln(im),den(im),t2,dt,astnap)
         ELSEIF(ipred.eq.3)THEN
            CALL tee(iun20,' BOTH ARCS=')
            CALL fobpre(iprob,initwo,covtwo,iun20,iun8,ok,
     +           titnam,filnam,tm,eq,hm,gmam,g,c,ids,iob1,t1,
     +           tut1,aln(im),den(im),t2,dt,titnam)
         ENDIF
c if some data are not available, this cannot be done
         IF(.not.ok)THEN
             WRITE(iun20,*)'    DATA NOT AVAILABLE'
             GOTO 56
         ENDIF
c =====================================================================
      ELSEIF(ifun.eq.7)THEN
c =====================================================================
c search for alternate solutions
         CALL tee(iun20,'MULTIPLE SOLUTIONS=')
c MENU 7: MULTIPLE SOLUTIONS
c choice of arc
 58      menunam='multisol'
         CALL menu(iarc,menunam,5,'which orbit?=',
     +      'arc 1=','arc 2=',
     +      'joint computed orbit=',
     +      'use already computed=',
     +      'input from file=',
     +      s6,s7,s8,s9,s10)
         IF(iarc.eq.0) GOTO 50
c =====================================================================
c check availability of initial conditions (also covariance for icov>1)
c and compute multiple solutions
c =====================================================================
c compute multiple solutions
         IF(iarc.eq.1)THEN
            CALL tee(iun20,'ARC 1=')
            tcmult=t0
            CALL mmulti(obs0,ini0,ok,cov0,tcmult,
     +           eq0,g0,c0,csino0,delno0,
     +           m,iobs,tau,idsta,aln,den,rmsa,rmsd,
     +           rmsmag,smag,h0,gma0,sel,
     +           iun20,eqm,gm,cm,hmu,csinom,delnom,
     +           sigma,imult,imim,imip,imi0)
         ELSEIF(iarc.eq.2)THEN
            CALL tee(iun20,'ARC 2=')
            tcmult=tp
            CALL mmulti(obsp,inip,ok,covp,tcmult,
     +           eqp,gp,cp,csinop,delnop,
     +           mp,iobs(m+1),tau(m+1),idsta(m+1),
     +           aln(m+1),den(m+1),rmsa(m+1),rmsd(m+1),
     +           rmsmag(m+1),smag(m+1),hp,gmap,sel(m+1),
     +           iun20,eqm,gm,cm,hmu,csinom,delnom,
     +           sigma,imult,imim,imip,imi0)
         ELSEIF(iarc.eq.3)THEN
            CALL tee(iun20,'BOTH ARCS=')
            tcmult=tm
            CALL mmulti(obstwo,initwo,ok,covtwo,tcmult,
     +           eq,g,c,csinor,delnor,
     +           mall,iobs,tau,idsta,aln,den,rmsa,rmsd,
     +           rmsmag,smag,hm,gmam,sel,
     +           iun20,eqm,gm,cm,hmu,csinom,delnom,
     +           sigma,imult,imim,imip,imi0)
         ELSEIF(iarc.eq.4)THEN
            CALL tee(iun20,'USE THE ALREADY COMPUTED ONES=')
            ok=imip-imim.gt.0
c input from file
         ELSEIF(iarc.eq.5)THEN
            CALL tee(iun20,'INPUT MULTIPLE SOL. FROM FILE=')
            WRITE(*,*)' File name?'
            READ(*,*)catname
            CALL mult_input(catname,eqm,cm,gm,hmu,tcmult,imim,imip,imi0
     +           ,vel_inf,ok)
         ENDIF
         IF(.not.ok)GOTO 58
         nmult=imip-imim+1
         WRITE(*,*)' number of multiple sol. available', nmult
         IF(nmult.le.0)THEN
            CALL tee(iun20,'FAILED MULTIPLE SOLUTIONS')
            GOTO 58
         ENDIF
c ======================================================
c how to use multiple solutions?
c ======================================================
 145     menunam='multiuse'
         CALL menu(ifff,menunam,5,'what to do?=',
     +      'a-e plot of multiple solutions=',
     +      'multiple predicted observations=',
     +      'adopt one of the above solution=',
     +      'propagate multiple solutions=',
     +      'close approach analysys=',
     +      s6,s7,s8,s9,s10)
         IF(ifff.eq.0) GOTO 50
c =================================================================
c setup title string for graphics output
         CALL titast(iarc,astna0,astnap,titnam,filnam,lnam)
c ======================================================
         IF(ifff.eq.1)THEN 
c ======================================================   
c a-e plot of multiple solutions
            IF(iarc.eq.1)THEN
               CALL fmuplo(eqm(1,imim),tcmult,nmult,eq0,titnam,sigma)
            ELSEIF(iarc.eq.2)THEN
               CALL fmuplo(eqm(1,imim),tcmult,nmult,eqp,titnam,sigma)
            ELSEIF(iarc.eq.3)THEN
               CALL fmuplo(eqm(1,imim),tcmult,nmult,eq,titnam,sigma)
            ELSE
               CALL fmuplo(eqm(1,imim),tcmult,nmult,
     +              eqm(1,imi0),titnam,sigma)
            ENDIF
c ======================================================
         ELSEIF(ifff.eq.2)THEN
c ====================================================== 
c multiple predicted observation     
            menunam='null'
            CALL menu(iff,menunam,2,'What is required?=',
     +      'observations (alpha, delta) only=',
     +      'compare with observations=',
     +       s3,s4,s5,s6,s7,s8,s9,s10)
            IF(iff.eq.0) GOTO 145
c =====================================================================
c assign observation time 
            CALL asstim(iff+2,iobs,tau,tut,idsta,m,mall,im,
     +           iob1,t1,tut1,ids)
c =====================================================================
c check availability of JPL ephemerides and ET-UT table
            CALL chetim(t1,t1,ok)
            IF(.not.ok) GOTO 145
c can be done
            CALL tee(iun20,'MULTIPLE PREDICTED OBSERVATIONS=')
c =================================================
c only alpha, delta, magnitude
            IF(iarc.eq.1)THEN
               gmag=gma0
            ELSEIF(iarc.eq.2)THEN
               gmag=gmap
            ELSEIF(iarc.eq.3)THEN
               gmag=gmam
            ENDIF
            CALL fmuobs(tcmult,gmag,iob1,ids,t1,tut1,sigma,eqm,hmu,
     +           aln(im),den(im),iff,imim,imip,imi0,titnam,filnam,iun20)
c ======================================================
         ELSEIF(ifff.eq.3)THEN
c ======================================================
c adopt alternate solution
            WRITE(*,*)' which one to keep? 0=none'
            READ(*,*) imi
            IF(imi.ge.imim.and.imi.le.imip)THEN
               CALL tee(iun20,'ALTERNATE SOLUTION ADOPTED=')
               WRITE(*,*)' NUMBER ',imi
               WRITE(iun20,*)' NUMBER ',imi
               WRITE(iun20,*)
               WRITE(*,144)imi,(eqm(j,imi),j=1,6),csinom(imi)
 144           FORMAT(i3,6f12.8,1p,e13.5,e12.3)
c copy state vector and matrices, norms
               icop=2
c handle case in which the multiple solutions do not come from arc
               IF(iarc.eq.4.or.iarc.eq.5)THEN
 147              WRITE(*,*)' for which arc? 1,2=arcs, 3=joint'
                  READ(*,*)iarm
                  IF(iarm.ge.1.and.iarm.le.3)THEN
                     iarc=iarm
                  ELSE
                     WRITE(*,*)' must be 1,2,3'
                     GOTO 147
                  ENDIF
               ENDIF
               IF(iarc.eq.1)THEN
                  t0=tcmult
                  CALL stacop(icop,eq0,g0,c0,csino0,delno0,
     +               eqm(1,imi),gm(1,1,imi),cm(1,1,imi),
     +               csinom(imi),delnom(imi))
                     h0=hmu(imi)
               ELSEIF(iarc.eq.2)THEN
                  tp=tcmult   
                  CALL stacop(icop,eqp,gp,cp,csinop,delnop,
     +               eqm(1,imi),gm(1,1,imi),cm(1,1,imi),
     +               csinom(imi),delnom(imi))
                     hp=hmu(imi)
               ELSEIF(iarc.eq.3)THEN
                  tm=tcmult   
                  CALL stacop(icop,eq,g,c,csinor,delnor,
     +               eqm (1,imi),gm(1,1,imi),cm(1,1,imi),
     +               csinom(imi),delnom(imi))
                     hm=hmu(imi)
               ELSE
                  WRITE(*,*)' iarm=',iarm,' not understood'
                  GOTO 145
               ENDIF
            ELSE
               CALL tee(iun20,'BACK TO THE ORIGINAL SOLUTION=')
            ENDIF
         ELSEIF(ifff.eq.4)THEN
c =====================================================================
c select propagation time
            CALL tee(iun20,'MULTIPLE PROPAGATION=')
            WRITE(*,*)' Current time is : ',tcmult,'(MJD).'
            WRITE(*,*)' propagate to epoch (MJD)?   '
            READ(*,*)trmult
c check availability of JPL ephemerides
            CALL chetim(tcmult,trmult,ok)
            IF(.not.ok)GOTO 145
c propagation
            IF(iarc.eq.1)THEN
               CALL fmupro(iun20,imim,imip,tcmult,eqm,hmu,
     +           gm,cm,trmult,eqm,gm,cm)
            ELSEIF(iarc.eq.2)THEN
               CALL fmupro(iun20,imim,imip,tcmult,eqm,hmu,
     +           gm,cm,trmult,eqm,gm,cm)
            ELSEIF(iarc.eq.3)THEN
               CALL fmupro(iun20,imim,imip,tcmult,eqm,hmu,
     +           gm,cm,trmult,eqm,gm,cm)
            ENDIF
            tcmult=trmult
c close approach analysis on multiple solutions (as in CLOMON2)
         ELSEIF(ifff.eq.5)THEN
c select interval in index (and sigma) space
 146        WRITE(*,*)' SELECT INTERVAL, between ',imim,' and ',imip
            READ(*,*)m1,m2
            IF(m1.lt.imim.or.m2.gt.imip.or.m1.ge.m2)THEN
               WRITE(*,*)'must be ',imim,' <= ',m1,' < ',m2,' <= ',imip
               WRITE(*,*)' try again'
               GOTO 146
            ENDIF
c select time of propagation
c select final time
            WRITE(*,*)' search for close approaches until time (MJD)?'
            READ(*,*) tmcla
            WRITE(*,*) ' this option not ready yet'
         ENDIF         
c stay inside the multiple orbits case for repeated use of the data
         GOTO 145
c =====================================================================
      ELSEIF(ifun.eq.8)THEN
c =====================================================================
c close approach analysis
c =====================================================================
         CALL tee(iun20,'CLOSE APPROACH ANALYSIS=')
 558     CALL orb_sel(ini0,inip,initwo,iclan)     
c menunam='closapan'

         IF(iclan.eq.0) GOTO 50
         WRITE(*,*)' search for close approaches until time (MJD)?'
         READ(*,*) tlim
c =====================================================================
c check availability of initial conditions (also covariance for icov>1)
c and compute multiple solutions
c =====================================================================
         icov=2
         batchcl=.false.
         IF(iclan.eq.1)THEN
            CALL chereq(icov,ini0,cov0,t0,iun20,iun8,ok)
            IF(ok)THEN
               CALL tee(iun20,'ARC 1=')
               CALL fclan2(batchcl,tlim,iun20,ok,siglim,
     +              eq0,t0,g0,c0,csino0,delno0,h0,gma0,astna0,
     +              obs0,m,objid,iobs,tau,tut,idsta,aln,den,
     +              rmsa,rmsd,sel)
            ENDIF
         ELSEIF(iclan.eq.2)THEN
            CALL chereq(icov,inip,covp,tp,iun20,iun8,ok)
            IF(ok)THEN
               CALL tee(iun20,'ARC 2=')
               CALL fclan2(batchcl,tlim,iun20,ok,siglim,
     +              eqp,tp,gp,cp,csinop,delnop,hp,gmap,astnap,
     +              obsp,mp,objid(m+1),iobs(m+1),tau(m+1),
     +              tut(m+1),idsta(m+1),
     +              aln(m+1),den(m+1),rmsa(m+1),rmsd(m+1),sel(m+1))
            ENDIF
         ELSEIF(iclan.eq.3)THEN
            CALL chereq(icov,initwo,covtwo,tm,iun20,iun8,ok)
            IF(ok)THEN
               CALL tee(iun20,'BOTH ARCS=')
               CALL fclan2(batchcl,tlim,iun20,ok,siglim,
     +              eq,tm,g,c,csinor,delnor,hm,gmam,astna0,
     +              obstwo,mall,objid,iobs,tau,tut,idsta,
     +              aln,den,rmsa,rmsd,sel)
            ENDIF
         ENDIF
c =====================================================================
      ELSEIF(ifun.eq.9)THEN
c =====================================================================
c show all the status flags
         WRITE(*,180)obs0,obsp,ini0,inip,inide,initwo,
     +         cov0,covp,covtwo
 180  FORMAT('   obs0',' obsp ',' ini0 ',' inip ','inide ','initwo',
     +       ' cov0 ',' covp ','covtwo'
     +           /10L6)
c give the epoch times for all the elements
         WRITE(*,181)t0,tp,tm,tide,tau(1),tau(m),tau(m+1),tau(mall)
 181     FORMAT('   t0   ','   tp   ','   tm   ','  tide  ',
     +          '  tin0  ','  tfi0  ','  tinp  ','  tfip  '/
     +          8f8.1)
c observational data available
         IF(obs0)THEN
            WRITE(*,*)' no obs =',m,' of asteroid ',astna0
            WRITE(*,*)'      observations used in fit=',iob0
         ENDIF 
         IF(obsp)THEN
            WRITE(*,*)' no obs =',mp,' of asteroid ',astnap
            WRITE(*,*)'      observations used in fit=',iobp
         ENDIF
         IF(obs0.and.obsp)THEN
            WRITE(*,*)' no obs =',mall,' of both asteroids '
            WRITE(*,*)'      observations used in fit=',iobm
         ENDIF 
c =====================================================================
      ELSEIF(ifun.eq.10)THEN
c =====================================================================
c calendar to MJD conversion
         WRITE(*,*)'calendar date, year, month, day, hour?'
         READ(*,*)iy,imo,iday,ihr
         sec=0.d0
         imin=0
         CALL julian(iy,imo,iday,ihr,imin,sec,jd)
         WRITE(*,777)jd-2400000.5d0
 777     FORMAT('MJD=',f13.6)
c =====================================================================
      ELSEIF(ifun.eq.11)THEN
c =====================================================================
         IF(iope.eq.0)THEN
            WRITE(*,*)' THIS FUNCTION IS NOT READY'
            GOTO 50
         ENDIF
c Check first and possibly second derivatives at the starting points
         WRITE(*,*)' test derivatives of order (1,2)?'
         READ(*,*) ider2
c        CALL twotes(m,t0,tau,idsta,eq0,ider2)
c =====================================================================
      ELSE 
c =====================================================================
c non existing option   
         WRITE(*,*)' This I cannot do'
      ENDIF
      IF(init)THEN
         init=.false.
         init2=.true.
      ELSEIF(init2)THEN
         init2=.false.
      ENDIF
      GOTO 50
      END





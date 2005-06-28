! Copyright (C) 1997,2003 ORBFIT consortium 
! version 3.1: A. Milani December 2003
! ===================================================================   
! PROGRAM FITOBS                                                        
! ==================================================================    
! Two arc differential corrector, state propagator,                     
!     observation predictor                                             
! ===================================================================   
PROGRAM fitobs 
  USE fund_const
  USE astrometric_observations
  USE orbit_elements
  USE output_control
  USE least_squares
  USE multiple_sol
  USE yark_pert
  USE util_suit
  USE force_model, ONLY: radar, norad_obs ! to report on radar data
  USE two_states
  USE attributable
  USE virtual_impactor
  USE tp_trace, ONLY: wri_tppoint,dtpde
  IMPLICIT NONE 
! ===== input observations ======================
  INTEGER ifobs ! menu index
! ===== input orbital elements===================
  INTEGER imeth,iele ! to select prelim orbit metod, which arc
! ===== differential corrections ==============                  
  LOGICAL succ ! success flag                             
  double precision meanti ! function for time of difcor 
  INTEGER idif, itsav ! menu index, save itmax when changed
  DOUBLE PRECISION peq(6) ! unit vector in weak direction
! ======== proposed identifications =======                             
! equinoctal orbital elements (guess from identif), epoch time
  TYPE(orbit_elem) :: elide          
  DOUBLE PRECISION enm ! mean motion for joint orbits  
  INTEGER ng, igue !  number of revolutions, menu index 
  LOGICAL ff  ! identification routine success flag
  LOGICAL inide, iniboth !  successful input flags, for identifications
! ============ propagation ============================
  DOUBLE PRECISION tr ! target time  
  DOUBLE PRECISION tf,step,interv ! generate ephemerides
  INTEGER numsav, iprop,iyes ! ephem records, menu index, yes/no interrogation
  CHARACTER*(8) ans 
! ===== predicted observations ===========                 
  DOUBLE PRECISION  tut1,tut2,t1,t2,dt ! time (TDT, UTC) 
  CHARACTER*(1) type1  ! obs. type
  INTEGER ids, iprob ! station code, menu index                
  INTEGER  im  ! confidence boundary, line of max variation 
! ======== multiple solutions ==========================                
  DOUBLE PRECISION sigma 
  INTEGER imult,ifff,iff,iffat,imi,nmult,iarm,marc,iscal
  INTEGER m1,m2 ! interval  
  CHARACTER*160 catname ! input multiple solution catalog
! propagation times, close appr. analysys time                         
  DOUBLE PRECISION trmult,tcmult,tmcla 
! minimum v_infty with respect to Earth (circular approx.)              
  DOUBLE PRECISION vel_inf 
! max value of sigma change allowed for newton step                     
  DOUBLE PRECISION siglim 
  INTEGER nvi, iunvi
  TYPE(orbit_elem) elop
  INTEGER, PARAMETER :: nj=20
  DOUBLE PRECISION :: del(2,nj)
! ===========close approach analysis======================
  INTEGER iclan 
  DOUBLE PRECISION tlim 
! ===============to change coord=================
  DOUBLE PRECISION, DIMENSION(6,6) :: dee !partials
  CHARACTER*3 cooy ! target type
  INTEGER fail_flag, icoord ! control coo_cha, menu index
! ===============atttributables ===================
  TYPE(attrib) attr0,attrp,attr,attrc !attributables
  DOUBLE PRECISION trou ! rounded time (at the integer MJD)
  DOUBLE PRECISION, PARAMETER :: sphx=2.d0 ! max arc span in degrees
  INTEGER iatt
  LOGICAL error ! only 1 obs and so on
  DOUBLE PRECISION r, rdot ! to complete an attributable to an ATT elem.
! ======== output moid =====================                            
  DOUBLE PRECISION moid0, dnp0, dnm0 
! ========= calendar to julian =================                        
  DOUBLE PRECISION jd,sec 
  INTEGER ihr,imin,iy,imo,iday 
! ========= input control ======================
  LOGICAL ok ! available data 
! file names depending upon run identifier                              
  CHARACTER*80 run 
  CHARACTER*80 titnam 
  CHARACTER*100 filnam,dummyfile 
  INTEGER le,lnam 
  CHARACTER*6 progna 
! logical units                                                         
!  INTEGER iunout,iuncov
! ======== controls and flags =============== 
  LOGICAL batch ! batch control
! main menus, choice of arc,covariance required?,copy to/from,test deriv.
  INTEGER ifun,iarc,icov,icop,ider2
  CHARACTER*20 menunam ! characters for menu
  INTEGER, PARAMETER:: iope=0 ! iope=1 experimental; iope=0 distribution vers.
! short circuit to force input of obs and orbit
  LOGICAL init,init2 
! asteroids with mass                                                   
  LOGICAL found 
  INTEGER nfound,lench 
! ======== loop indexes =====================                           
!  ii=subset of 1,6 etc,j=1,6,i=1,m                                     
  INTEGER j,ii,i 
! ===================================================================== 
! Run name                                                              
  WRITE(*,*) 'Run name =' 
  READ(*,100) run 
100 FORMAT(a) 
  IF(run.eq.'')stop 'No run specified.' 
! input options                                                         
  progna='fitobs' 
  CALL finopt(progna,run,astna0,astnap,error_model) 
  CALL errmod_set(error_model)
! check for asteroid masses                                             
  IF(lench(astnap).ne.0)THEN 
     CALL selpert2(astna0,astnap,nfound) 
  ELSE 
     CALL selpert(astna0,found) 
  ENDIF
! ===================================================================== 
! initialisations   
  batch=.false. ! batch control
! asteroid name for messages
  name_obj=astna0
  CALL rmsp(name_obj,lobjnam)  
! ===================================================================== 
! verbosity levels for an interactive program                           
  verb_pro=10 
  verb_clo=10 
  verb_dif=20 
  verb_mul=10 
  verb_rej=20
  verb_io=10
  verb_moid=20 
! setting of logical flags: nothing is available at the beginning
  CALL set_state_def
  inide=.false. 
  elide=undefined_orbit_elem
  elide%t=0.d0
! intiialization of times, even for non defined elements
  el0=undefined_orbit_elem
  el0%t=0.d0
  elp=undefined_orbit_elem
  elp%t=0.d0    
! multiple solutions not yet computed                                   
  imip=0 
  imim=0 
! unit numbers for orbital elements;zero if not opened                  
  iunel0=0 
  iunelp=0 
  iunelt=0 
! ================SHORT CIRCUIT ==============================          
  init=.true. 
  init2=.false. 
! ================MAIN MENU===================================          
! Choice of function                                                    
50 CONTINUE 
  IF(init2)THEN 
     ifun=2 
     GOTO 60 
  ENDIF
  IF(init)THEN 
     ifun=1 
     GOTO 60 
  ENDIF
  menunam='mainmenu' 
  CALL menu(ifun,menunam,11,'What would you like?=',                &
     &   'input of observational data=',                                &
     &   'acquire orbital elements=',                                   &
     &   'differential corrections=',                                   &
     &   'first guess for identification=',                             &
     &   'state propagation=',                                          &
     &   'predictions of observations=',                                &
     &   'multiple solutions=',                                         &
     &   'coordinate change=',                                          &
     &   'attributables=',                                              &
     &   'status=',                                                     &
     &   'date conversion=')                                            
60 IF(ifun.eq.0)THEN 
! ==========TERMINATE CLEANLY=========================                  
     CALL filclo(iun_log,' ') 
     CALL filclo(iun_covar,' ') 
! close close approach file                                             
     IF(numcla.gt.0)THEN 
        CALL filclo(iuncla,' ') 
     ELSE 
        CALL filclo(iuncla,'DELETE') 
     ENDIF
! close propagator parameters file                                      
     CALL filclo(ipirip,'  ') 
! close error file                                                      
     IF(numerr.gt.0)THEN 
        CALL filclo(ierrou,' ') 
     ELSE 
        CALL filclo(ierrou,'DELETE') 
     ENDIF
     STOP 
  ELSEIF(ifun.eq.1)THEN 
! ================MENU 1: INPUT OBS============================         
     IF(init)THEN 
        ifobs=3 
        GOTO 61 
     ENDIF
     WRITE(*,*)' INPUT OF OBSERVATIONAL DATA' 
51   menunam='inputobs' 
     CALL menu(ifobs,menunam,5,' which data to input?=',            &
     &      'first arc=','second arc=','both=',                         &
     &      'undiscard outliers, arc 1=',                               &
     &      'undiscard outliers, arc 2=')
61   IF(ifobs.eq.0) GOTO 50 
! ===================================================================== 
! input data, according to request                                      
     IF(ifobs.eq.1.or.ifobs.eq.3)THEN 
! ===================================================================== 
! Arc 1 input                                                           
        CALL tee(iun_log,' INPUT OF OBSERVATIONAL DATA, ARC 1=') 
        CALL finobs(progna,1,astna0,obs0,nobx,m,obs,obsw, &
      &          rwofi0,error_model)
        IF(obs0)THEN 
! compose elements  file full name  
           IF(iunel0.eq.0)THEN
              elefi0=astna0//'.fel' 
              CALL rmsp(elefi0,le) 
              CALL filopn(iunel0,elefi0(1:le),'unknown') 
! output header                                                         
              CALL wromlh (iunel0,'ECLM','J2000')
           ENDIF 
        ELSE
           m=0
        ENDIF
     ENDIF
     IF(ifobs.eq.2.or.ifobs.eq.3)THEN 
! ===================================================================== 
! Arc 2 input                                                           
        CALL tee(iun_log,' INPUT OF OBSERVATIONAL DATA, ARC 2=') 
        CALL finobs(progna,2,astnap,obsp,nobx-m,mp,obs(m+1:nobx), &
     &   obsw(m+1:nobx),rwofip,error_model) 
        IF(obsp)THEN 
! compose elements  file full name  
           IF(iunelp.eq.0)THEN
              elefip=astnap//'.fel' 
              CALL rmsp(elefip,le) 
              CALL filopn(iunelp,elefip(1:le),'unknown') 
! output header                                                         
              CALL wromlh (iunelp,'ECLM','J2000')
           ENDIF 
        ENDIF
     ENDIF
! ===================================================================== 
! reintroduce outliers                                                  
! Arc 1                                                                 
     IF(ifobs.eq.4)THEN 
        obsw(1:m)%sel_coord=1
! Arc 2                                                                 
     ELSEIF(ifobs.eq.5)THEN 
        obsw(m+1:m+mp)%sel_coord=1
     ENDIF
! ===================================================================== 
!  end input observational data: summary of input status
     obstwo=obs0.and.obsp 
     mall=m+mp 
     IF(mp.eq.0)obs(m+1)%time_tdt=0.d0 
! find if there are radar data                                          
     CALL radar_ob(obs(1:mall)%type,mall) 
! weights and elements files for identification                         
     IF(obstwo.and.iunelt.eq.0)THEN 
        CALL titast(3,astna0,astnap,titnam,rwofil,le) 
        CALL rmsp(rwofil,le) 
        rwotwo=rwofil(1:le)//'.rwo' 
! compose elements file full name                                       
        eletwo=rwofil(1:le)//'.fel' 
        CALL rmsp(eletwo,le) 
        CALL filopn(iunelt,eletwo(1:le),'unknown') 
! output header                                                         
        CALL wromlh (iunelt,'ECLM','J2000') 
     ENDIF
! ===================================================================== 
  ELSEIF(ifun.eq.2)THEN 
! ===================================================================== 
! input orbital elements                                                
     IF(init2)THEN 
        iele=3 
        GOTO 62 
     ENDIF
     WRITE(*,*)' INPUT OF ORBITAL ELEMENTS' 
! ================MENU 2: INPUT ELEMENTS======================          
52   menunam='inputele' 
     CALL menu(iele,menunam,7,                                      &
     &      ' Which orbital elements to input/compute?=',               &
     &      ' input arc 1=',' input arc 2=',                            &
     &      ' input both arcs=',                                        &
     &      ' compute arc 1 by Gauss/Vaisala method=',                  &
     &      ' compute arc 2 by Gauss/Vaisala method=',                  &
     &      ' compute both arcs by Gauss/Vaisala=',                     &
     &      ' give to arc 2 ele of arc 1=')
62   IF(iele.eq.0)GOTO 50 
     IF(iele.eq.1.or.iele.eq.3)THEN 
! ====================================================================  
! initial conditions for arc 1                                          
        CALL tee(iun_log,' INPUT OF ORBITAL ELEMENTS, ARC 1=') 
        CALL finele(progna,1,astna0,el0,ini0,cov0,unc0) 
        IF(.not.ini0)THEN 
           IF(.not.init2)GOTO 52 
        ENDIF
     ENDIF
     IF(iele.eq.2.or.iele.eq.3)THEN 
! ===================================================================   
! initial conditions for arc 2                                          
        CALL tee(iun_log,' INPUT OF ORBITAL ELEMENTS, ARC 2=') 
        CALL finele(progna,2,astnap,elp,inip,covp,uncp)
        IF(.not.inip)THEN 
           IF(.not.init2)GOTO 52 
        ENDIF
     ENDIF
     IF(iele.ge.4.and.iele.le.6)THEN
! ===================================================================   
! use Gauss/Vaisala method for preliminary orbit
        menunam='prelimet' 
        CALL menu(imeth,menunam,3,' Which method to use?=',         &
     &            ' Automatic=',                                        &
     &            ' Gauss=',                                            &
     &            ' Vaisala=')                                
        IF(imeth.eq.1)THEN 
           CALL tee(iun_log,' AUTO SELECT METHOD=') 
        ELSEIF(imeth.eq.2)THEN 
           CALL tee(iun_log,' GAUSS METHOD=') 
        ELSEIF(imeth.eq.3)THEN 
           CALL tee(iun_log,' VAISALA METHOD=') 
        ENDIF
        IF(iele.eq.4.or.iele.eq.6)THEN 
! use Gauss/Vaisala method for preliminary orbit, arc1
           IF(.not.obs0)THEN 
              WRITE(*,*)'missing observations for arc 1' 
              GOTO 52 
           ENDIF
           CALL tee(iun_log,' PRELIM. ORB. ARC 1=') 
           CALL f_gauss(iunel0,astna0,ini0,cov0,         &
     &           rwofi0,obs,obsw,m,error_model,imeth,el0)        
        ENDIF
        IF(iele.eq.5.or.iele.eq.6)THEN 
! use Gauss/Vaisala method for preliminary orbit, arc 2                 
           IF(.not.obsp)THEN 
              WRITE(*,*)'missing observations for arc 2' 
              GOTO 52 
           ENDIF
           CALL tee(iun_log,' PRELIM. ORB. ARC 2=') 
           CALL f_gauss(iunelp,astnap,inip,covp,rwofip,   &
     &        obs(m+1:m+mp),obsw(m+1:m+mp),mp,error_model,imeth,elp)
        ENDIF
     ENDIF
     IF(iele.eq.7)THEN 
! ===================================================================== 
! copy elements of arc 1 into elements of arc 2                         
        IF(ini0)THEN 
           elp=el0 
           inip=.true. 
! covariance is not copied                                              
           covp=.false. 
        ELSE 
           WRITE(*,*)' initial conditions for arc 1 not available' 
        ENDIF
     ENDIF
! initialisation of Yarkovski after acquiring elements                  
! (the physical model of the first asteroid is assumed)                  
       IF(ini0) CALL yarkinit(astna0,el0) 
! ===================================================================== 
    ELSEIF(ifun.eq.3)THEN 
       CALL tee(iun_log,' DIFFERENTIAL CORRECTIONS=') 
! =================MENU 3: DIFFERENTIAL CORRECTIONS==================== 
53     CALL orb_sel2(.false.,iarc)
       IF(iarc.eq.0)GOTO 50
       CALL obs_cop(1,iarc) ! copy observations to obsc, obswc
       IF(iarc.eq.1)THEN
          rwofic=rwofi0 ! ARC 1
          iunelc=iunel0
       ELSEIF(iarc.eq.2)THEN
          rwofic=rwofip ! ARC 2
          iunelc=iunelp
       ELSEIF(iarc.eq.3)THEN
          rwofic=rwotwo ! ARCS IDENTIFIED
          iunelc=iunelt
       ENDIF
! =========== select mode======================
       menunam='difcomod' 
       CALL menu(idif,menunam,5,' select correction and reject mode=',&
     &      'correct all, autoreject=',                                 &
     &      'correct all, no rejections=',                              &
     &      'constrained solution on LOV (no rejections)=',             &
     &      'correct only some elements (no rejections)=',              &
     &      'compute residuals and covariance (no correction)=')
       IF(idif.eq.0) GOTO 50 
       itsav=itmax 
       IF(idif.eq.5)itmax=0 
! auto/manual reject                                                    
       IF(idif.eq.1)THEN 
          autrej=.true. 
       ELSE 
          autrej=.false. 
       ENDIF
! which elements to correct                                             
       IF(idif.ne.4)THEN 
          interactive=0 
       ELSE 
          interactive=1 
       ENDIF
! ====================constrained solution=========================== 
       IF(idif.eq.3)THEN
! check availability of observations and initial condition              
          CALL cheobs(obsflag,inic,ok) 
          IF(.not.ok) GOTO 53
! check availability of JPL ephemerides and ET-UT table                 
          CALL chetim(obsc(1)%time_tdt,obsc(m)%time_tdt,ok) 
          IF(.not.ok) GOTO 53
          IF(iope.eq.1)THEN 
681          WRITE(*,*)' use scaling, 1=yes, 0=no?'
             READ(*,*)iscal
             IF(iscal.eq.0)THEN
                scaling_lov=.false.
             ELSEIF(iscal.eq.1)THEN
                scaling_lov=.true.
             ELSE
                WRITE(*,*)' answer ', iscal, ' not understood'
                GOTO 681
             ENDIF
682          WRITE(*,*)' which LOV, 1=largest eigenv., 2=second'
             READ(*,*)iscal
             IF(iscal.eq.1)THEN
                second_lov=.false.
             ELSEIF(iscal.eq.2)THEN
                second_lov=.true.
             ELSE
                WRITE(*,*)' answer ', iscal, ' not understood'
                GOTO 682
             ENDIF
          ENDIF
          CALL tee(iun_log,' CONSTRAINED DIFFERENTIAL CORRECTIONS=') 
          CALL constr_fit(mc,obsc,obswc,elc,peq,elc,uncc,csinoc,delnoc,rmshc,iobc,succ)
! =========================full solution=============================
       ELSE 
          CALL tee(iun_log,' FULL DIFFERENTIAL CORRECTIONS=') 
          CALL fdiff_cor(batch,1,obsflag,inic,ok,covc,elc,mc,obsc,obswc,iobc,    &
      &         rwofic,elc,uncc,csinoc,delnoc,rmshc,succ)
       ENDIF
! availability of observations, initial condition, JPL and ET-UT data   
       IF(.not.ok.or..not.succ) GOTO 50 
! output new elements  
       CALL write_elems(elc,astnac,'ML',dummyfile,iunelc,uncc)
       CALL nomoid(elc%t,elc,moid0,dnp0,dnm0) 
       write(*,199)moid0,0,dnp0,dnm0 
199    format('orb.dist.      dist.n+  dist.n-'/                &
            &              f8.5,1x,i4,1x,f8.5,1x,f8.5)                         
       write(*,*) 
       write(iunelc,198)moid0,0,dnp0,dnm0 
198    format('!MOID ',f8.5,1x,i4/'!NODES ',f8.5,1x,f8.5) 
! restore state 
       icop=2
       CALL sta_cop(icop,iarc)
       CALL obs_cop(icop,iarc)
       itmax=itsav
! ===================================================================== 
    ELSEIF(ifun.eq.4)THEN 
! ===================================================================== 
! check availability of JPL ephemerides and ET-UT table                 
       CALL chetim(obs(1)%time_tdt,obs(mall)%time_tdt,ok) 
       IF(.not.ok) GOTO 50 
! ok, go on with arc                                                    
       CALL tee(iun_log,' FIRST GUESS FOR IDENTIFICATION=') 
! ==================MENU 4: FIRST GUESS=======================          
54     menunam='firstgue' 
       CALL menu(igue,menunam,4,' which initial guess?=',             &
     &      'use averages, fit longitudes=',                            &
     &      'use elements of arc 1=',                                   &
     &      'use elements of arc 2=',                                   &
     &      'recompute ident=')
       IF(igue.eq.0)GOTO 50 
! ===================================================================== 
       IF(igue.eq.1)THEN 
! ===================================================================== 
! first guess of a unique solution for both arcs, with averages         
!   and fit to longitudes                                               
!  =====================================================================
! check availability of observations and initial condition              
          iniboth=ini0.and.inip 
          CALL cheobs(obstwo,iniboth,ok) 
          IF(.not.ok) GOTO 54 
! can be done                                                           
          CALL tee(iun_log,' GUESS FROM LONGITUDE FIT=') 
! ===================================================================== 
! epoch time in the middle, unless they are too close                   
          IF(abs(el0%t-elp%t).lt.1.d0)THEN 
             WRITE(*,*)' initial times ',el0%t,elp%t,' too close' 
             GOTO 54 
          ENDIF
! ===================================================================== 
! magnitude is the one of the first (not knowing any better)
          el=el0
          el%t=(el0%t+elp%t)/2.d0 
          WRITE(iun_log,123) el%t 
          WRITE(*,123) el%t 
123       FORMAT(1x,'tm =',f8.1) 
          WRITE(iun_log,*) 
! ===================================================================== 
! use as first guess linear interpolation for h,k,p,q                   
          el%coord(2:5)=(el0%coord(2:5)+elp%coord(2:5))/2.d0 
! Estimate of the mean motion (hence semimajor axis)                    
! allowing to combine two arcs                                          
          CALL start(el0,elp,1,ng,enm,el%coord(1),el%coord(6),ok)
          IF(ok)THEN 
             initwo=.true. 
             CALL wriequ(iun_log,astna0,el%t,el%coord) 
             WRITE(*,*)' INITIAL COND. AVAILABLE, NOW USE DIFF. CORR.'
          ENDIF 
! ===================================================================== 
       ELSEIF(igue.eq.2)THEN 
! ===================================================================== 
          CALL cheobs(obstwo,ini0,ok) 
          IF(.not.ok) GOTO 54 
! can be done                                                           
          CALL tee(iun_log,' USE ELEMENTS OF ARC 1=')  
          el=el0 ! magnitude from first arc 
          initwo=.true. 
          WRITE(*,*)' INITIAL COND. AVAILABLE, NOW USE DIFF. CORR.' 
! ===================================================================== 
       ELSEIF(igue.eq.3)THEN 
! ===================================================================== 
          CALL cheobs(obstwo,inip,ok) 
          IF(.not.ok) GOTO 54 
! can be done                                                           
          CALL tee(iun_log,' USE ELEMENTS OF ARC 2=') 
          el=elp ! magnitude from second arc 
          initwo=.true. 
          WRITE(*,*)' INITIAL COND. AVAILABLE, NOW USE DIFF. CORR.' 
! ===================================================================== 
       ELSEIF(igue.eq.4)THEN 
! ===================================================================== 
! recompute first guess with the identif algorithm                      
! ===================================================================== 
          IF(el0%coo.ne.'EQU'.or.elp%coo.ne.'EQU')THEN
             WRITE(*,*)'not possible unless EQU, given ', el0%coo, elp%coo
             GOTO 54
          ENDIF
          CALL fident(6,cov0,covp,el0,elp,unc0,uncp,elide,ff)
          IF(ff)THEN 
             el=elide 
             initwo=.true. 
             WRITE(*,*)' INITIAL COND. AVAILABLE, NOW USE DIFF. CORR.' 
          ELSE 
             GOTO 54 
          ENDIF
! ===================================================================== 
       ENDIF
! ===================================================================== 
    ELSEIF(ifun.eq.5)THEN 
! ===================================================================== 
! Computation of asteroid elements at a required time                   
       CALL tee(iun_log,' PROPAGATION OF ELEMENTS=') 
! ================= MENU 5: PROPAGATION =======================         
55     menunam='propagat' 
       CALL menu(iprop,menunam,8,'what to propagate?=',               &
     &      'propagate arc 1, variable time=',                          &
     &      'propagate arc 2, variable time=',                          &
     &      'propagate identified orbit, variable time=',               &
     &      'propagate arc 1, weighted center of observations=',        &
     &      'propagate arc 2, weighted center of observations=',        &
     &      'generate time history of orbital elements, arc 1=',        &
     &      'generate time history of orbital elements, arc 2=',        &
     &   'generate time history of orbital elements, identified orbit=')
       IF(iprop.eq.0) GOTO 50 
! ===================================================================== 
! state vector only? also covariance? (not meaningful for ephemerides)  
! ===================================================================== 
       IF(iprop .le. 5)THEN 
          menunam='null' 
          CALL menu(icov,menunam,2,'What is required?=',              &
     &           'orbital elements only=',                              &
     &           'also covariance matrix=')                              
          IF(icov.eq.0) GOTO 55 
       ELSE 
          icov=1 
       ENDIF
! ===================================================================== 
! selection of target epoch                                             
! ===================================================================== 
       IF(iprop.eq.4)THEN 
          tr= meanti(obs(1:m)%time_tdt,obsw(1:m)%rms_coord(1),        &
     &           obsw(1:m)%rms_coord(2),m) 
       ELSEIF(iprop.eq.5)THEN 
          tr= meanti(obs(m+1:m+mp)%time_tdt,obsw(m+1:m+mp)%rms_coord(1),&
     &           obsw(m+1:m+mp)%rms_coord(2),mp) 
       ELSEIF(iprop.ge.6.and. iprop.le.8)THEN 
          WRITE(*,*)' Current time is : ',el0%t,'(MJD).' 
          WRITE(*,*)' begin ephemerides from epoch (MJD)?   ' 
          READ(*,*)tr 
          WRITE(*,*)' end ephemerides from epoch (MJD)?   ' 
          READ(*,*)tf 
          WRITE(*,*)' time step in days?' 
          READ(*,*)step 
          WRITE(*,*) 'Is data correct? (y/n)' 
          READ(*,*)ans 
          IF(ans(1:1).eq.'n' .or. ans(1:1).eq.'N') GOTO 55 
!   determine number of steps before t0                         
          IF(tf .lt. el0%t)THEN 
             interv=tf-tr 
          ELSEIF(tr .gt. el0%t)THEN 
             interv=0 
          ELSE 
             interv=el0%t-tr 
          ENDIF
          numsav=interv/step+10 
! determine type of elements output                                     
          menunam='coord' 
          call menu(icoord,menunam,5,'What type of elements?=',         &
     &           'KEPlerian=', 'EQUinoctial=',                          &
     &           'CARtesian=', 'COMetary=','ATTributable=')              
          IF(icoord.eq.0)THEN 
             goto 55 
          ELSE
             cooy=cootyp(icoord)
          ENDIF
! warning: funny result if mp=0; check obsp?                            
       ELSE 
          WRITE(*,*)' propagate to epoch (MJD)?   ' 
          READ(*,*)tr 
       ENDIF
! ===================================================================== 
! propagation                                                           
! ===================================================================== 
       IF(iprop.eq.1.or.iprop.eq.4)THEN
          iarc=1
          iunelc=iunel0
       ELSEIF(iprop.eq.2.or.iprop.eq.5)THEN 
          iarc=2
          iunelc=iunelp
       ELSEIF(iprop.eq.3)THEN 
          iarc=3
          iunelc=iunelt
       ENDIF
       CALL orb_sel2(.true.,iarc)
       CALL fstpro(.false.,icov,inic,covc,iun_log,iun_covar,ok,       &
     &         elc,uncc,tr,elc,uncc)                               
       IF(ok)THEN 
! problem with name for identification; used arc 1, but...  
          CALL write_elems(elc,astnac,'ML',dummyfile,iunelc,uncc)
          CALL nomoid(elc%t,elc,moid0,dnp0,dnm0) 
          write(*,199)moid0,0,dnp0,dnm0 
          write(iunelc,198)moid0,0,dnp0,dnm0 
! copy into the right state
          CALL sta_cop(2,iarc)
       ENDIF
       GOTO 50
! ===================================================================== 
! ephemerides (orbital elements) generation                             
! ===================================================================== 
       IF(iprop.eq.6)THEN 
          iarc=1
       ELSEIF(iprop.eq.7)THEN 
          iarc=2
       ELSEIF(iprop.eq.8)THEN 
          astnaj=astna0//'joint' 
          iarc=3
       ENDIF
       CALL orb_sel2(.true.,iarc)
       CALL fsteph(astnac,'.',inic,ok,elc,                &
     &           tr,tf,step,numsav,.true.,cooy,.true.) 
! if some data are not available, this cannot be done                   
       IF(.not.ok)THEN 
          WRITE(iun_log,*)'    DATA NOT AVAILABLE' 
          GOTO 55 
       ENDIF
! ===================================================================== 
    ELSEIF(ifun.eq.6)THEN 
! ===================================================================== 
! prediction of observations                                            
       CALL tee(iun_log,' PREDICTION OF OBSERVATIONS=') 
! MENU 6: PREDICTIONS                                                   
56     CALL orb_sel2(.false.,iarc)                                            
       IF(iarc.eq.0) GOTO 50 
! setup title string for graphics output                                
       CALL titast(iarc,astna0,astnap,titnam,filnam,lnam)
! ===================================================================== 
! observations vector only? also covariance?                            
! ==================================================================   
! move to fobpre???
       menunam='predicbd' 
       CALL menu(iprob,menunam,6,'What is required?=',                &
     &      'observations (alpha, delta) only=',                      &
     &      'use simulated observation=',                             &
     &      'also covariance matrix=',                                &
     &      'confidence boundary=',                                   &
     &      'compare CB with observations=',                          &
     &      'ephemerides (on the sky)=')
       ok=.true. 
       IF(iprob.eq.0)GOTO 50 
! ===================================================================   
! assign observation time                                               
556    IF(iprob.le.5)THEN 
          CALL asstim(iprob,obs(1:mall)%type,obs(1:mall)%time_tdt,    &
    &        obs(1:mall)%time_utc,obs(1:mall)%obscod_i,m,mall,im,     &
     &       type1,t1,tut1,ids)                                          
       ELSE 
          CALL seleph(tut1,t1,tut2,t2,dt,ids) 
          im=1 
       ENDIF
! ===================================================================   
! predict                                                               
! ===================================================================== 
       CALL fobpre(iprob,inic,covc,ok,                &
     &           titnam,filnam,elc,uncc,ids,type1,t1,        &
     &           tut1,obs(im)%coord(1),obs(im)%coord(2),t2,dt,astnac)
! if some data are not available, this cannot be done                   
       IF(.not.ok)THEN 
          WRITE(iun_log,*)'    DATA NOT AVAILABLE' 
          GOTO 56 
       ENDIF
       IF(iprob.eq.2)THEN
          WRITE(*,*)'more simulated observations? 0=no'
          READ(*,*)iyes
          IF(iyes.ne.0)GOTO 556
       ENDIF
! ===================================================================== 
    ELSEIF(ifun.eq.7)THEN 
! ===================================================================== 
! search for alternate solutions                                        
       CALL tee(iun_log,'MULTIPLE SOLUTIONS=') 
! MENU 7: MULTIPLE SOLUTIONS                                            
! choice of arc                                                         
58     menunam='multisol' 
       CALL menu(marc,menunam,5,'which orbit?=',                      &
     &      'arc 1=','arc 2=',                                          &
     &      'joint computed orbit=',                                    &
     &      'use already computed=',                                    &
     &      'input from file=')
       IF(marc.eq.0) GOTO 50
       IF(iope.eq.1)THEN 
581    WRITE(*,*)' use scaling, 1=yes, 0=no?'
       READ(*,*)iscal
       IF(iscal.eq.0)THEN
          scaling_lov=.false.
       ELSEIF(iscal.eq.1)THEN
          scaling_lov=.true.
       ELSE
          WRITE(*,*)' answer ', iscal, ' not understood'
          GOTO 581
       ENDIF 
582    WRITE(*,*)' which LOV, 1=largest eigenv., 2=second'
       READ(*,*)iscal
       IF(iscal.eq.1)THEN
          second_lov=.false.
       ELSEIF(iscal.eq.2)THEN
          second_lov=.true.
       ELSE
          WRITE(*,*)' answer ', iscal, ' not understood'
          GOTO 582
       ENDIF 
       ENDIF
! ===================================================================== 
! check availability of initial conditions (also covariance for icov>1) 
! and compute multiple solutions                                        
! ===================================================================== 
! compute multiple solutions                                            
       IF(marc.ge.1.and.marc.le.3)THEN
          CALL orb_sel2(.true.,marc)
          CALL obs_cop(1,marc) ! copy observations to obsc, obswc
          CALL tee(iun_log,'COMPUTE MULTIPLE SOL=') 
          CALL f_multi(batch,obsflag,inic,ok,covc, &
     &           elc,uncc,csinoc,delnoc,mc,obsc,obswc,sigma,imult)         
       ELSEIF(marc.eq.4)THEN 
          CALL tee(iun_log,'USE THE ALREADY COMPUTED ONES=') 
          ok=imip-imim.gt.0 
! input from file                                                       
       ELSEIF(marc.eq.5)THEN 
          CALL tee(iun_log,'INPUT MULTIPLE SOL. FROM FILE=') 
          WRITE(*,*)' File name?' 
          READ(*,*)catname 
          CALL mult_input(catname,ok)
       ENDIF
       IF(.not.ok)GOTO 58 
       nmult=imip-imim+1 
       WRITE(*,*)' number of multiple sol. available', nmult 
       IF(nmult.le.0)THEN 
          CALL tee(iun_log,'FAILED MULTIPLE SOLUTIONS') 
          GOTO 58 
       ENDIF
! ======================================================                
! how to use multiple solutions?                                        
! ======================================================                
145    menunam='multiuse' 
       CALL menu(ifff,menunam,5,'what to do?=',                       &
     &      'plot of multiple solutions=',                          &
     &      'multiple predicted observations=',                         &
     &      'adopt one of the above solution=',                         &
     &      'propagate multiple solutions=',                            &
     &      'close approach analysys=')
       IF(ifff.eq.0) GOTO 50 
! ================================================================= 
       CALL orb_sel2(.false.,iarc) 
! setup title string for graphics output                                
       CALL titast(iarc,astna0,astnap,titnam,filnam,lnam) 
! ======================================================                
       IF(ifff.eq.1)THEN 
! ======================================================                
! a-e plot of multiple solutions                                        
          IF(marc.ge.4) elc=elm(imi0)
          CALL fmuplo(titnam,sigma)
! ======================================================                
       ELSEIF(ifff.eq.2)THEN 
! ======================================================                
! multiple predicted observation                                        
          menunam='null' 
          CALL menu(iff,menunam,2,'What is required?=',               &
     &      'observations (alpha, delta) only=',                        &
     &      'compare with observations=')
          IF(iff.eq.0) GOTO 145 
! ===================================================================== 
! assign observation time 
          IF(iff.eq.2)iffat=5 
          IF(iff.eq.1)iffat=3                                             
          CALL asstim(iff+2,obs%type,obs%time_tdt,obs%time_utc,           &
     &           obs%obscod_i,m,mall,im,type1,t1,tut1,ids)
! ===================================================================== 
! check availability of JPL ephemerides and ET-UT table                 
          CALL chetim(t1,t1,ok) 
          IF(.not.ok) GOTO 145 
! can be done                                                           
          CALL tee(iun_log,'MULTIPLE PREDICTED OBSERVATIONS=') 
! =================================================                     
! only alpha, delta, magnitude                                          
          CALL fmuobs(type1,ids,t1,tut1,sigma,     &
     &      obs(im)%coord(1),obs(im)%coord(2),iff, &
     &            titnam,filnam,iun_log)
! ======================================================                
       ELSEIF(ifff.eq.3)THEN 
! ======================================================                
! adopt alternate solution                                              
          WRITE(*,*)' which one to keep? 0=none' 
          READ(*,*) imi 
          IF(imi.ge.imim.and.imi.le.imip)THEN 
             CALL tee(iun_log,'ALTERNATE SOLUTION ADOPTED=') 
             WRITE(*,*)' NUMBER ',imi 
             WRITE(iun_log,*)' NUMBER ',imi 
             WRITE(iun_log,*) 
             WRITE(*,144)imi,elm(imi)%coord,csinom(imi) 
144          FORMAT(i3,6f12.8,1p,e13.5,e12.3) 
! copy state vector and matrices, norms                                 
             icop=2 
! handle case in which the multiple solutions do not come from arc      
             IF(iarc.eq.4.or.iarc.eq.5)THEN 
147             WRITE(*,*)' for which arc? 1,2=arcs, 3=joint' 
                READ(*,*)iarm 
                IF(iarm.ge.1.and.iarm.le.3)THEN 
                   iarc=iarm 
                ELSE 
                   WRITE(*,*)' must be 1,2,3' 
                   GOTO 147 
                ENDIF
             ENDIF
             IF(iarc.eq.1)THEN 
                CALL stacop(icop,el0,unc0,csino0,delno0,             &
     &               elm(imi),unm(imi),csinom(imi),delnom(imi))
             ELSEIF(iarc.eq.2)THEN 
                CALL stacop(icop,elp,uncp,csinop,delnop,             &
     &               elm(imi),unm(imi),csinom(imi),delnom(imi))
             ELSEIF(iarc.eq.3)THEN 
                CALL stacop(icop,el,unc,csinor,delnor,                &
     &               elm(imi),unm(imi),csinom(imi),delnom(imi)) 
             ELSE 
                WRITE(*,*)' iarm=',iarm,' not understood' 
                GOTO 145 
             ENDIF
          ELSE 
             CALL tee(iun_log,'BACK TO THE ORIGINAL SOLUTION=') 
          ENDIF
       ELSEIF(ifff.eq.4)THEN 
! ===================================================================== 
! select propagation time                                               
          CALL tee(iun_log,'MULTIPLE PROPAGATION=') 
          tcmult=elm(imi0)%t
          WRITE(*,*)' Current time is : ',tcmult,'(MJD).' 
          WRITE(*,*)' propagate to epoch (MJD)?   ' 
          READ(*,*)trmult 
! check availability of JPL ephemerides                                 
          CALL chetim(tcmult,trmult,ok) 
          IF(.not.ok)GOTO 145 
! propagation                                                           
          CALL fmupro(iun_log,trmult)
! close approach analysis on multiple solutions (as in CLOMON2)         
       ELSEIF(ifff.eq.5)THEN 
          CALL tee(iun_log,'CLOSE APPROACH ANALYSIS=')
          filnam=run//'.vis' 
          CALL rmsp(filnam,le) 
          call filopn(iunvi,filnam(1:le),'UNKNOWN')           
! select interval in index (and sigma) space                            
146       WRITE(*,*)' SELECT INTERVAL, between ',imim,' and ',imip 
          READ(*,*)m1,m2 
          IF(m1.lt.imim.or.m2.gt.imip.or.m1.ge.m2)THEN 
             WRITE(*,*)'must be ',imim,' <= ',m1,' < ',m2,' <= ',imip 
             WRITE(*,*)' try again' 
             GOTO 146 
          ENDIF 
          WRITE(iun_log,*)' VA interval between m1=',m1,' m2=',m2
! select final time                                                     
          WRITE(*,*)' search for close approaches until time (MJD)?' 
          READ(*,*) tmcla 
          WRITE(iun_log,*)' search until time',tmcla,' MJD'
          CALL fclomon2(progna,mc,obsc,obswc,m1,m2,tmcla)          
          IF(num_vi.le.0)THEN
             CALL filclo(iunvi,'DELETE')
             GOTO 145
          ENDIF
          DO nvi=1,num_vi
             CALL wri_tppoint(vis(nvi)%tp,iunvi,.true.)
             CALL write_elems(vis(nvi)%ele,astnac,'ML',dummyfile,iunvi,vis(nvi)%unc)
          ENDDO
          CALL filclo(iunvi,' ')
! what to do with VIs
          IF(iope.eq.0) GOTO 145
 345      CONTINUE
          DO nvi=1,num_vi
             WRITE(*,*)nvi,vis(nvi)%tp%tcla, vis(nvi)%tp%b, vis(nvi)%tp%txi, vis(nvi)%tp%tze
          ENDDO
          WRITE(*,*)' which VI to analyse? 0=none'
          READ(*,*) nvi
          IF(nvi.le.0) GOTO 145
          curr_vi=vis(nvi)
          CALL vi_draw(del)
          GOTO 345
       ENDIF
! stay inside the multiple orbits case for repeated use of the data     
       GOTO 145 
! ===================================================================== 
    ELSEIF(ifun.eq.8)THEN 
! ===================================================================== 
! coordinate changes
       CALL tee(iun_log,'COORDINATE CHANGES=') 
         menunam='coord' 
         CALL menu(icoord,menunam,5,'Coordinates?=',          &
     &      'KEPlerian=',                                   &
     &      'EQUinoctal=',                                  &
     &      'CARtesian=',                                   &
     &      'COMetary=',                                    &
     &      'ATTributables=')
         ok=.true. 
         IF(icoord.eq.0)GOTO 50
         cooy=cootyp(icoord)
         menunam='cooarc' 
         CALL menu(iarc,menunam,4,                               &
     &      ' Which orbital elements to convert?=',          &
     &      ' arc 1=',' arc 2=',' both arcs=',               &
     &      ' identification=')
         IF(iarc.eq.1.or.iarc.eq.3)THEN
            IF(ini0)THEN
               IF(cov0)THEN
                  CALL coo_cha(el0,cooy,el0,fail_flag,dee)
                  CALL convertunc(unc0,dee,unc0)
               ELSE
                  CALL coo_cha(el0,cooy,el0,fail_flag)
                  unc0%succ=.false.
               ENDIF
               CALL write_elems(el0,astna0,'ML',dummyfile,iunel0,unc0)
               WRITE(*,*)' elements for arc 1', el0
            ELSE
               WRITE(*,*)' initial conditions not available for ',astna0
            ENDIF
         ENDIF
         IF(iarc.eq.2.or.iarc.eq.3)THEN
            IF(inip)THEN
               IF(covp)THEN
                  CALL coo_cha(elp,cooy,elp,fail_flag,dee)
                  CALL convertunc(uncp,dee,uncp)
               ELSE
                  CALL coo_cha(elp,cooy,elp,fail_flag)
                  uncp%succ=.false.
               ENDIF
               CALL write_elems(elp,astnap,'ML',dummyfile,iunelp,uncp)
               WRITE(*,*)' elements for arc 2', elp
            ELSE
               WRITE(*,*)' initial conditions not available for ',astnap
            ENDIF
         ENDIF
         IF(iarc.eq.4)THEN
            IF(initwo)THEN
               IF(covtwo)THEN
                  CALL coo_cha(el,cooy,el,fail_flag,dee)
                  CALL convertunc(unc,dee,unc)
               ELSE
                  CALL coo_cha(el,cooy,el,fail_flag)
                  unc%succ=.false.
               ENDIF
               CALL write_elems(el,astna0,'ML',dummyfile,iunelt,unc)
               WRITE(*,*)' elements for both arcs', el
            ELSE
               WRITE(*,*)' initial conditions not available for '  &
                    &                ,astna0//'='//astnap
            ENDIF
         ENDIF
! ===================================================================== 
      ELSEIF(ifun.eq.9)THEN 
! ===================================================================== 
! compute attributables
         CALL tee(iun_log,'COMPUTE ATTRIBUTABLES=')
         CALL orb_sel(.false.,iarc)
         IF(iarc.eq.0) GOTO 50 
         CALL obs_cop(1,iarc) ! copy observations to obsc, obswc
         CALL attri_comp(mc,obsc,obswc,attrc,error)
         IF(error) GOTO 50
         trou=nint(attrc%tdtobs)
         IF(attrc%sph*radeg.gt.sphx)THEN
            WRITE(*,*)' arc too wide ', attrc%sph*radeg
         ELSE
            CALL wri_attri(0,0,astnac,attrc,trou)
         ENDIF
         IF(iarc.eq.1)THEN
            attr0=attrc
            elc=el0
         ELSEIF(iarc.eq.2)THEN
            attrp=attrc
            elc=elp
         ELSEIF(iarc.eq.3)THEN
            attr=attrc
            elc=el
         ENDIF
         menunam='dummy' 
         CALL menu(iatt,menunam,3,'What to do with attributable?=', &
 &                    'Assign r, rdot to form ATT elements=',       &
 &                    'Compare with predictions=',                  &
 &                    'Output to file=')
         IF(iatt.eq.0) GOTO 50
         IF(iatt.eq.1)THEN
 599        WRITE(*,*) 'assign r (AU) '
            READ(*,*,ERR=599) r
 598        WRITE(*,*) 'assign rdot (AU/day) '
            READ(*,*,ERR=598) rdot
            CALL attelements(attrc,r,rdot,elc,uncc)
            IF(iarc.eq.1)THEN
               el0=elc
               ini0=.true.
               cov0=.true.
               unc0=uncc
               WRITE(*,*)' ATT elements for arc 1'
               WRITE(*,*) el0
            ELSEIF(iarc.eq.2)THEN
               elp=elc
               inip=.true.
               covp=.true.
               uncp=uncc
               WRITE(*,*)' ATT elements for arc 2'
               WRITE(*,*) elp
            ELSEIF(iarc.eq.3)THEN
               el=elc
               initwo=.true.
               covtwo=.true.
               unc=uncc
               WRITE(*,*)' ATT elements for identification'
               WRITE(*,*) elc
            ENDIF
         ELSE
            WRITE(*,*)iatt, ' option not operational'
         ENDIF
! ===================================================================== 
      ELSEIF(ifun.eq.10)THEN 
! ===================================================================== 
! show all the status flags                                             
         WRITE(*,180)obs0,obsp,ini0,inip,inide,initwo,                  &
     &         cov0,covp,covtwo                                         
  180 FORMAT('   obs0',' obsp ',' ini0 ',' inip ','inide ','initwo',    &
     &       ' cov0 ',' covp ','covtwo'                                 &
     &           /10L6)                                                 
! give the epoch times for all the elements                            
         WRITE(*,181)el0%t,elp%t,el%t,elide%t, &
     & obs(1)%time_tdt,obs(m)%time_tdt,obs(m+1)%time_tdt,obs(mall)%time_tdt
  181    FORMAT('   t0   ','   tp   ','   tm   ','  tide  ',            &
     &          '  tin0  ','  tfi0  ','  tinp  ','  tfip  '/            &
     &          8f8.1)                                                  
! observational data available                                          
         IF(obs0)THEN 
            WRITE(*,*)' no obs =',m,' of asteroid ',astna0 
            IF(cov0)THEN
               WRITE(*,*)'      obs used in fit=',iob0, ' RMS =',csino0
            ELSE
               WRITE(*,*)'      obs used in fit=',iob0
            ENDIF 
         ENDIF 
         IF(obsp)THEN 
            WRITE(*,*)' no obs =',mp,' of asteroid ',astnap 
            IF(covp)THEN
               WRITE(*,*)'      obs used in fit=',iobp, ' RMS =',csinop
            ELSE
               WRITE(*,*)'      observations used in fit=',iobp
            ENDIF 
         ENDIF 
         IF(obs0.and.obsp)THEN 
            WRITE(*,*)' no obs =',mall,' of both asteroids ' 
            IF(covtwo)THEN
               WRITE(*,*)'      obs used in fit=',iobtwo, ' RMS =',csinor
            ELSE
               WRITE(*,*)'      observations used in fit=',iobtwo
            ENDIF
         ENDIF 
         IF(ini0)WRITE(*,*)' coord. orbit of ',astna0,' are ',el0%coo
         IF(inip)WRITE(*,*)' coord. orbit of ',astnap,' are ',elp%coo
         IF(initwo)WRITE(*,*)' coord. orbit of ',astna0//'='//astnap,' are ',el%coo
! ===================================================================== 
      ELSEIF(ifun.eq.11)THEN 
! ===================================================================== 
! calendar to MJD conversion                                            
         WRITE(*,*)'calendar date, year, month, day, hour?' 
         READ(*,*)iy,imo,iday,ihr 
         sec=0.d0 
         imin=0 
         CALL julian(iy,imo,iday,ihr,imin,sec,jd) 
         WRITE(*,777)jd-2400000.5d0 
  777    FORMAT('MJD=',f13.6) 
! ===================================================================== 
      ELSEIF(ifun.eq.13)THEN 
! ===================================================================== 
         IF(iope.eq.0)THEN 
            WRITE(*,*)' THIS FUNCTION IS NOT READY' 
            GOTO 50 
         ENDIF 
! Check first and possibly second derivatives at the starting points    
         WRITE(*,*)' test derivatives of order (1,2)?' 
         READ(*,*) ider2 
!        CALL twotes(m,t0,tau,idsta,eq0,ider2)                          
! ===================================================================== 
      ELSE 
! ===================================================================== 
! non existing option                                                   
         WRITE(*,*)' This I cannot do' 
      ENDIF 
      IF(init)THEN 
         init=.false. 
         init2=.true. 
      ELSEIF(init2)THEN 
         init2=.false. 
      ENDIF 
      GOTO 50 
    END PROGRAM fitobs

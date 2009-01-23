! Copyright (C) 1998 by OrbFit Consortium                               
! Version: December 14, 1998 Mario Carpino                              
! ================================================================      
! RMODEL                                                                
! ================================================================      
! this subroutine reads from file 'propag.def' all the propagator option
!  WARNING: we advise the user against changing the file propag.def;    
!                do it at your risk !....                               
! input options for phisical model! version 1.2 of propagator           
!       ilun=0!  0=no moon  1= yes                                      
!       imerc=1! 0=no mercury 1=yes (recommended =1 for asteroids)      
!       iplut=0! 0=no pluto 1=yes (recommended =0 for asteroids, =1 for 
!                                   transneptunian)                     
!       irel=0!  0=newtonian 1=gen. relativity (recommended =0 for main 
!                                   =1 for Earth crossing)              
!       iast=0!  0=no asteroids with mass n=no. of massive asteroids    
!       filbe='CPV'! name of the asteroid ephemerides file (written by B
!       iclap=1! 0=no close approach control 1=yes (recommended =1)     
!       iaber=1! aberration 0=no 1=yes (recommended =1 always)          
!       istat=1! 0=no topocentric corr 1= yes (recommended =1 always)   
!    iclap.h90                                                               
!       iclap = 1 means for Radau alg. to calculate time pos. and vel.  
!                at closest approach; for multistep alg. the propagator 
!                detects the time of close-appr.                        
!      iclap =0 the subroutine force does not check for possible close  
!                approaches to major planets andor massive asteroids.
!   The minimal distance to define a close-approch is controlled by:    
!              dmea (for the Earth only), dter (for Mercury, Venus and M
!              dmoon (for the moon), dmjup (for giant planets),         
!              dmast (for massive asteroids).     
! in module yark_pert   
!       iyark=1 ! Yarkovski force                                       
!       iyarpt=0! partials of yarkovski are generally not computed      
! for rhs=2
! options applicable only to satellite case
! .ites=2      ! max harmonic degree
! .irad=0      ! radiation pressure 1=spher.sat
! .itide=0     ! tidal perturbation 1=k2 no lag
! .ipla=2      ! 0=2-body 2=Sun+Moon
! ================================================================      
! ========INTERFACE============================                         
SUBROUTINE rmodel(rhs0)
! add a lot of ONLY here 
! it needs only these
  USE least_squares, ONLY : difini, rejini, output_old_rwo
  USE propag_state, ONLY : inipro ! calls inipro
  USE tp_trace, ONLY: tpplane
  USE close_app ! sets former closapl options
  USE yark_pert !
  USE non_grav
  USE perturbations
  USE planet_masses ! 
  USE force_model
  USE force_sat
  USE spher_harm
  USE astrometric_observations, ONLY: radius ! to set default value
  USE fund_const
  USE output_control
  USE reference_systems
  IMPLICIT NONE 
  INTEGER, INTENT(IN) :: rhs0
! ========END INTERFACE========================
! no headers:
! close approach options in force_model.mod
! radar quality flag in force_model.mod
! default asteroid radius in astrometric_observations.mod
  CHARACTER*80 filbe                                               
  LOGICAL fail,fail1,found 
  INTEGER ll,ia 
! controls for bizarre orbits                                           
  DOUBLE PRECISION ecclim, samin,samax,phmin,ahmax,qmax 
! initialization of rhs
  rhs=rhs0
!****************                                                       
!   static memory not required (used only once)                         
!****************                                                       
! read time-range (initialize JPL ephemerides) 
  CALL trange 
! setup the rotation from equatorial to ecliptic and back
  CALL rotpn(roteqec,'EQUM','J2000',0.d0,'ECLM','J2000',0.d0) 
  roteceq=TRANSPOSE(roteqec) 
! restart flag set to default value                                     
  CALL set_restart(.TRUE.) 
! read options
  IF(rhs.EQ.1)THEN                                                          
! force model flags for asteroids OrbFit
     fail=.FALSE. 
! all wrong: reqrd should be false!!!
     CALL rdnint('propag.','ilun',ilun,.TRUE.,found,fail1,fail) 
     CALL rdnint('propag.','imerc',imerc,.TRUE.,found,fail1,fail) 
     CALL rdnint('propag.','iplut',iplut,.TRUE.,found,fail1,fail) 
     CALL rdnint('propag.','irel',irel,.TRUE.,found,fail1,fail) 
! in case selmet is not called                                          
     icrel=irel 
     CALL rdnint('propag.','iast',iast,.TRUE.,found,fail1,fail) 
     CALL rdncha('propag.','filbe',filbe,.TRUE.,found,fail1,fail) 
     CALL rmsp(filbe,ll) 
     IF(ll.LE.0)STOP '**name of asteroid ephem file is wrong**' 
     filbep=filbe(1:ll)//'.bep' 
     filbec=filbe(1:ll)//'.bai' 
     DO ia=1,iast 
        astid(ia)=ia 
     ENDDO
     iatrue=iast
! SEP violation
     CALL rdnlog('propag.','sep_viol',sep_viol,.true.,found,fail1,fail)
     CALL rdnrea('propag.','eta_sep',eta_sep,.true.,found,fail1,fail)
! non gravitational perturbations                                       
     CALL rdnint('propag.','iyark',iyark,.TRUE.,found,fail1,fail) 
     CALL rdnint('propag.','iyarpt',iyarpt,.TRUE.,found,fail1,fail) 
     CALL rdncha('propag.','yardir',yardir,.TRUE.,found,fail1,fail) 
     yarfil=.FALSE. 
     yarini=.FALSE. 
     CALL rdnint('propag.','inongrav',inongrav,.TRUE.,found,fail1,fail)
     CALL rdnrea('propag.','dadt',dadt,.TRUE.,found,fail1,fail)
     CALL rdnrea('propag.','a1ng',a1ng,.TRUE.,found,fail1,fail)
     CALL rdnrea('propag.','a2ng',a2ng,.TRUE.,found,fail1,fail)
     CALL rdnrea('propag.','a3ng',a3ng,.TRUE.,found,fail1,fail)
     CALL rdnrea('propag.','dtdelay',dtdelay,.TRUE.,found,fail1,fail)
  ELSEIF(rhs.EQ.2)THEN
     fail=.false.
     CALL rdnint('propag.','irad',irad,.TRUE.,found,fail1,fail)
     CALL rdnint('propag.','itide',itide,.TRUE.,found,fail1,fail)
     CALL rdnint('propag.','ites',ites,.TRUE.,found,fail1,fail)
     CALL rdnint('propag.','ipla',ipla,.TRUE.,found,fail1,fail)     
     CALL rdncha('propag.','modfile',modfile,.TRUE.,found,fail1,fail)
     IF(ites.ge.2)THEN
        CALL geopot(ites)
     ENDIF
     CALL eamoon_mass
  ELSEIF(rhs.EQ.3)THEN
     WRITE(*,*)'rmodel: for rhs=3 work in progress'
     STOP
  ENDIF 
! close approach control                                                
  CALL rdnint('propag.','iclap',iclap,.TRUE.,found,fail1,fail) 
  fix_mole=.FALSE.
  kill_propag=.FALSE.
  tpplane=.FALSE. ! use MTP unless otherwise stated
  min_dist=.TRUE.
  eprdot=1.d-10 
!                                                                       
  IF(iclap.NE.0)THEN 
! close approach control if requeststed                                 
     CALL rdnint('propag.','npoint',npoint,.TRUE.,found,fail1,fail) 
     CALL rdnrea('propag.','dmea',dmea,.TRUE.,found,fail1,fail) 
     CALL rdnrea('propag.','dmoon',dmoon,.TRUE.,found,fail1,fail) 
     CALL rdnrea('propag.','dmjup',dmjup,.TRUE.,found,fail1,fail) 
     CALL rdnrea('propag.','dmast',dmast,.TRUE.,found,fail1,fail) 
     CALL rdnrea('propag.','dter',dter,.TRUE.,found,fail1,fail) 
  ENDIF
! optical observations
  CALL rdnint('propag.','iaber',iaber,.TRUE.,found,fail1,fail) 
  CALL rdnint('propag.','istat',istat,.TRUE.,found,fail1,fail) 
! radar flag defaults to false          
  radar=.FALSE. 
  radius=-1.d0
! numerical integrator options                                          
  CALL inipro 

! Options for difcor (including outlier rejection)                      
  CALL difini 
  CALL rejini 
! output new/old format for rwo files (only applies in fdiff_cor)
  CALL rdnlog('propag.','output_old_rwo',output_old_rwo,.TRUE.,found,fail1,fail)
! Options for stopping difcor at bizarre orbits; use default            
  IF(rhs.eq.1)THEN
     ecclim=0.d0 
     samin=0.d0 
     samax=0.d0 
     phmin=0.0d0 
     ahmax=0.d0 
     qmax=0.d0
     CALL bizset(ecclim,samin,samax,phmin,ahmax,qmax)
  ELSEIF(rhs.eq.2)THEN
     WRITE(iun_log,*)' rmodel: bizarre controls set by fitdeb.def'
  ELSE
     WRITE(*,*)' rmodel: not supported rhs=', rhs
     STOP
  ENDIF 
! availaility of covarinace matrix is false at start                   
  CALL cov_not_av 
! verbosity is set at the minimum level by default                      
  verb_clo=1 
  verb_pro=1 
  verb_obs=1
  verb_dif=1 
  verb_mul=1 
  verb_rej=1 
  verb_io=1
  verb_moid=1
  verb_prelim=1
  verb_covariance=1
  verb_matrix=1
  verb_2body=1
!                                                                       
  IF(fail)STOP '**** rmodel: abnormal end ****' 
  RETURN 
END SUBROUTINE rmodel

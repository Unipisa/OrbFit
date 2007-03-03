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
!                                                                       
! ================================================================      
! ========INTERFACE============================                         
SUBROUTINE rmodel 
! add a lot of ONLY here
  USE least_squares, ONLY : difini, rejini, output_old_rwo ! it needs only these
  USE propag_state, ONLY : inipro ! calls inipro
  USE tp_trace, ONLY: tpplane
  USE close_app ! sets former closapl options
  USE yark_pert !
  USE planet_masses ! 
  USE force_model
  USE astrometric_observations, ONLY: radius ! to set default value
  USE output_control
  USE reference_systems
  implicit none 
! ========END INTERFACE========================
! no headers:
! close approach options in force_model.mod
! radar quality flag in force_model.mod
! default asteroid radius in astrometric_observations.mod
  character*80 filbe                                               
  logical fail,fail1,found 
  integer ll,ia 
! controls for bizarre orbits                                           
  DOUBLE PRECISION ecclim, samin,samax,phmin,ahmax,qmax 
!****************                                                       
!   static memory not required (used only once)                         
!****************                                                       
! read time-range                                                       
  call trange 
! setup the rotation from equatorial to ecliptic and back
  call rotpn(roteqec,'EQUM','J2000',0.d0,'ECLM','J2000',0.d0) 
  roteceq=TRANSPOSE(roteqec) 
! restart flag set to default value                                     
  CALL set_restart(.true.) 
! read options                                                          
! force model flags                                                     
  fail=.false. 
! all wrong: reqrd should be false!!!
  call rdnint('propag.','ilun',ilun,.true.,found,fail1,fail) 
  call rdnint('propag.','imerc',imerc,.true.,found,fail1,fail) 
  call rdnint('propag.','iplut',iplut,.true.,found,fail1,fail) 
  call rdnint('propag.','irel',irel,.true.,found,fail1,fail) 
! in case selmet is not called                                          
  icrel=irel 
  call rdnint('propag.','iast',iast,.true.,found,fail1,fail) 
  call rdncha('propag.','filbe',filbe,.true.,found,fail1,fail) 
  call rmsp(filbe,ll) 
  if(ll.le.0)STOP '**name of asteroid ephem file is wrong**' 
  filbep=filbe(1:ll)//'.bep' 
  filbec=filbe(1:ll)//'.bai' 
  DO ia=1,iast 
     astid(ia)=ia 
  ENDDO
  iatrue=iast 
! close approach control                                                
  call rdnint('propag.','iclap',iclap,.true.,found,fail1,fail) 
  iorbfit=11
  fix_mole=.false.
  kill_propag=.false.
  tpplane=.false. ! use MTP unless otherwise stated
  min_dist=.true.
  call rdnint('propag.','iaber',iaber,.true.,found,fail1,fail) 
  call rdnint('propag.','istat',istat,.true.,found,fail1,fail) 
  eprdot=1.d-10 
! non gravitational perturbations                                       
  call rdnint('propag.','iyark',iyark,.true.,found,fail1,fail) 
  call rdnint('propag.','iyarpt',iyarpt,.true.,found,fail1,fail) 
  call rdncha('propag.','yardir',yardir,.true.,found,fail1,fail) 
  yarfil=.false. 
  yarini=.false. 
! radar flag defaults to false                                          
  radar=.false. 
  radius=-1.d0
! numerical integrator options                                          
  call inipro 
!                                                                       
  if(iclap.ne.0) then 
! close approach control if requeststed                                 
     call rdnint('propag.','npoint',npoint,.true.,found,fail1,fail) 
     call rdnrea('propag.','dmea',dmea,.true.,found,fail1,fail) 
     call rdnrea('propag.','dmoon',dmoon,.true.,found,fail1,fail) 
     call rdnrea('propag.','dmjup',dmjup,.true.,found,fail1,fail) 
     call rdnrea('propag.','dmast',dmast,.true.,found,fail1,fail) 
     call rdnrea('propag.','dter',dter,.true.,found,fail1,fail) 
  endif
! Options for difcor (including outlier rejection)                      
  call difini 
  CALL rejini 
! output new/old format for rwo files (only applies in fdiff_cor)
  CALL rdnlog('propag.','output_old_rwo',output_old_rwo,.true.,found,fail1,fail)
! Options for stopping difcor at bizarre orbits; use default            
  ecclim=0.d0 
  samin=0.d0 
  samax=0.d0 
  phmin=0.0d0 
  ahmax=0.d0 
  qmax=0.d0
  CALL bizset(ecclim,samin,samax,phmin,ahmax,qmax) 
! availability of covarinace matrix is false at start                   
  CALL cov_not_av 
! verbosity is set at the minimum level by default                      
  verb_clo=1 
  verb_pro=1 
  verb_dif=1 
  verb_mul=1 
  verb_rej=1 
  verb_io=1
!                                                                       
  if(fail) stop '**** rmodel: abnormal end ****' 
  return 
END SUBROUTINE rmodel

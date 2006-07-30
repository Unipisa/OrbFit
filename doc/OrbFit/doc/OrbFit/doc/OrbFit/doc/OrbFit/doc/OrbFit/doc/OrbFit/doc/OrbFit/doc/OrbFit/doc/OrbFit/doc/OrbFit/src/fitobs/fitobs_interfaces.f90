!                                                                       
! ====================================================                  
! FCLOMON2 impact monitoring from FITOBS
! ====================================================  
SUBROUTINE fclomon2(progna,m,obs,obsw,mm1,mm2,tmcla,sigma)
  USE propag_state
  USE tp_trace
  USE output_control
  USE orbit_elements
  USE ret_analysistp
  USE multiple_sol
  USE obssto
  USE force_model, ONLY: masjpl
  USE planet_masses, ONLY: dmea
  USE close_app, ONLY: fix_mole
  USE eval_risk, ONLY: massgiven,givenmass
  IMPLICIT NONE
! =================INPUT=========================================
  character*(6), INTENT(IN) :: progna ! program name
  integer, INTENT(IN) :: m ! observation number
! new data types
  TYPE(ast_obs),DIMENSION(m), INTENT(IN) :: obs
  TYPE(ast_wbsr),DIMENSION(m), INTENT(IN) :: obsw
! min and max index of alternate orbits                                
  INTEGER, INTENT(IN) :: mm1
  INTEGER, INTENT(IN) :: mm2
! search until time (MJD)
  DOUBLE PRECISION, INTENT(IN) :: tmcla 
  DOUBLE PRECISION, INTENT(IN) :: sigma
! workspace por pro_ele 
  TYPE(orbit_elem) el1
  TYPE(orb_uncert) unm1
  CHARACTER*10 mulname
  INTEGER nm,j,le  ! number of virtual objects, index, string length
  INTEGER no                                    ! number of close 
                                                ! approaches found
  INTEGER, PARAMETER :: nox=10000               ! max no close app
  TYPE(tp_point), DIMENSION(nox) :: vas_tr      ! array of tp_traces(global)
  TYPE(tp_point), DIMENSION(nox) :: vas_trl     ! array of tp_traces(of return)
  DOUBLE PRECISION :: dist(nox)
  DOUBLE PRECISION :: dminpos(nox) 
  LOGICAL :: is_min, newflag                    ! local minimum, low  MOID
  INTEGER :: no_risk                   ! number of Virtual impactors found
  INTEGER, PARAMETER :: nshx=10000              ! maximum number of showers
  INTEGER, PARAMETER :: nretx=50000             ! maximum number of returns 
  INTEGER            :: nsho                    ! showers number
  INTEGER            :: isho(nshx)              ! showers index
  INTEGER            :: iret(nretx)             ! returns (trails) index
  INTEGER            :: nret                    ! returns (trails) number
  INTEGER            :: lre,ire,nr,ns
  DOUBLE PRECISION   :: tcat                    ! time of initial conditions
  CHARACTER(LEN=15)  :: despla  ! desired planet 
!=============================OPTIONS=====================================  
  DOUBLE PRECISION :: dt             ! length of showers 
  DOUBLE PRECISION :: tgap           ! time gap desired
  DOUBLE PRECISION :: dmin_used      ! radius of TP
  DOUBLE PRECISION :: dnewton        ! control on distance
  LOGICAL ireq,found
  INTEGER vdifold, vmultold
  CHARACTER*60 comment 
! =======================================================================
! modification August 16, 2005 to input mass from physical observations
  INCLUDE 'parlib.h90' 
  CHARACTER*200 filmass
  LOGICAL ok
  INTEGER iunmass
! ===============================================================
! input options
! ================shower length time============================  
  ireq=.false.
  dt=90.d0  
  comment='shower length time'
  CALL input_rea_opt(progna,'dt',dt,ireq,found,comment,iun_log)
! ================shower gap time============================           
  ireq=.false.
  tgap=30  
  comment='shower gap time'
  CALL input_rea_opt(progna,'tgap',tgap,ireq,found,comment,iun_log)
! ================radius of target plane used===========================
  dmin_used=dmea ! must be equal to the one used in the propagation
! ================control for startup of Newton's method================
! in earth radii                                                        
  ireq=.false.
  dnewton=60.d0 
  comment='newton control value, in Earth radii'
  CALL input_rea_opt(progna,'dnewton',dnewton,ireq,found,comment,iun_log)
! ================control factor on beta================
! to identify entangled                                
  ireq=.false. 
  beta_factor=20.d0
  comment='control factor on beta to identify entangled'
  CALL input_rea_opt(progna,'beta_factor',beta_factor,ireq,found,comment,iun_log)
! ================fix mole by using 2-body inside Earth?==============
  ireq=.false.
  fix_mole=.true.            ! default for fix_mole     
  comment='fix mole by using 2-body inside Earth'
  CALL input_log_opt(progna,'fix_mole',fix_mole,ireq,found,comment,iun_log)
  vdifold=verb_dif
  verb_dif=1
  vmultold=verb_mul
  verb_mul=1
! ============get mass known from other sources, if available=====================
  filmass=dlibd//'/'//name_obj//'.mass'
  CALL rmsp(filmass,le)
  INQUIRE(file=filmass,exist=ok)
  IF(ok)THEN
     CALL filopn(iunmass,filmass(1:le),'old')
     READ(iunmass,*) givenmass
     massgiven=.true.
     CALL filclo(iunmass,' ')
  ELSE
     massgiven=.false.
  ENDIF
! END modification August 16, 2005 to input mass from physical observations
! ===================================================================
  deltasig=delta_sigma ! align stepsize in multiple_sol and tp_trace
  imul0=imi0 ! align index of nominal solution in multiple_sol and tp_trace
  smax=sigma ! to use in computation of minimum possible distance on TP
! copy observations 
  m_m=m
  obs_m(1:m)=obs(1:m)
  obsw_m(1:m)=obsw(1:m)
! setup close app. output
  tpplane=.true. ! use TP plane, not MTP
  REWIND(iuncla) ! cleanup previous close app. records
! propagation and output close approach files
  nm=mm2-mm1+1
  DO j=mm1,mm2
     WRITE(mulname,110)j
 110 FORMAT('mult_',I5)
     CALL rmsp(mulname,le)
     WRITE(iuncla,*)mulname
     CALL cov_avai(unm(j),elm(j)%coo,elm(j)%coord)
     CALL pro_ele(elm(j),tmcla,el1,unm(j),unm1)
     WRITE(*,*) ' VA number ',j,' propagated'
  ENDDO 
  REWIND(iuncla)
! now input tp records
  despla='EARTH'
  CALL masjpl  ! should not matter, unless fclomon2 is executed too early....
  CALL inclolinctp(iun_log,iuncla,despla,vas_tr,no,nox)
  IF(no.le.0)THEN
     WRITE(*,*) ' no close approach to ', despla, ' found up to time MJD ',tmcla
     RETURN
  ENDIF
! shower analysis
  CALL showret3tp(iun_log,no,vas_tr,dt,tgap,isho,nsho,iret,nret)
! main loop on returns                                                  
  ns=1 
  CALL header_rep(iun_log)
  CALL header_rep(0)
  DO 1 nr=1,nret 
! shower counter                                                        
     if(iret(nr).ge.isho(ns+1))ns=ns+1 
! reopen files for filament report                                      
     ire=iret(nr) 
! warning: remember to define iret(nr+1)                                
     lre=iret(nr+1)-iret(nr) 
! analyse return, finding minimum distance etc.                         
     WRITE(iun_log,197)ns,nr,lre,vas_tr(ire)%rindex,vas_tr(ire+lre-1)%rindex 
197  FORMAT('shower no. ', i4,' return no. ',i5,' length ',i5,' from ',f6.1,' to ',f6.1)
! vas_trace  copied in a "return record' vas_traceloc                                  
     CALL arrcut(vas_tr,ire,lre,nox,iun_log,vas_trl) 
     dist(1:lre)=vas_trl(1:lre)%b
     dminpos(1:lre)=vas_trl(1:lre)%minposs
     newflag=.false. 
     IF(lre.gt.1)THEN 
! finding local minima                                                  
        DO 2 j=1,lre 
           is_min=.false. 
           IF(j.eq.1)THEN 
              IF(dist(j).lt.dist(j+1))is_min=.true. 
           ELSEIF(j.eq.lre)THEN 
              IF(dist(j).lt.dist(j-1))is_min=.true. 
           ELSE 
              IF(dist(j).lt.dist(j+1).and.dist(j).lt.dist(j-1))     &
     &                 is_min=.true.                                    
           ENDIF
! if local minimum, check minimum possible                              
           IF(is_min)THEN
              CALL header_rep(iun_log)
              CALL header_rep(0)
              CALL wrireptp(vas_trl(j),vas_trl(1)%rindex,vas_trl(lre)%rindex,-iun_log)   
              IF(dminpos(j).lt.dnewton)THEN 
                 newflag=.true. 
              ENDIF
           ENDIF
2       ENDDO
     ELSE 
! singletons: check if moid is small 
        CALL header_rep(iun_log)
        CALL header_rep(0)  
        CALL wrireptp(vas_trl(1),vas_trl(1)%rindex,vas_trl(lre)%rindex,-iun_log)
        DO j=1,lre 
           IF(dminpos(j).lt.dnewton)newflag=.true.
        ENDDO
     ENDIF
! selected for further analysis?                                        
     IF(newflag)THEN 
! reopen files for newton report                                        
        WRITE(*,*)' return number nr=',nr, ' falsi/newton ' 
!        nnew=nnew+1 
        CALL ret_min(lre,vas_trl,tdt_cat,dnewton,no_risk)
! test output                                                           
!        nrisk=nrisk+no_risk 
     ELSE 
        WRITE(*,*)' return number nr=',nr, 'NO  falsi/newton ' 
     ENDIF
! end loop on returns 
1 ENDDO
  verb_dif=vdifold
  verb_mul=vmultold
END SUBROUTINE fclomon2

! ==================================================================    
! FINOBS                                                                
! ==================================================================    
! Arc  observations input                                               
subroutine finobs(progna,iar,astna0,obs0,nlef,m,obs,obsw,rwofi0  &
     &     ,error_model) 
  USE astrometric_observations
  USE output_control
  implicit none 
! =======INPUT==================================                        
! ========= file names, i/o control ============ 
  character*(6), INTENT(IN) :: progna ! program name
  character*(9), INTENT(IN) :: astna0 ! asteroid name (9 characters)
  integer, INTENT(IN) :: iar ! flag for arcs 1-2 
  CHARACTER*20, INTENT(IN) ::  error_model ! weighing model
! =======OUTPUT================================== 
! ===== main output: observational data =========================== 
  integer, INTENT(IN) :: nlef ! observation number
  integer, INTENT(OUT) :: m ! observation number
! new data types
  TYPE(ast_obs),DIMENSION(nlef), INTENT(OUT) :: obs
  TYPE(ast_wbsr),DIMENSION(nlef), INTENT(OUT) :: obsw
 ! auxiliary output                     
  character*60 rwofi0 ! residuals and weights file name  
  logical obs0 ! successful input flag
  logical change ! change flag
! ============END INTERFACE========================== 
  CHARACTER*60 comment ! for input-options
  logical found,ireq ! input flags                            
  character*60 obsdir0 ! file names, , for observations   
  CHARACTER*18 nam0                                                          
  INTEGER lnam, le, lench  ! string length (after blank removal)  
  integer ll ! loop indexes                                              
  INCLUDE 'sysdep.h90' ! directory char  
! precedence to .obs file with respect to .rwo file: false for use in fitobs
  LOGICAL precob 
  precob=.false.                                             
! ====================================================   
  obs0=.false. ! nothing found yet
  m=0 
! find observations file name
  ireq=.false. 
  IF(iar.eq.1)THEN
     nam0=astna0 
     comment='first object obs. file name' 
     CALL input_cha_opt(progna,'nam0',nam0,ireq,found,comment,iun_log)
  ELSEIF(iar.eq.2 )THEN
     nam0=astna0
     comment='second object obs. file name' 
     CALL input_cha_opt(progna,'namp',nam0,ireq,found,comment,iun_log)
  ENDIF
  CALL rmsp(nam0,lnam)
  IF(lnam.eq.0) RETURN ! blank object does not exist
! find directory of observations    
  ireq=.false.
  obsdir0='obsdata' ! locally called obsdir0 even for second arc    
  if(iar.eq.1)then 
     comment='first object obs. directory name'
     CALL input_cha_opt(progna,'obsdir0',obsdir0,ireq,found,comment,iun_log)
  elseif(iar.eq.2)then 
     comment='second object obs. directory name'
     CALL input_cha_opt(progna,'obsdirp',obsdir0,ireq,found,comment,iun_log)
  endif
  CALL rmsp(obsdir0,le) 
! input data 
  CALL input_obs(obsdir0,nam0,precob,error_model,obs0,obs,obsw,m,   &
     &    iun_log,change)
! compute name of .rwo file                                             
  ll=lench(obsdir0) 
  rwofi0=obsdir0(1:ll)//dircha 
  ll=lench(rwofi0) 
  rwofi0=rwofi0(1:ll)//nam0//'.rwo' 
  CALL rmsp(rwofi0,ll) 
END SUBROUTINE finobs

SUBROUTINE f_gaussdeg8(uniele,name,deforb,defcov,          &
     &     rwfil,obs,obsw,n,error_model,el)
  USE fund_const, ONLY: gms  
  USE astrometric_observations 
  USE least_squares, ONLY: rms_compute
  USE orbit_elements
  USE output_control
  IMPLICIT NONE 
! INPUT observations: new data types
  INTEGER, INTENT(IN) :: n  !number of observations
  TYPE(ast_obs),DIMENSION(n),INTENT(IN) :: obs
  TYPE(ast_wbsr),DIMENSION(n),INTENT(INOUT) :: obsw
  CHARACTER*20, INTENT(IN) ::  error_model ! weighing model
! other input   
  INTEGER,INTENT(IN) :: uniele ! unit for output
  CHARACTER*18, INTENT(IN) :: name !asteroid name 
  CHARACTER*60, INTENT(IN) :: rwfil ! output file for residuals
! OUTPUT
  LOGICAL, INTENT(OUT) :: deforb,defcov ! orbit state on output 
  TYPE(orbit_elem), INTENT(OUT) :: el ! elements  
! END INTERFACE
  DOUBLE PRECISION, DIMENSION(3) :: tobs,alpha3,delta3
  DOUBLE PRECISION tr,dt,dt1 ! mean time, half interval, dist. to center
  INTEGER :: obscod(3), isel(3),j,nselect
! for call to gaussdeg8
  INTEGER nroots,nsol
  LOGICAL fail, debug
  CHARACTER*(20) msg 
  TYPE(orbit_elem), DIMENSION(3) :: elv
  DOUBLE PRECISION rms
! select obs: first, last closest to mean time
     defcov=.false.
  IF(n.le.2)THEN
     deforb=.false.
     RETURN
  ENDIF
  isel(1)=1
  isel(3)=n
  tr= (obs(1)%time_tdt+obs(n)%time_tdt)/2.d0
  dt=(obs(n)%time_tdt-obs(1)%time_tdt)/2.d0
  isel(2)=1
  DO j=2,n
    dt1=abs(obs(j)%time_tdt-tr)
    IF(dt1.lt.dt)THEN
       isel(2)=j
       dt=dt1
    ENDIF
  ENDDO
  IF(isel(2).eq.1.or.isel(2).eq.n)THEN
     WRITE(*,*) 'f_gaussdeg8: logical error in times'
     WRITE(*,*) obs(1:n)%time_tdt
     STOP
  ENDIF     
! copy in array
  DO j=1,3
    tobs(j)=obs(isel(j))%time_tdt
    alpha3(j)=obs(isel(j))%coord(1)
    delta3(j)=obs(isel(j))%coord(2)
    obscod(j)=obs(isel(j))%obscod_i
  ENDDO
! find solution(s)
  debug=.true.
  CALL gaussdeg8(tobs,alpha3,delta3,obscod,elv,nroots,nsol,fail,msg,debug)
  IF(msg.ne.' ')WRITE(*,*)' error message from gaussdeg8=',msg
! assess solutions
  IF(nsol.eq.0)THEN
     deforb=.false.
     RETURN 
  ENDIF
!  DO j=1,nsol
!     rms=rms_compute(obs,obsw,n)
!     WRITE(iun_log,222)j, rms 
!     WRITE(*,222)j, rms 
!222  FORMAT(' preliminary orbit no. ',i3,' RMS of residuals=',f10.4)  
!  ENDDO
! select solution, if more than one
  IF(nsol.eq.1)THEN
     nselect=1
  ELSE
22   WRITE(*,*)' select one solution between 1 and ', nsol
     READ(*,*) nselect
     IF(nselect.lt.1.or.nselect.gt.nsol) GOTO 22
  ENDIF
  deforb=.true.
  el=elv(nselect)  
  WRITE(*,*) ' selected preliminary orbit ', el 

END SUBROUTINE f_gaussdeg8

! Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version OFINOD: December 1, 1997                                      
! Adapted for FITOBS, AM/ZK March 1998                                  
! Adapted for use of VAISALA, ZK November 1998                          
! IOBS added (Feb 10, 1999) MC   
! Fortran 90 version 3.0 A. Milani November 2002   
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                         F G A U S S                           *    
!  *                                                               *    
!  *                Auxiliary routine for FITOBS:                  *    
!  *                 initial orbit determination                   *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! In the present version, only Gauss' method is supported               
!                                                                       
! INPUT:    UNIELE    -  FORTRAN unit for elements                      
!           NAME      -  Object name                                    
!           DEFORB    -  Orbit definition flag                          
!           DEFCOV    -  Orbit covariance definition flag               
!           RWFIL     -  Names of residual/weight files
!           OBS, OBSW -  Input observations and weights                 
!           OBSW%SEL_COORD  -  Selection index (0=don't use; 1=use for fit;   
!                        2=use for fit & Gauss method)                  
!           N         -  Number of observations for each object         
!           IMETH     -  Method to be used: 1=Gauss, 2=Vaisala          
!                                                                       
! OUTPUT:   EL        -  Orbital elements                    
! NOT OUTPUT, but available: residulas inside obsw
!                                                                       
! Input variables DEFORB, DEFCOV, OBSW%SEL_COORD are modified by the routine
!                                                                       
SUBROUTINE f_gauss(uniele,name,deforb,defcov,          &
     &     rwfil,obs,obsw,n,error_model,imeth,el)
  USE fund_const, ONLY: gms  
  USE astrometric_observations 
  USE least_squares, ONLY: rms_compute
  USE orbit_elements
  USE output_control
  IMPLICIT NONE 
! INPUT observations: new data types
  INTEGER, INTENT(IN) :: n  !number of observations
  TYPE(ast_obs),DIMENSION(n),INTENT(IN) :: obs
  TYPE(ast_wbsr),DIMENSION(n),INTENT(INOUT) :: obsw
  CHARACTER*20, INTENT(IN) ::  error_model ! weighing model
! other input   
  INTEGER,INTENT(IN) :: uniele, imeth ! unit for output, method
  CHARACTER*18, INTENT(IN) :: name !asteroid name 
  CHARACTER*60, INTENT(IN) :: rwfil ! output file for residuals
! OUTPUT
  LOGICAL, INTENT(OUT) :: deforb,defcov ! orbit state on output 
  TYPE(orbit_elem), INTENT(OUT) :: el ! elements  
! END INTERFACE

! HEADERS (to be revised)
  INCLUDE 'parcmc.h90'
  INCLUDE 'pariod.h90' ! initial orbit determination parameters
  INCLUDE 'comiod.h90' ! initial orbit determination common
!
  INTEGER ln,j, nused
  INTEGER, EXTERNAL :: lench 
  LOGICAL fail 
  CHARACTER*3 eltype 
  CHARACTER*60 comele 
  DOUBLE PRECISION enne,rms 
  DOUBLE PRECISION, EXTERNAL :: snormd 
  DOUBLE PRECISION elem(6),h,g,gg(6,6),cc(6,6) 
  DOUBLE PRECISION eq(6), telem
  ln=lench(name) 
  IF(imeth.eq.1)THEN 
     WRITE(iun_log,203) name(1:ln) 
203  FORMAT('Automatic method for object ',A,':') 
     iodnm=2 
     iodmet(1)=1 
     iodmen(1)='GAUSS' 
     iodmet(2)=2 
     iodmen(2)='VAISALA' 
  ELSEIF(imeth.eq.2)THEN 
     WRITE(iun_log,201) name(1:ln) 
201  FORMAT('Gauss method for object ',A,':') 
     iodnm=1 
     iodmet(1)=1 
     iodmen(1)='GAUSS' 
  ELSEIF(imeth.eq.3)THEN 
     WRITE(iun_log,202) name(1:ln) 
202  FORMAT('Vaisala method for object ',A,':') 
     iodnm=1 
     iodmet(1)=2 
     iodmen(1)='VAISALA' 
  ENDIF
  CALL io_det(iun_log,rwfil,name,obs,obsw,n,error_model,1,elem,   &
     &     telem,eltype,comele,fail)                          
! **********************************************************************
  deforb=(.NOT.fail) 
  defcov=.false. 
  IF(fail) THEN 
     WRITE(*,204) 
204  FORMAT(2X,'INITIAL ORBIT DETERMINATION FAILED') 
  ELSE 
! change to equinoctal                                                  
     CALL coocha(elem,eltype,gms,eq,'EQU',enne) 
     WRITE(*,220) telem,eq 
220  FORMAT(' preliminary orbit elements for epoch=',f12.4/6f13.7)
     el=undefined_orbit_elem
     el%coo='EQU'
     el%coord=eq
     el%t=telem 
     h=el%h_mag !default undefined values
     g=el%g_mag
! output new elements, multiline format, no header, in .fel file.       
     CALL wro1lh(uniele,'ECLM','J2000','KEP') 
     
     CALL wromlr(uniele,name,elem,'KEP',telem,gg,.false.,cc,.false. &
     &        ,h,g,0.d0)  
     rms=rms_compute(obs,obsw,n)
     WRITE(iun_log,222)rms 
     WRITE(*,222)rms 
222  FORMAT(' RMS of residuals, preliminary orbit=',f10.4) 

  END IF
END SUBROUTINE f_gauss

! ====================================================                  
! FOBPRE predict observations                                           
! ====================================================                  
SUBROUTINE fobpre(icov,ini00,cov00,ok,titnam,filnam,  &
     &   el00,unc00,ids,type1,t1,tut,aobs0,dobs0,t2,dt,astnam)
  USE multiple_sol, ONLY: outmul 
  USE pred_obs
  USE orbit_elements
  USE util_suit
  USE output_control
  USE astrometric_observations
  USE two_states
  USE fund_const
  USE station_coordinates, ONLY: codestat,obscoo
  USE reference_systems, ONLY: observer_position
  IMPLICIT NONE 
! =================INPUT=========================================       
  INTEGER,INTENT(IN) :: icov ! requirements on covariance
  LOGICAL ini00,cov00 ! availability of initial conditions, covariance
  DOUBLE PRECISION, INTENT(IN) :: t1,t2,dt ! observation time, 
! also beginning of ephemerides time, end of ephemerides time, step
  DOUBLE PRECISION tut ! UTC of observation 
  INTEGER ids ! station code
  CHARACTER*(1) type1 ! observation type  
  TYPE(orbit_elem), INTENT(IN) :: el00 ! elements
  TYPE(orb_uncert), INTENT(IN) :: unc00 !covariance and normal matrix 
  CHARACTER*80, INTENT(INOUT) :: titnam ! asteroid name etc. 
  CHARACTER*60, INTENT(IN) :: filnam                     
  DOUBLE PRECISION, INTENT(IN) :: aobs0,dobs0 ! actual obs. for comparison 
  CHARACTER*(*), INTENT(IN) :: astnam ! asteroid name
! ONLY DIRECT OUTPUT
  LOGICAL, INTENT(OUT) :: ok ! necessary data available
! ===== predicted observations ===========                              
! angles (best fit prediction), apparent magnitude
  DOUBLE PRECISION alpha,delta,hmagn
! covariance of the observations                                        
  DOUBLE PRECISION gamad(2,2),axes(2,2),sig(2),gamad1(2,2)
! noise in observations
  DOUBLE PRECISION rmssec, rmsmag
! confidence boundary, line of max variation 
  INTEGER npo, ibv, npo1, npop 
  DOUBLE PRECISION sigma
  DOUBLE PRECISION :: aobs,dobs,adot,ddot ! alpha, delta, proper motion
  DOUBLE PRECISION :: pha,dis,dsun,elo,gallat ! phase, dist. Earth, dist. Sun
  INTEGER  inl ! menu: handling of nonlinearity
  CHARACTER*20 menunam ! menu                                     
  CHARACTER*100 file,fields ! ephemerides output 
  CHARACTER*3 scale 
  INTEGER ln,iuneph,le 
! ===================================================================== 
! options                                                               
! ===================================================================   
! ===================================================================   
! chose handling of nonlinearity                                        
57 IF(icov.ge.3)THEN 
     menunam='prednonl' 
     CALL menu(inl,menunam,3,'How to handle nonlinearity?=',        &
     &         'linear map=',                                           &
     &         '2-body nonlinearity=',                                  &
     &         'full n-body nonlinearity=')
     IF(inl.eq.0)GOTO 57 
  ELSE ! icov=1 for simple obs, icov=2 for use simulated obs
     inl=-1
  ENDIF
  ibv=0
! ===================================================================== 
! check availability of initial conditions (also covariance for icov=2) 
  CALL chereq(icov,ini00,cov00,el00%t,iun_log,ok) 
  IF(.not.ok)RETURN 
! ===================================================================== 
! check availability of JPL ephemerides and ET-UT table                 
  CALL chetim(t1,t1,ok) 
  IF(.not.ok)RETURN 
! ===================================================================== 
! compute prediction; without and with covariance                       
! ===================================================================== 
  IF(icov.eq.1.or.icov.eq.2)THEN 
! ===================================================================== 
! only alpha, delta, magnitude
     inl=1                                          
     CALL predic_obs(el00,ids,t1,type1,             &
     &        alpha,delta,hmagn,inl,                                      &
     &        ADOT0=adot,DDOT0=ddot,PHA0=pha,DIS0=dis,                  &
     &        DSUN0=dsun,ELO0=elo,GALLAT0=gallat)
     IF(icov.eq.1)THEN
        CALL outobc(iun_log,type1,ids,tut,alpha,delta,hmagn,adot,ddot,    &
     &     elo,dis,icov,gamad,sig,axes)
     ELSE
! add second arc, if not there
        IF(.not.obsp)THEN
           obsp=.true.
           astnap='sim'
        ENDIF
        IF(obs0)obstwo=.true.
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
! select noise level
         WRITE(*,*)' give RMS of simulated observation'
         IF(type1.eq.'O')THEN
            WRITE(*,*)' for both alpha, delta in arcsec'
            READ(*,*) rmssec
         ELSEIF(type1.eq.'R')THEN
            WRITE(*,*)' for distance, in km'
            READ(*,*) rmssec
         ELSEIF(type1.eq.'V')THEN
            WRITE(*,*)' for range rate, in km/day'
            READ(*,*) rmssec
         ENDIF
         rmsmag=0.7d0
! add observation to second arc
         mall=mall+1
         mp=mp+1   
         CALL obs_simul(type1,t1,tut,astnam,ids,rmssec,rmsmag,alpha,delta,   &
      &                hmagn,obs(mall),obsw(mall))
         IF(type1.eq.'R'.or.type1.eq.'V')CALL radar_ob(obs(1:mall)%type,mall)
      ENDIF
   ELSEIF(icov.eq.3)THEN 
! ===================================================================== 
! alpha, delta, magnitude, covariance and ellipse of confidence         
     CALL predic_obs(el00,ids,t1,type1,             &
     &        alpha,delta,hmagn,inl,                                    &
     &        UNCERT=unc00,GAMAD=gamad,SIG=sig,AXES=axes,          &
     &        ADOT0=adot,DDOT0=ddot,PHA0=pha,DIS0=dis,                  &
     &        DSUN0=dsun,ELO0=elo,GALLAT0=gallat)
     CALL outobc(iun_log,type1,ids,tut,alpha,delta,hmagn,adot,ddot,    &
     &     elo,dis,icov,gamad,sig,axes)                                 
! ===================================================================   
! generation of sky epehemrides                                         
  ELSEIF(icov.eq.6)THEN 
! check availability of JPL ephemerides and ET-UT table for entire time 
     CALL chetim(t1,t2,ok) 
     IF(.not.ok)RETURN 
     IF(nint(abs(t2-t1)/dt).gt.500)THEN 
        write(*,*)' Too many ephemerides points:',                  &
     &           nint(abs(t2-t1)/dt)                                    
        write(*,*)'Select a time interval and span to ',            &
     &           'ensure that there are fewer than 500 points.'         
     ELSE 
! open ephemerides file in current directory                            
        file=astnam//'.eph' 
        CALL rmsp(file,ln) 
        CALL filopn(iuneph,file(1:ln),'unknown') 
        fields='cal,mjd,coord,mag,elong,glat,r,delta,appmot,skyerr' 
        scale='UTC' 
        CALL ephemc(iuneph,el00,unc00,.true.,t1,t2,dt,ids,scale,fields)
        CALL filclo(iuneph,' ') 
        WRITE(*,*)' Generated ephemeris in file: ',file(1:ln) 
     ENDIF
! ===================================================================== 
  ELSEIF(icov.eq.4.or.icov.eq.5)THEN 
! ===================================================================== 
! alpha, delta, magnitude, covariance and confidence boundary;          
! input specification of set of points                                  
     CALL asscbd(iun_log,npoinx,npo,sigma,ibv) 
! ===================================================================== 
! compute prediction, boundary                                          
     CALL predic_obs(el00,ids,t1,type1,             &
     &        alpha,delta,hmagn,inl,                                    &
     &        unc00,sigma,npo,ibv,gamad,sig,axes,npo1,                   &
     &        adot,ddot,pha,dis,dsun,elo,gallat)
     CALL outobc(iun_log,type1,ids,tut,alpha,delta,hmagn,adot,ddot,    &
     &        elo,dis,icov,gamad,sig,axes)                                 
     IF(npo1.le.0)THEN 
        WRITE(*,*)'fobpre: no elliptic orbits ',npo1 
        RETURN 
     ENDIF
     IF(ibv.eq.1)THEN 
! confidence boundary; one point added to close line                    
        al_m(npo1+1)=al_m(1) 
        de_m(npo1+1)=de_m(1) 
        hmag_m(npo1+1)=hmag_m(1) 
        el_m(1:6,npo1+1)=el_m(1:6,1) 
!            IF(inl.eq.3)disv(npo+1)=disv(1) 
        npop=npo1+1 
     ELSEIF(ibv.eq.2)THEN 
! line of variations                                                    
        npop=npo1 
     ENDIF
! if no observation is given, use the nominal marked with a cross       
     IF(icov.eq.4)THEN 
        aobs=alpha 
        dobs=delta 
     ELSE
        aobs=aobs0
        dobs=dobs0
     ENDIF
! ===================================================================== 
! output observation, apparent motion, confidence boundary              
     CALL outmul(titnam,filnam,tut,sigma,alpha,delta,               &
     &              al_m,de_m,hmag_m,1,npop,1,icov-3,aobs,dobs,type1)
  ENDIF
END SUBROUTINE fobpre                                      
! ===================================================================   
! FINELE                                                                
! ===================================================================   
! input of initial conditions                                           
SUBROUTINE finele(progna,iar,astna0,el0,ini0,cov0,unc0)
  USE fund_const                                
  USE orbit_elements
  USE output_control           
  IMPLICIT NONE 
! ==========INPUT===================                                    
! file names, i/o control                                      
  CHARACTER*18 astna0 ! asteroid full name                               
  CHARACTER*6 progna ! program name 
  INTEGER iar ! flag for arcs 1-2 
! ==========OUTPUT==================                                    
! epoch time (MJD), elements (equinoctal), absolute magnitude, opp.effec
  TYPE(orbit_elem), INTENT(INOUT) :: el0 !remains defined!
  LOGICAL ini0 ! successful input flag 
  TYPE(orb_uncert), INTENT(OUT) :: unc0 ! covariance matrices 
  LOGICAL cov0 ! covariance available
! =========END INTERFACE=============                                   
! asteroid name for elements (18 characters)                            
  CHARACTER*18 namel0(1) 
  INTEGER lnam                                                   
  LOGICAL ireq, found ! logical input flags 
  CHARACTER*60 comment                     
  INTEGER le,lench ! length of names
! variables for reading routines  
  CHARACTER*60 elefi0(1)   ! file name for elements    
  DOUBLE PRECISION elem(6) 
  DOUBLE PRECISION mass(1) 
  CHARACTER*(80) comele(1)                                     
  LOGICAL ini(1),cov(1) 
  CHARACTER*3 coox,coo(1) 
  DOUBLE PRECISION t(1),h(1),gma(1) 
  DOUBLE PRECISION enne, eq0(6), gam(6,6), c(6,6) !elements, covariances 
! ====================================                                  
  lnam=lench(astna0) 
! this arc does not exist                                               
  IF(lnam.eq.0)RETURN 
! ====================================================                  
! check that the asteroid name has been read 
  ireq=.false. 
  namel0(1)=astna0
  IF(iar.eq.1)THEN 
     comment='first arc asteroid name'
     CALL input_cha_opt(progna,'namel0',namel0(1),ireq,found,comment,iun_log)
  ELSEIF(iar.eq.2)THEN 
     comment='second arc asteroid name'
     CALL input_cha_opt(progna,'namelp',namel0(1),ireq,found,comment,iun_log)
  ENDIF
  CALL rmsp(namel0(1),lnam) 
! read the name of the elements file, inquire                           
  ireq=.false. 
  elefi0(1)='ast.cat' 
  IF(iar.eq.1)THEN 
     comment='first arc elements file'
     CALL input_cha_opt(progna,'elefi0',elefi0(1),ireq,found,comment,iun_log)
  ELSEIF(iar.eq.2)THEN 
     comment='second arc elements file'
     CALL input_cha_opt(progna,'elefip',elefi0(1),ireq,found,comment,iun_log)
  ENDIF
  CALL rmsp(elefi0(1),le) 
  INQUIRE(file=elefi0(1),exist=found) 
  IF(found)THEN 
     ini(1)=.false. 
     CALL rdelem(iun_log,namel0(1),1,elefi0(1),1,ini,cov,             &
     &        coo,t,elem,gam,c,mass,h,gma,comele)                       
! error case      
     ini0=ini(1)                                      
     IF(.not.ini0)THEN 
        write(*,*)'asteroid ',namel0,' not found in ',elefi0 
        write(iun_log,*)'asteroid ',namel0,' not found in ',elefi0 
        RETURN 
     ENDIF
     cov0=cov(1) 
     coox=coo(1)
     el0=undefined_orbit_elem
     el0%coo=coo(1) 
     el0%t=t(1) 
     IF(h(1).gt.0.d0)THEN
        el0%mag_set=.true.
     ELSE
        el0%mag_set=.false.
     ENDIF
     el0%h_mag=h(1) 
     el0%ndim=6
     IF(gma(1).lt.-1.d6)THEN 
        el0%g_mag=0.15d0
     ELSE 
        el0%g_mag=gma(1)
     ENDIF
     el0%coord=elem
! initial conditions found                                              
     WRITE(*,*)el0
  ELSE 
     WRITE(*,*) 'File ',elefi0(1)(1:le),' not found!' 
  ENDIF
  IF(cov0)THEN
     unc0=undefined_orb_uncert
     unc0%g=gam
     unc0%c=c
     unc0%succ=.true.
  ENDIF
END SUBROUTINE finele
                                         
! ===================================================================== 
! WRIEQU (write initial conditions, equinoctal)                         
! ===================================================================== 
subroutine wriequ(iun,astna0,t0,eq0) 
  implicit none 
  double precision t0,eq0(6) 
  integer iun 
  character*18 astna0 
! initial conditions found                                              
  write(*,108)astna0,t0 
108 format(1x,a18,' initial elem (a,h,k,p,q,lam), epoch=',f8.1) 
  write(iun,104) eq0 
  write(*,104) eq0 
104 format(6f13.7) 
  write(iun,*)' ' 
  return 
END SUBROUTINE wriequ

! ========================================                              
! FIDENT compute identification norm                                    
! =======================================                               
SUBROUTINE fident(id_dim,cov0,covp,el0,elp,unc0,uncp,elid,ff)  
  USE orbit_elements 
  USE output_control
  IMPLICIT NONE 
  LOGICAL, INTENT(OUT):: ff !failure of identification proposal
  LOGICAL fail 
  INTEGER,INTENT(IN) :: id_dim ! method control (normally 6) 
  TYPE(orbit_elem), INTENT(IN) :: el0,elp
  LOGICAL,INTENT(IN) :: cov0,covp
  TYPE(orb_uncert), INTENT(IN) :: unc0,uncp
  TYPE(orbit_elem), INTENT(OUT) :: elid
! two input elements, their uncertainties
  DOUBLE PRECISION eq0(6),eqp(6)
  DOUBLE PRECISION c0(6,6),cp(6,6),g0(6,6),gp(6,6) 
! identification elements                                          
  DOUBLE PRECISION eqf(6) 
! similarity norms   
  DOUBLE PRECISION d2,da2,dista,dq,dqalt 
! determinants and eigenvalues                                          
  DOUBLE PRECISION detc2,detc5,detc6,eigen5(5),eigen6(6) 
! ========================================                              
! check requirements 
  ff=.false.                                                   
  IF(.not.cov0)THEN 
     WRITE(*,*)' covariance matrix for arc 1 not ready' 
     RETURN 
  ENDIF
  IF(.not.covp)THEN 
     WRITE(*,*)' covariance matrix for arc 2 not ready' 
     RETURN 
  ENDIF
  IF(el0%t.ne.elp%t)THEN 
     WRITE(*,*)' to=',el0%t,' tp=', elp%t
     WRITE(*,*)' elements and covariances must be available' 
     WRITE(*,*)' for the same time, use propagation first' 
     RETURN 
  ENDIF
!  IF(el0%coo.ne.'EQU'.or.elp%coo.ne.'EQU')THEN
!     WRITE(*,*)' identification guess only in EQU, not in ', el0%coo, elp%coo
!     RETURN
!  ENDIF
  eq0=el0%coord
  eqp=elp%coord
  g0=unc0%g
  c0=unc0%c
  gp=uncp%g
  cp=uncp%c
! selection of algorithm                                                
  IF(id_dim.eq.5)THEN 
! identification by 5x5 matrix                                          
     CALL idno5(eq0,eqp,g0,c0,gp,cp,                                &
     &    d2,da2,dq,dqalt,dista,detc2,detc5,eqf,fail,eigen5)            
  ELSEIF(id_dim.eq.6)THEN 
! identification by 6x6 matrix                                          
     CALL idno6(eq0,eqp,g0,c0,gp,cp,                                &
     &    d2,da2,dq,dqalt,dista,detc2,detc6,eqf,fail,eigen6)            
  ELSE 
     WRITE(*,*)' fident: id_dim=',id_dim,' not understood' 
     RETURN 
  ENDIF
! output                                                                
  IF(fail)THEN 
     CALL tee(iun_log,' FAILED IDENTIFICATION ALGORITHM=') 
     ff=.false. 
  ELSE 
     WRITE(*,*)' INITIAL COND. AVAILABLE, NOW USE DIFF. CORR.' 
     ff=.true. 
! store as proposed identification 
     elid=el0 ! magnitude from the first...
     elid%coord=eqf
  ENDIF
! penalties                                                             
  WRITE(*,*)' d2, d2alt,  dq,   dqalt,    dista' 
  WRITE(*,194)d2,da2,dq,dqalt,dista 
194 FORMAT(4(f13.4,1x),f10.6) 
  WRITE(iun_log,*)' d2, d2alt,  dq,   dqalt,    dista' 
  WRITE(iun_log,194)d2,da2,dq,dqalt,dista 
! proposed elements                                                     
  WRITE(iun_log,208)elid%t 
208 FORMAT(' ident. elem (a,h,k,p,q,lam), epoch=',f8.1) 
  WRITE(*,208)elid%t 
  WRITE(*,104)eqf 
  WRITE(iun_log,104)eqf 
104 FORMAT(6f13.7) 
! store as proposed identification                                      
END SUBROUTINE fident








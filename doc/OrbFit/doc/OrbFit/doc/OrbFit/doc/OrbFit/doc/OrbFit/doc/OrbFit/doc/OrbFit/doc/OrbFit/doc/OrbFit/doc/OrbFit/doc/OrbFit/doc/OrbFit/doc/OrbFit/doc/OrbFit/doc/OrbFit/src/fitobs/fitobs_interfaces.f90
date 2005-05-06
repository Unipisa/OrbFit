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
  USE fund_const    
  USE astrometric_observations 
  USE least_squares
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
SUBROUTINE fobpre(icov,ini0,cov0,iun8,ok,titnam,filnam,  &
     &   el0,unc0,ids,type1,t1,tut,aobs0,dobs0,t2,dt,astnam)
  USE multiple_sol, ONLY: outmul 
  USE pred_obs
  USE orbit_elements
  USE util_suit
  USE output_control
  IMPLICIT NONE 
! =================INPUT=========================================       
  INTEGER,INTENT(IN) :: icov ! requirements on covariance
  LOGICAL ini0,cov0 ! availability of initial conditions, covariance
  INTEGER, INTENT(IN) ::  iun8 ! output unit for covariance 
  DOUBLE PRECISION, INTENT(IN) :: t1,t2,dt ! observation time, 
! also beginning of ephemerides time, end of ephemerides time, step
  DOUBLE PRECISION tut ! UTC of observation 
  INTEGER ids ! station code
  CHARACTER*(1) type1 ! observation type  
  TYPE(orbit_elem), INTENT(IN) :: el0 ! elements
  TYPE(orb_uncert), INTENT(IN) :: unc0 !covariance and normal matrix 
  CHARACTER*80, INTENT(IN) :: titnam ! asteroid name etc. 
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
! confidence boundary, line of max variation 
  INTEGER npo, ibv, npo1, npop 
  DOUBLE PRECISION sigma
  DOUBLE PRECISION :: aobs,dobs,adot,ddot ! alpha, delta, proper motion
  DOUBLE PRECISION :: pha,dis,dsun,elo,gallat ! phase, dist. Earth, dist. Sun
  INTEGER  inl ! menu: handling of nonlinearity
  CHARACTER*20 menunam ! menu                                     
  CHARACTER*100 file,fields ! ephemerides output 
  CHARACTER*3 scale 
  INTEGER ln,iuneph 
! ===================================================================== 
! options                                                               
! ===================================================================   
! ===================================================================   
! chose handling of nonlinearity                                        
57 IF(icov.ge.2)THEN 
     menunam='prednonl' 
     CALL menu(inl,menunam,3,'How to handle nonlinearity?=',        &
     &         'linear map=',                                           &
     &         '2-body nonlinearity=',                                  &
     &         'full n-body nonlinearity=')
     IF(inl.eq.0)GOTO 57 
  ELSE
     inl=-1
  ENDIF
  ibv=0
! ===================================================================== 
! check availability of initial conditions (also covariance for icov=2) 
  CALL chereq(icov,ini0,cov0,el0%t,iun_log,iun8,ok) 
  IF(.not.ok)RETURN 
! ===================================================================== 
! check availability of JPL ephemerides and ET-UT table                 
  CALL chetim(t1,t1,ok) 
  IF(.not.ok)RETURN 
! ===================================================================== 
! compute prediction; without and with covariance                       
! ===================================================================== 
  IF(icov.eq.1)THEN 
! ===================================================================== 
! only alpha, delta, magnitude
     inl=1                                          
     CALL predic_obs(el0,ids,t1,type1,             &
     &        alpha,delta,hmagn,inl,                                      &
     &        ADOT0=adot,DDOT0=ddot,PHA0=pha,DIS0=dis,                  &
     &        DSUN0=dsun,ELO0=elo,GALLAT0=gallat)
     CALL outobc(iun_log,type1,ids,tut,alpha,delta,hmagn,adot,ddot,    &
     &     elo,dis,icov,gamad,sig,axes)                                 
  ELSEIF(icov.eq.2)THEN 
! ===================================================================== 
! alpha, delta, magnitude, covariance and ellipse of confidence         
     CALL predic_obs(el0,ids,t1,type1,             &
     &        alpha,delta,hmagn,inl,                                    &
     &        UNCERT=unc0,GAMAD=gamad,SIG=sig,AXES=axes,          &
     &        ADOT0=adot,DDOT0=ddot,PHA0=pha,DIS0=dis,                  &
     &        DSUN0=dsun,ELO0=elo,GALLAT0=gallat)
     CALL outobc(iun_log,type1,ids,tut,alpha,delta,hmagn,adot,ddot,    &
     &     elo,dis,icov,gamad,sig,axes)                                 
! ===================================================================   
! generation of sky epehemrides                                         
  ELSEIF(icov.eq.5)THEN 
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
        CALL ephemc(iuneph,el0,unc0,.true.,t1,t2,dt,ids,scale,fields)
        CALL filclo(iuneph,' ') 
        WRITE(*,*)' Generated ephemeris in file: ',file(1:ln) 
     ENDIF
! ===================================================================== 
  ELSEIF(icov.eq.3.or.icov.eq.4)THEN 
! ===================================================================== 
! alpha, delta, magnitude, covariance and confidence boundary;          
! input specification of set of points                                  
     CALL asscbd(iun_log,npoinx,npo,sigma,ibv) 
! ===================================================================== 
! compute prediction, boundary                                          
     CALL predic_obs(el0,ids,t1,type1,             &
     &        alpha,delta,hmagn,inl,                                    &
     &        unc0,sigma,npo,ibv,gamad,sig,axes,npo1,                   &
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
     IF(icov.eq.3)THEN 
        aobs=alpha 
        dobs=delta 
     ELSE
        aobs=aobs0
        dobs=dobs0
     ENDIF
! ===================================================================== 
! output observation, apparent motion, confidence boundary              
     CALL outmul(titnam,filnam,tut,sigma,alpha,delta,               &
     &              al_m,de_m,hmag_m,1,npop,1,icov-2,aobs,dobs,type1)
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
  IF(el0%coo.ne.'EQU'.or.elp%coo.ne.'EQU')THEN
     WRITE(*,*)' identification guess only in EQU, not in ', el0%coo, elp%coo
     RETURN
  ENDIF
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








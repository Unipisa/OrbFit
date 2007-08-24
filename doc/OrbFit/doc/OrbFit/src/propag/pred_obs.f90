! ==============MODULE PRED_OBS======================
! replaces both predict_obs.f90 and obs_compute.f90 
! (but also alph_del.f90 is required, may be added here)
! headers npoint and phase are eliminated 
! OUT OF MODULE
!                  aber1      aberration 
!                  outobc     observation output
! ROUTINES   
!                  predic_obs  replaces preobX
!                  predic_obs2 for attributables
!                  pre_obs_att for attribute only
!                  alph_del2    optical obs. with proper motion   public       
!                    oss_dif2       "   
!                  alph-del     optical obs.                      public
!                    oss_dif
!                  r_rdot      radar observations                 public       
!                    deltau        "                                    
!                    deldop1       "                                    
!                    deldop2       "   
! HEADERS and MODULES
!     pred_obs.o: \
!	../suit/ASTROMETRIC_OBSERVATIONS.mod \
!	../suit/FUND_CONST.mod \
!	../suit/OUTPUT_CONTROL.mod \
!	../suit/REFERENCE_SYSTEMS.mod \
!	../suit/STATION_COORDINATES.mod \
!	force_model.o \
!	propag_state.o 

MODULE pred_obs
USE fund_const
USE output_control
IMPLICIT NONE
PRIVATE

! PUBLIC ROUTINES
PUBLIC predic_obs, alph_del2, r_rdot, alph_del, pre_obs_att, predic_obs2

! PUBLIC DATA

! ===============================================================
! former phase header 
! galactic pole alpha=12h 51.3m, delta=27d 7.0min
double precision algal,degal,gax,gay,gaz
parameter (algal=3.36543113015807d0)
parameter (degal=0.47327511549913d0)
parameter (gax=-0.86787505968543d0)
parameter (gay=-0.19757464383054d0)
parameter (gaz=0.45580384036475d0)
! from http://www.seds.org/~spider/spider/ScholarX/coords.html#galactic
! The galactic north pole is at RA = 12:51.4, Dec = +27:07 (2000.0)
! the galactic center at RA = 17:45.6, Dec = -28:56 (2000.0). 
! The inclination of the galactic equator to Earth''s equator is thus 62.9 deg. 
! The intersection, or node line of the two equators is at 
! RA = 18:51.4, Dec = 0:00 (2000.0), and at l = 33 deg, b=0.

! ===============================================================
! former npoint header
! common to handle data on multiple observations
integer,parameter :: npoinx=4000
DOUBLE PRECISION :: al_m(npoinx),de_m(npoinx),hmag_m(npoinx) ! observation
DOUBLE PRECISION :: el_m(6,npoinx) ! line of elements 
PUBLIC npoinx,al_m,de_m,hmag_m,el_m
! phase, distance to Earth, distance to Sun (to compute magnitude)
! elongation, galactic latitude, apparent motion
double precision phav(npoinx),disv(npoinx),dsunv(npoinx)
double precision elov(npoinx)
double precision gallav(npoinx),adotv(npoinx),ddotv(npoinx)
! but the above need to be public???
PUBLIC disv

CONTAINS
! ===================================================================== 
! PREDIC_OBS   fortran90 version, A. Milani, January 2003
! ===================================================================== 
! method: call with minimum set of arguments results in simple
!         nominal prediction; with optional arguments for linear
!         confidence boundary output in optional args.
!         for nonlinear boundary, the multiple data structure is
!         in the public data of the module, but is not accessed if the 
!         USE statement is restricted to the routine pred_obs  
! problems:
!          1: adot,ddot and other data from former phase header are 
!              part of the output, supplied by alph_del2 or similar
!  input: minimum
!          el    = orbital elements, including
!                     coordinate type EQU, KEP, COM, CAR, ATT
!                     epoch time (MJD,TDT)
!                     orbital elements vector 
!                     absolute magnitude                                   
!                     opposition effect coefficient
!          idsta = station code (integer)                                      
!          tobs  = prediction time (MJD,TDT)                                
!          type  = observation type (O for optical, R=V for radar; S not handled)
! input optional
!          uncert = orbital uncertainty, including 
!                   covariance at epoch time 
!                    normal matrix at epoch time (should be the inverse)
!          sigma = level of the confidence boundary in RMS values       
!          npo   = number of points in the boundary                     
!          ibv   = Type of depiction                                    
!                       =0 for automatic selection                      
!                       =1 for confidence boundary                      
!                       =2 for line of maximum variation                
!          inl   = handling of nonlinearity                             
!                       =0 for automatic selection                      
!                       =1 for linear ellipse                           
!                       =2 for 2-Body nonlinear propagation of covarianc
!                       =3 for n-body nonlinear propagation of covarianc
!                                                                       
!  output: minimum
!          alpha = right ascension (equatorial J2000), radians          
!          delta = declination (equatorial J2000), radians              
!          hmagn = apparent magnitude, as predicted, from h and g given 
!  output optional :
!          gamad = covariance matrix of observations alpha, delta       
!          sig   = sqrt(eigenvalues) of gamad                           
!          axes  = the eigenvectors of gamad are the columns of this mat
!  hidden output :
!          npo1  = number of output dta points (could be less than npo) 
!          al(npo1),de(npo1) points on the confidence boundary          
!                   (difference with respect to best prediciton, radians
!          elm(npo1) alternate elements for observation time            
!                                                                       
!  In the linear approximation, the ellipse of confidence has principal axes 
!          along axes; the semiaxes lenghts are sig                     
!  In the nonlinear approximation, the boundary is a map of the         
!          confidence ellipse in the elements space                     
!                                                                       
!  WARNING: the magnitudes are often very poorly predictable            
! ============INTERFACE=================================================
SUBROUTINE predic_obs(el,idsta,tobs,type,       &
     &    alpha,delta,hmagn,inl,                                    &
     &    uncert,sigma,npo,ibv,gamad,sig,axes,npo1,               &  
     &    adot0,ddot0,pha0,dis0,dsun0,elo0,gallat0,ecllat0,twobo,elev0,elsun0)
  USE station_coordinates 
  USE orbit_elements                 
! ============= input ==================================================
  TYPE(orbit_elem), INTENT(IN) :: el ! elements
  DOUBLE PRECISION, INTENT(IN) :: tobs ! elements, epoch time, obs.time
  INTEGER,INTENT(IN) ::  idsta ! station code   
  CHARACTER*1, INTENT(IN) :: type ! observation type
! flag for 2-body approximation; must be .false. for full n-body computa
  LOGICAL,INTENT(IN),OPTIONAL :: twobo 
! ============= output =================================================
  DOUBLE PRECISION, INTENT(OUT) :: alpha,delta,hmagn ! best fit obs., apparent magnitude
! ======optional input ================================
  TYPE(orb_uncert), INTENT(IN),OPTIONAL :: uncert ! covariance, normal matrices
  INTEGER, INTENT(INOUT) :: inl ! nonlinearity; used in output when 0 in input 
  INTEGER, INTENT(IN),OPTIONAL :: npo ! no points
  INTEGER, INTENT(INOUT),OPTIONAL ::  ibv ! conf.bd/LOV; used in output when 0 in input
  DOUBLE PRECISION, INTENT(IN),OPTIONAL :: sigma ! sigma value for the boundary
! ======optional output =================================================
  DOUBLE PRECISION,INTENT(OUT),OPTIONAl :: adot0,ddot0 ! proper motion
  DOUBLE PRECISION,INTENT(OUT),OPTIONAl :: gamad(2,2),axes(2,2),sig(2) ! covariance on sky plane
  DOUBLE PRECISION,INTENT(OUT),OPTIONAl :: pha0,dis0,dsun0,elo0,gallat0,ecllat0 ! phase, distance to Earth, distance to Sun 
                                             ! elongation, galactic latitude
  DOUBLE PRECISION,INTENT(OUT),OPTIONAl :: elev0,elsun0 ! elevation of obj, of Sun
! points on the confidence boundary (difference w.r. to alpha,delta)    
! WARNING! the output number of points is npo1.le.npo; this beacuse hyperbolic points are discarded 
  INTEGER,INTENT(OUT),OPTIONAL :: npo1 
! ============END INTERFACE=============================================
  DOUBLE PRECISION :: gameq(6,6),ceq(6,6) ! covariance, normal matrices 
! partial derivatives of alpha, delta, w.r. to elements (by columns)    
  DOUBLE PRECISION daddet(6,2),dummy(6) 
! second derivatives of alpha, delta, w.r. to elements (not used)       
!   DOUBLE PRECISION ddade(6,6),dddde(6,6) 
! ===================================================================   
! orthonormal basis, matrix defining the plane of the ellipse           
  DOUBLE PRECISION v(6,6),ceicel(4,2) 
! transformation matrix between the two planes                          
  DOUBLE PRECISION b(2,2)
  TYPE(orbit_elem) :: elv
! number of full revolutions around the sky                             
  INTEGER ng,nrev 
! functions                                                             
  DOUBLE PRECISION appmag
! for astronomical unit in km from fund_const
! elongation,distance to Earth, distance to Sun (to compute magnitude)  
  DOUBLE PRECISION adot,ddot,pha,dis,rdot,dsun,elo,gallat, ecllat, elev, elsun
  double precision obs4(4) ! observations: alpha, delta in RAD  alphadot, deltadot in RAD/day
  double precision dobde(4,6) ! partial derivatives of obs, w.r. to asteroid elements
! for r_rdot
  DOUBLE PRECISION posr(3) ! B-F position of reciever, assumed to be aldo transmitter
  CHARACTER*1 tech ! assumed center of mass correction applied
  CHARACTER(LEN=16)  :: stname 
  INTEGER idstarad  ! with encoded trasnmitter and receiver    
! ===================================================================   
! constant of gravitation, trigonometric constants from fund_const.mod 
! temporaries, indexes                                                  
  DOUBLE PRECISION dal,ddl,maxsig,minsig,allin,delin 
  INTEGER n,inlu,ibvu,ider
  LOGICAL cov ! true if covariance computations are possible
  LOGICAL twobo1 
!****************                                                       
!   static memory not required                                          
!****************                                                       
  IF(PRESENT(twobo))THEN
     twobo1=twobo
  ELSE
     twobo1=.false.
  ENDIF
! ===================================================================== 
  if((type.eq.'R'.or.type.eq.'V').and.(inl.eq.2.or.twobo1))then 
     WRITE(*,*)' predic_obs: mixing of radar and two-body '//      &
          &        'approximation not permitted'                             
     RETURN 
  endif
! simple prediction without derivatives
  cov=PRESENT(uncert).and.PRESENT(sig).and.PRESENT(axes)
  IF(cov)cov=(cov.and.uncert%succ)
  IF(cov)THEN
    gameq=uncert%g
    ceq=uncert%c
    ider=1
  ELSE
    ider=0
  ENDIF
! ===================================================================== 
! compute observation; derivatives (of order 1) if required                
  if(type.eq.'O')then 
     CALL alph_del2 (el,tobs,idsta,obs4,ider,dobde,      &
            &   pha,dis,rdot,dsun,elo,gallat,ecllat,twobo1,elev,elsun) 
     alpha=obs4(1)
     delta=obs4(2)
     IF(ider.gt.0)THEN
        daddet(1:6,1)=dobde(1,1:6)
        daddet(1:6,2)=dobde(2,1:6)
     ENDIF
! store true apparent motion, etc.  
     IF(PRESENT(adot0))adot0=obs4(3)
     IF(PRESENT(ddot0))ddot0=obs4(4)    
     IF(PRESENT(pha0))pha0=pha 
     IF(PRESENT(dis0))dis0=dis 
     IF(PRESENT(dsun0))dsun0=dsun 
     IF(PRESENT(elo0))elo0=elo 
     IF(PRESENT(gallat0))gallat0=gallat 
     IF(PRESENT(ecllat0))ecllat0=ecllat
     IF(PRESENT(elev0))elev0=elev
     IF(PRESENT(elsun0))elsun0=elsun
! compute apparent magnitude at time of observation                     
     hmagn=appmag(el%h_mag,el%g_mag,dsun,dis,pha) 
  elseif(type.eq.'R'.or.type.eq.'V')then 
     tech='c'  ! center of mass correction assumed
     CALL obscoo(idsta,posr,stname) ! transmitter and receiver assumed both idsta
     idstarad=idsta*10000+idsta
     CALL r_rdot (el,tobs,idstarad,tech,posr,posr,alpha,delta,       &
     &       daddet(1,1), daddet(1,2),ider) 
     hmagn=0.d0 
  ELSE
! ouput adot...... 
     WRITE(*,*)' predic_obs: this observation type not supported ',type 
     RETURN                                                         
  ENDIF
  IF(.not.cov) RETURN
! ===================================================================== 
! *******************fix for infamous polar bug*********************
  IF(type.eq.'O'.or.type.eq.'S')THEN
     daddet(1:6,1)=daddet(1:6,1)*cos(delta)
  ENDIF
! compute ellipse of covariance of alpha*cos(delta),delta
! *******************end fix polar bug******************************
! compute ellipse of covariance of alpha,delta                          
  CALL ellips(daddet,gameq,sig,axes,gamad)
! confidence boundary?
  IF(.not.PRESENT(ibv).or..not.PRESENT(npo))RETURN
  IF(.not.PRESENT(uncert))THEN
     WRITE(*,*)' predic_obs: error in arguments, c0 missing, inl=',inl
     RETURN
  ENDIF
! If inl=0 then use automatic selection method                          
  IF(inl.eq.0)THEN                                                  
     maxsig=max(sig(1),sig(2))*degrad 
     if(maxsig.le.1.0d0)then                                        
        inlu=1                                                       
! Is it safe to use two body if we may have close approach? NO           
!         elseif(maxsig .le. 5.d0)then                                  
!            inl=2                                                      
     else                                                           
        inlu=3                                                       
     endif
  else
     inlu=inl                                                        
  endif
  inl=inlu ! to output choice done when in auto mode
! If ibv=0 then use automatic selection method                          
  IF(ibv.eq.0)THEN                                                  
     maxsig=max(sig(1),sig(2))                                      
     minsig=min(sig(1),sig(2))                                      
     if(maxsig/minsig.le.200.d0)then                                
        ibvu=1                                                       
     else                                                           
        ibvu=2                                                       
     endif
  else
     ibvu=ibv! left as it was
  endif
  ibv=ibvu ! to output choice done when in auto mode
  if(inlu.eq.2)then                                                  
! 2-body aproximation for the central point, no derivatives   
     CALL alph_del2 (el,tobs,idsta,obs4,ider,dobde,TWOBO=.true.) 
     allin=obs4(1)
     delin=obs4(2)                                   
  endif
! ===================================================================== 
! compute ellipse in the elements space                                 
  CALL slinel(daddet,gameq,ceq,ceicel,b,v)                          
! ===========================================================           
! compute line of orbital elements                                      
  CALL linobs(ibvu,npo,el,axes,sig,b,v,sigma,ceicel,el_m,npo1)        
! ===========================================================           
  ng=0                                                              
  DO 7 n=1,npo1   
! chose method to handle nonlinearity                                   
     IF(inlu.eq.1)THEN                                                
! linear map from ellipse
        dal=DOT_PRODUCT(el_m(1:6,n),daddet(1:6,1))
        ddl=DOT_PRODUCT(el_m(1:6,n),daddet(1:6,2))
        al_m(n)=dal                                                    
        de_m(n)=ddl                                                    
! apparent magnitude is the one of the nominal orbit                    
        hmag_m(n)=hmagn  
! compute elements in plane corresponding to sky plane
        el_m(1:6,n)=el%coord+el_m(1:6,n)                   
     ELSEIF(inlu.eq.2)THEN  
! compute elements in plane corresponding to sky plane
        el_m(1:6,n)=el%coord+el_m(1:6,n) 
        elv=el
        elv%coord=el_m(1:6,n)
! 2-body propagation from ellipse, no derivatives      
        CALL alph_del2 (elv,tobs,idsta,obs4,ider,dobde,&
     &   pha,dis,rdot,dsun,elo,gallat,TWOBO=.true.) 
! compute apparent magnitude at time of observation                     
        hmag_m(n)=appmag(elv%h_mag,elv%g_mag,dsun,dis,pha)
! difference is with respect to 2-body approx., used w.r. to true orbit 
        al_m(n)=obs4(1)-allin                                            
        de_m(n)=obs4(2)-delin  
! should we store elongation etc????                                          
     ELSEIF(inlu.eq.3)THEN   
! compute elements in plane corresponding to sky plane
        el_m(1:6,n)=el%coord+el_m(1:6,n)
        elv=el
        elv%coord=el_m(1:6,n)                       
! full n-body propagation from ellipse, no derivatives   
        if(type.eq.'O')then  
           CALL alph_del2 (elv,tobs,idsta,obs4,ider,dobde, &
     &   pha,dis,rdot,dsun,elo,gallat) 
           al_m(n)=obs4(1)
           de_m(n)=obs4(2)                     
! other prediction data stored in common                                
           phav(n)=pha                                               
           disv(n)=dis                                               
           dsunv(n)=dsun                                             
           elov(n)=elo                                               
           gallav(n)=gallat                                          
           adotv(n)=obs4(3)                                             
           ddotv(n)=obs4(4)                                            
! compute apparent magnitude at time of observation                     
           hmag_m(n)=appmag(elv%h_mag,elv%g_mag,dsun,dis,pha)
        elseif(type.eq.'R'.or.type.eq.'V')then
! tech, posr,idstarad already set           
           CALL r_rdot (el,tobs,idstarad,tech,posr,posr,al_m(n),de_m(n), & 
                  &             dummy,dummy,0)  
           hmag_m(n)=0.d0                                             
        ELSE                                                         
           WRITE(*,*) 'predic_obs: type ',type,' not handled'  
           RETURN                         
        ENDIF
        al_m(n)=al_m(n)-alpha                                            
        de_m(n)=de_m(n)-delta                                            
     ELSE                                                            
        WRITE(*,*)' predic_obs: this we have not invented yet ', inl     
        RETURN                                                       
     ENDIF
! keep count of lost revolutions
     IF(type.eq.'O')THEN                                        
        IF(n.eq.1)THEN                                                  
           IF(al_m(n).gt.pig)al_m(n)=al_m(n)-dpig          
        ELSE                                                            
           CALL angupd(al_m(n),al_m(n-1),ng)                                
        ENDIF
     ENDIF
! temporary output
     if(type.eq.'O'.or.type.eq.'S')then         
        write(*,*)n,', RA/DEC (deg)',al_m(n)*degrad,de_m(n)*degrad,ng    
     elseif(type.eq.'R'.or.type.eq.'V')then        
        write(*,*)n,', R/RDOT (km,km/day)',al_m(n)*aukm,de_m(n)*aukm,ng
     endif
7 ENDDO
! ===================================================================== 
! ensure that LOV is consistent with nominal point                      
! first find midpoint of LOV, assume npo is even                        
  if(ibvu.eq.2)then                                                  
     nrev=nint((al_m(npo/2)+al_m(npo/2+1))/2.d0/dpig)                   
     write(*,*)'debug: nrev:',nrev                                  
     if(nrev.ne.0)then                                              
        do n=1,npo1                                                 
           al_m(n)=al_m(n)-nrev*dpig                                    
        enddo
     endif
  endif
END SUBROUTINE predic_obs
! ======================================================================
! PREDIC_OBS2
! generates attributable with its covariance
! ONLY to be used for ./src/triang directory
! =====================================================================
SUBROUTINE predic_obs2(el,idsta,tobs,att,uncert,rr,pha,dsun,twobo,dobde)
  USE station_coordinates 
  USE orbit_elements                 
  USE attributable
! ============= input ==================================================
  TYPE(orbit_elem), INTENT(IN) :: el ! elements
  DOUBLE PRECISION, INTENT(IN) :: tobs ! obs.time TDT
  INTEGER,INTENT(IN) ::  idsta ! station code   
! ============= output =================================================
  TYPE(attrib), INTENT(OUT) :: att ! best fit obs., apparent magnitude
! observations: alpha, delta in RAD  alphadot, deltadot in RAD/day
! ======optional input ================================
  TYPE(orb_uncert), INTENT(IN),OPTIONAL :: uncert ! covariance, normal matrices
! flag for 2-body approximation; must be .false. for full n-body computa
  LOGICAL,INTENT(IN),OPTIONAL :: twobo
! ======optional output =================================================
  DOUBLE PRECISION, INTENT(OUT),OPTIONAL :: rr(2) ! predicted r,rdot
  DOUBLE PRECISION, INTENT(OUT),OPTIONAL :: pha,dsun !phase, dist. to sun
! partial derivatives of obs, w.r. to asteroid elements (= dadde)
  DOUBLE PRECISION, INTENT(OUT), OPTIONAL:: dobde(4,6) 
! ============END INTERFACE=============================================
  DOUBLE PRECISION :: gameq(6,6) ! covariance of elements
  DOUBLE PRECISION :: att4(4), gamad(4,4) ! angles,covariance on att. 4-space
! partial derivatives of alpha, delta, w.r. to elements (by columns, by rows)    
  DOUBLE PRECISION daddet(6,4),dadde(4,6)
  DOUBLE PRECISION axes(4,4),sig(4) ! axes, eigenvalues (NOT to be used:polar bug)
  DOUBLE PRECISION :: pha0,dis0,dsun0,rdot0
! ===================================================================   
! functions                                                             
  DOUBLE PRECISION appmag
! elongation, galactic latitude (not really used)
  DOUBLE PRECISION elo,gallat,ecllat 
! TIME CONVERSION
  integer MJD1,MJD2
  double precision SEC1,SEC2
! ===================================================================   
! temporaries, indexes
  CHARACTER*3 obscod                                                  
  LOGICAL cov ! true if covariance computations are possible
  LOGICAL twobo1 
  INTEGER ider,i,j
!****************                                                       
!   static memory not required                                          
!****************                                                       
  IF(PRESENT(twobo))THEN
     twobo1=twobo
  ELSE
     twobo1=.false.
  ENDIF
! simple prediction without derivatives
  cov=PRESENT(uncert)
  IF(cov)cov=(cov.and.uncert%succ)
  IF(cov)THEN
    gameq=uncert%g
    ider=1
  ELSE
    ider=0
  ENDIF
! ===================================================================== 
! compute observation; derivatives (of order 1) if required                
  CALL alph_del2 (el,tobs,idsta,att4,ider,dadde,      &
            &   pha0,dis0,rdot0,dsun0,elo,gallat,ecllat,twobo1)
! optional arguments for magnitude
  IF(PRESENT(pha))pha=pha0
  IF(PRESENT(dsun))dsun=dsun0 
! also computation of r, rdot
  IF(PRESENT(rr))THEN
     rr(1)=dis0
     rr(2)=rdot0
  ENDIF
  att=undefined_attrib
  DO i=1,4
     att%angles(i)=att4(i)
  ENDDO
  att%tdtobs=tobs
! find UT of observation
  mjd1=FLOOR(att%tdtobs)
  sec1=(att%tdtobs-mjd1)*86400.d0
  CALL cnvtim(mjd1,sec1,'TDT',mjd2,sec2,'UTC')
  att%tutobs=mjd2+sec2/86400.d0
! station code
  CALL codestat(idsta,obscod)
  att%obscod=obscod
! compute apparent magnitude at time of observation                     
  att%apm=appmag(el%h_mag,el%g_mag,dsun0,dis0,pha0) 
  IF(.not.cov) RETURN
  daddet=TRANSPOSE(dadde)
  IF(PRESENT(dobde))THEN
     dobde=dadde
  ENDIF
! *******************fix for infamous polar bug*********************
! NOT DONE 
!  daddet(1:6,1)=daddet(1:6,1)*cos(att(2))
!  daddet(1:6,3)=daddet(1:6,3)*cos(att(2))
! *******************end fix polar bug******************************
! compute ellipsoid of covariance of alpha,delta,alphadot,deltadot
  CALL ellipsoid(daddet,gameq,sig,axes,gamad)
  att%g=gamad
END SUBROUTINE predic_obs2

! ===================================================================== 
! ALPH_DEL2
! ===================================================================== 
! Computation of alpha, delta, alphadot, deltadot  and their derivatives
! ===================================================================== 
!
! Input                                                                 
!    t0: epoch time                                                     
!    tobs: observation time                                             
!    iobscod: station code (integer) 
!    el: orbital elements 
!    ider: flag for derivatives options:                                
!        0 no derivatives                                               
!        1 only first deriv. of alpha,delta w.r. to coord
! Optional input      
!    twobo: logical flag for 2-body approximation; if .true., 2-body    
!        approximation (but Earth position comes from JPL ephem); if .fa
!        all orbit propagations are full n-body                         
! Output      
!   obs4: attributable, including
!           alpha,delta  computed at time tobs                        
!           adot,ddot their time derivatives
!   dobde: their partial derivatives with respect to 
!         the coordinates of el
! Optional output: pha0,dis0,rdot0,dsun0,elo0,gallat0
!
! ==============INTERFACE============================================   
SUBROUTINE alph_del2 (el,tobs,iobscod,obs4,ider,   &
     &   dobde,pha0,dis0,rdot0,dsun0,elo0,gallat0,eclat0,twobo,elev0,elsun0) 
  USE propag_state
  USE orbit_elements                 
! ==============INPUT==========================
  TYPE(orbit_elem), INTENT(IN) :: el ! elements
  double precision, intent(in) :: tobs ! observation time (MJD, TDT)
  integer,intent(in) :: iobscod ! observatory code
  integer, intent(IN) :: ider ! flag to control computation of derivatives
! ======optional input =================================================
  logical,intent(in), optional :: twobo ! two-body approx: if .false, full n-body  
! ============OUTPUT==============================                      
  double precision obs4(4) ! observations: alpha, delta in RAD  
! alphadot, deltadot in RAD/day
  double precision, intent(OUT), OPTIONAL :: dobde(4,6) ! partial derivatives 
! of obs, w.r. to asteroid elements
! ======optional output =================================================
! phase, distance to Earth, distance to Sun, elongation, galactic and ecliptic latitude
  DOUBLE PRECISION,INTENT(OUT),OPTIONAl :: pha0,dis0,dsun0,rdot0,elo0,gallat0,eclat0, elev0, elsun0
! =============END INTERFACE=========================================   
  LOGICAL twobo1 ! control of 2-body approximation
  double precision xast(6) ! asteroid cartesian coordinates 
  double precision xea(6) ! cartesian coordinates of the Earth 
  double precision dobdx(4,6) ! partial derivatives of alpha, delta, 
                        !  adot, ddot  w.r. elements
! derivatives of cart. coord. w. r. elements 
  double precision dxde(6,6),ddxde(3,6,6) 
! elongation,distance to Earth, distance to Sun (to compute magnitude)  
  DOUBLE PRECISION pha,dis,rdot,dsun,elo,gallat,eclat,elev,elsun
! **********************************************************************
!****************                                                       
!   static memory not required                                          
!****************             
! Orbit propagation:                                                    
  IF(PRESENT(twobo))then
     twobo1=twobo
  ELSE
     twobo1=.false. ! default is full n-body
  ENDIF
! propagation to time t2
  CALL propag(el,tobs,xast,xea,ider,dxde,twobo1) 
! Computation of observations                                           
  call oss_dif2(xast,xea,tobs,iobscod,obs4,ider,dobdx,      &
     &   pha,dis,rdot,dsun,elo,gallat,eclat,elev,elsun) 
 ! store true phase, etc.  
  IF(PRESENT(pha0))pha0=pha 
  IF(PRESENT(dis0))dis0=dis 
  IF(PRESENT(dsun0))dsun0=dsun
  IF(PRESENT(rdot0))rdot0=rdot 
  IF(PRESENT(elo0))elo0=elo 
  IF(PRESENT(gallat0))gallat0=gallat 
  IF(PRESENT(eclat0))eclat0=eclat
  IF(PRESENT(elev0))elev0=elev
  IF(PRESENT(elsun0))elsun0=elsun
  if(ider.lt.1) return 
! derivatives with respect to equinoctal elements                       
  IF(PRESENT(dobde))THEN
     dobde=MATMUL(dobdx,dxde)
  ENDIF 
END SUBROUTINE alph_del2
! ===================================================================== 
! OSS_DIF2                                                               
! ===================================================================== 
! Corrections to observations                                           
! ===================================================================== 
! Input                                                                 
! xast asteroid cartesian coordinates at time tauj                      
! xea Earth cartesian coordinates at time tauj                          
!     both in the ecliptic system                                       
! tauj observation time                                                 
! ioc station code                                                      
! ider flag for derivatives options:                                    
!       =0 no derivatives                                               
!       <2 only first deriv. of $\alpha,\delta$ w.r.                    
! to positions                                                          
!       great equal 2 also second partial derivatives                   
! (approximation with 2-body case)                                      
!                                                                       
! Output                                                                
! obs=alpha,delta,adot,ddot computed at time tauj (in the equatorial sys
! (if required)                                                         
! dobdx matrix of first derivatives w. r. to ecliptic positions         
! (if required)                                                         
! ===================================================================== 
SUBROUTINE oss_dif2(xast,xea,tobs,iobscod,obs4,ider,dobdx,      &
     &   pha0,dis0,rdot0,dsun0,elo0,gallat0,eclat0,elev,elsun) 
  USE reference_systems, ONLY: pvobs  
  USE force_model
! INPUT
  DOUBLE PRECISION, INTENT(IN) :: tobs ! observation time (MJD, TDT)
  DOUBLE PRECISION, INTENT(IN) :: xast(6),xea(6) ! ast and Earth cartesian coord
  INTEGER, INTENT(IN) :: iobscod,ider ! obs.code(integer), control on no. derivatives
! OUTPUT
  DOUBLE PRECISION, INTENT(OUT) :: obs4(4) ! observations alpha, delta, aodt,ddot
  DOUBLE PRECISION, INTENT(OUT) :: dobdx(4,6) ! partials of obs  w.r.t. cartesian ecliptic
! phase, distance to Earth, distance to Sun elongation, galactic latitude, elevation, elev. Sun
  DOUBLE PRECISION,INTENT(OUT) :: pha0,dis0,rdot0,dsun0,elo0,gallat0,eclat0, elev, elsun
! ===================================================================== 
  double precision d(6) ! vector difference of cartesian coordinates, ecliptic 
! WARNING: even if the velocities are not always propagated, they are available
! at the time of the observations                                       
!     note that alpha, delta are equatorial
  double precision xo(3),vo(3) ! topocentric position of the observatory 
  double precision xo2(3),vo2(3)
  double precision cospha,coselo ! phase and elongation cosine
  double precision sinelo,vvv(3)
  double precision dz ! auxiliary var
  double precision vsize,prscal ! real functions 
  double precision deq(6),tmp(3),ddd(3) ! rotation matr, rotated vectors
  double precision alpha,delta,adot,ddot ! scalar variables to preserve old code 
! partials of equatorial alpha, delta w.r. to cartesian ecliptic coord. 
  double precision dadx(3),dddx(3),ddadx(3,3),ddddx(3,3) 
! aux. var for computation of second derivatives                        
  double precision x,y,z 
  double precision den,x2,y2,z2,x4,y4 
  double precision x2y2,x2z2,y2z2 
  DOUBLE PRECISION coscoelev, coscoelsun
! ===================================================================== 
! Difference vector
  d=xast-xea  ! difference in ecliptic coordinates gives a geocentric vector
! ===================================================================== 
! Displacement of the station with respect to the center of the Earth   
  if(istat.gt.0.and.iobscod.ne.500)then 
     call pvobs(tobs,iobscod,xo,vo) 
     d(1:3)=d(1:3)-xo  !   topocentric correction in ecliptic coordinates
     d(4:6)=d(4:6)-vo  
! elevation on horizon 
     coscoelev=prscal(xo,d)/(vsize(d)*vsize(xo))
     IF(coscoelev.ge.0.d0.and.coscoelev.lt.1.d0)THEN
        elev=pig/2.d0-acos(coscoelev)
     ELSE
        elev=0.d0
!        WRITE(ierrou,*) 'oss_dif2: xo,d,dis0, coscoelev',xo,d,dis0,coscoelev
        WRITE(ierrou,*) 'oss_dif2: coscoelev',coscoelev
        numerr=numerr+1
     ENDIF
  endif
! ===================================================================== 
! Aberration (only time delay)                                          
  if(iaber.gt.0)then 
     call aber1(d,xast(4),d) 
  endif
! ===================================================================== 
! Computation of solar distance, earth distance, phase, elongation      
  dsun0=vsize(xast) 
  dis0=vsize(d) 
  rdot0=prscal(d(1:3),d(4:6))/dis0
  cospha=prscal(d,xast)/(dis0*dsun0) 
  pha0=acos(cospha) 
  coselo=-prscal(d,xea)/(dis0*vsize(xea)) 
  CALL prvec(d,xea,vvv)
  sinelo=-vvv(3)
  elo0=acos(coselo)
  IF(sinelo.lt.0.d0)elo0=-elo0
  eclat0=asin(d(3)/dis0)
! =====================================================================
! illumination angle at the station
  coscoelsun=prscal(d,-xea(1:3))/(dis0*dsun0)
  elsun=pig/2.d0-acos(coscoelsun)
! ===================================================================== 
! rotation to the equatorial reference system   
  deq(1:3)=MATMUL(roteceq,d(1:3)) 
  deq(4:6)=MATMUL(roteceq,d(4:6)) 
  d=deq                          
! ===================================================================== 
! galactic latitude                                                     
  gallat0=pig/2d0-acos((d(1)*gax+d(2)*gay+d(3)*gaz)/dis0) 
! ===================================================================== 
! Computation of observation: right ascension (radians)                 
  dz=d(1)**2+d(2)**2 
  if (dz.le.100.d0*epsilon(1.d0)) then
! remove singularity at poles 
     alpha=0.d0 
  else 
     alpha=atan2(d(2),d(1)) 
     if (alpha.lt.0.d0) then 
        alpha=alpha+dpig 
     endif
  endif
! Computation of observation: declination (radians)                     
  delta=asin(d(3)/dis0) 
! ===================================================================== 
! Computation of first derivatives of $\alpha$ and $\delta$ w.r. to posi
! (if required): we derived eq. (2.20)                                  
  dadx(1)=-d(2)/dz 
  dadx(2)=d(1)/dz 
  dadx(3)=0.d0 
  dddx(1)=-d(3)*(d(1)/(sqrt(dz)*dis0**2)) 
  dddx(2)=-d(3)*(d(2)/(sqrt(dz)*dis0**2)) 
  dddx(3)=sqrt(dz)/dis0**2 
! ===================================================================== 
! Apparent motion:                                                      
  adot=prscal(dadx,d(4)) 
  ddot=prscal(dddx,d(4)) 
! store into obs vector                                                 
  obs4(1)=alpha 
  obs4(2)=delta 
  obs4(3)=adot 
  obs4(4)=ddot 
! check if observation partials are required                            
  if(ider.eq.0)RETURN 
! ===================================================================== 
! partials of alpha, delta have already been computed                   
! partials of adot,ddot with respect to velocities are the same         
! rotation to the equatorial reference system                           
  tmp=MATMUL(roteqec,dadx) 
  dobdx(1,1:3)=tmp 
  dobdx(1,4:6)=0.d0 
  dobdx(3,4:6)=tmp 
  tmp=MATMUL(roteqec,dddx)
  dobdx(2,1:3)=tmp 
  dobdx(2,4:6)=0.d0 
  dobdx(4,4:6)=tmp
! ===================================================================== 
! partials of adot,ddot with respect to positions require the           
! second derivatives of alpha, delta with respect to the                
! equatorial reference system                                           
! ===================================================================== 
! Computation of second derivatives of $\alpha$ w.r. to positions       
  ddadx(1,1)=2.d0*d(1)*d(2)/dz**2 
  ddadx(1,2)=(d(2)**2-d(1)**2)/dz**2 
  ddadx(2,1)=ddadx(1,2) 
  ddadx(2,2)=-ddadx(1,1) 
  ddadx(3,1)=0.d0 
  ddadx(3,2)=0.d0 
  ddadx(3,3)=0.d0 
  ddadx(2,3)=0.d0 
  ddadx(1,3)=0.d0 
! Computation of second derivatives of $\delta$ w.r. to positions       
  den=1.d0/(dis0**4*dz*sqrt(dz)) 
  x=d(1) 
  y=d(2) 
  z=d(3) 
!                                                                       
  x2=x*x 
  y2=y*y 
  z2=z*z 
  x4=x2*x2 
  y4=y2*y2 
!                                                                       
  x2y2=x2*y2 
  x2z2=x2*z2 
  y2z2=y2*z2 
!                                                                       
  ddddx(1,1)=z*(2.d0*x4+x2y2-y2z2-y4)*den 
  ddddx(2,2)=z*(2.d0*y4+x2y2-x2z2-x4)*den 
  ddddx(1,2)=x*y*z*(z2+3.d0*x2+3.d0*y2)*den 
  ddddx(2,1)=ddddx(1,2) 
  ddddx(3,3)=-2.d0*z*dz**2*den 
  ddddx(1,3)=x*dz*(z2-x2-y2)*den 
  ddddx(3,1)=ddddx(1,3) 
  ddddx(2,3)=y*dz*(z2-x2-y2)*den 
  ddddx(3,2)=ddddx(2,3) 
! =======================================================               
! chain rule for derivatives of adot                                    
  ddd=MATMUL(ddadx,d(4:6))
! =======================================================               
! rotation to the equatorial reference system                           
  dobdx(3,1:3)=MATMUL(roteqec,ddd)
! =======================================================               
! chain rule for derivatives of adot                                    
  ddd=MATMUL(ddddx,d(4:6))
! =======================================================               
! rotation to the equatorial reference system                           
  dobdx(4,1:3)=MATMUL(roteqec,ddd)
END SUBROUTINE oss_dif2
! ===================================================================== 
! PRE_OBS_ATT   fortran90 version, A. Milani, May 2003
! ===================================================================== 
! method: call with minimum set of arguments for attribute only
!  input: minimum
!          el    = orbital elements, including
!                  epoch time (MJD,TDT)                                     
!          tobs  = prediction time (MJD,TDT)                                
!          eq    = equinoctial orbital elements vector at time t0
!  output: minimum
!          alpha = right ascension (equatorial J2000), radians          
!          delta = declination (equatorial J2000), radians
!          optional adot,ddot  proper motion            
! ============INTERFACE=================================================
SUBROUTINE pre_obs_att(el,tobs,alpha,delta,adot,ddot,twobo)
!  USE reference_systems 
  USE propag_state
  USE orbit_elements
! ============= input ==================================================
  TYPE(orbit_elem), INTENT(IN) :: el
  DOUBLE PRECISION, INTENT(IN) :: tobs ! obs.time
!  INTEGER, INTENT(IN) ::  idsta ! station code   
! ============= output =================================================
  DOUBLE PRECISION, INTENT(OUT) :: alpha,delta ! best fit obs.
  DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: adot,ddot ! proper motion
  LOGICAL, INTENT(IN), OPTIONAL :: twobo
! ============END INTERFACE============================================= 
! constant of gravitation, trigonometric constants from fund_const.mod 
  double precision xast(6) ! asteroid cartesian coordinates 
  double precision xea(6) ! cartesian coordinates of the Earth 
  double precision dxde(6,6) ! derivatives of cart. coord. w.r. to elements 
  INTEGER ider
  double precision dadx(3),dddx(3) ! partial of angles w.r. cart. coord
  double precision d(6),deq(6) ! diff. of cartesian coord., ecliptic, equatorial 
  double precision dis0,vsize,dz,prscal ! distance from earth, from z axis
  double precision rot(3,3),rotinv(3,3) ! rotation matr,inverse
  LOGICAL twobo1
! ===================================================================== 
! Orbit propagation:                                                    
  IF(PRESENT(twobo))then
     twobo1=twobo
  ELSE
     twobo1=.false. ! default is full n-body
  ENDIF
! compute observation; no derivatives
  ider=0
! full n-body numerical integration                                     
  call propag (el,tobs,xast,xea,ider,dxde,twobo1) 
! Difference vector
  d=xast-xea 
! rotation to the equatorial reference system   
  deq(1:3)=MATMUL(roteceq,d(1:3)) 
  deq(4:6)=MATMUL(roteceq,d(4:6))
  d=deq                         
! Computation of observation: right ascension (radians)                 
  dz=d(1)**2+d(2)**2 
  if (dz.le.100*epsilon(1.d0)) then
! remove singularity at poles 
     alpha=0.d0 
  else 
     alpha=atan2(d(2),d(1)) 
     if (alpha.lt.0.d0) then 
        alpha=alpha+dpig 
     endif
  endif
! Computation of observation: declination (radians) 
  dis0=vsize(d)
  delta=asin(d(3)/dis0) 
  IF(.not.(PRESENT(adot).or.PRESENT(ddot)))RETURN
! Computation of first derivatives of $\alpha$ and $\delta$ w.r. to posi
! (if required): we derived eq. (2.20)                                  
  dadx(1)=-d(2)/dz 
  dadx(2)=d(1)/dz 
  dadx(3)=0.d0 
  dddx(1)=-d(3)*(d(1)/(sqrt(dz)*dis0**2)) 
  dddx(2)=-d(3)*(d(2)/(sqrt(dz)*dis0**2)) 
  dddx(3)=sqrt(dz)/dis0**2 
! Apparent motion:
  IF(PRESENT(adot))THEN                                                      
     adot=prscal(dadx,d(4))
  ENDIF
  IF(PRESENT(ddot))THEN  
     ddot=prscal(dddx,d(4))
  ENDIF  
END SUBROUTINE pre_obs_att
! ===================================================================== 
! ALPH_DEL  vers. 3.0
! A. Milani, November 2002                    
! ===================================================================== 
! Computation of alpha and delta and their derivatives                  
! ===================================================================== 
!
! Input
!    el : orbital elements, including epoch time
!    tauj: observation time
!    ioj: station code; pos, vel position and velocity as computed in input    
!    ider: flag for derivatives options:                                
!        0 no derivatives  1 derivatives  
!    twobo: logical flag for 2-body approximation; if .true., 2-body    
!        approximation (but Earth position comes from JPL ephem); if .fa
!        all orbit propagations are full n-body                         
! Output
!   alj,dej: alpha,delta  computed at time tauj                         
!   dade,ddde:  matrices of first derivatives (if required)             
!
! ==============INTERFACE============================================   
SUBROUTINE alph_del (el,tauj,iocj,pos,vel,ider,twobo,alj,dej,dade,ddde, &
     &       adot0,ddot0,pha0,dis0,dsun0,elo0,gallat0)
  USE propag_state
  USE orbit_elements
  USE close_app, ONLY : kill_propag  
! kill_propag e' gia' nelle definizioni.....
!=========================================================
! INPUT
! flag to control computation of derivatives                            
  INTEGER, INTENT(IN) ::  ider
! flag to control two-body approximation: if .false., full n-body       
  LOGICAL, INTENT(IN) :: twobo 
! times: epoch time for asteroid elements, observation time (MJD)       
  DOUBLE PRECISION, INTENT(IN) :: tauj 
! asteroid elements                                          
  TYPE(orbit_elem), INTENT(IN) ::  el 
! observatory code and position, velocity
  INTEGER, INTENT(in) :: iocj 
  DOUBLE PRECISION, INTENT(IN) :: pos(3),vel(3)
! OUTPUT
! observations: alpha (right ascension) delta (eclination), in RAD      
  DOUBLE PRECISION, INTENT(OUT) ::  alj,dej
! partial derivatives of alpha, delta, w.r. to asteroid coordinates     
  DOUBLE PRECISION, INTENT(OUT) :: dade(6),ddde(6) 
! ======optional output =================================================
  DOUBLE PRECISION,INTENT(OUT),OPTIONAl :: adot0,ddot0 ! proper motion
  DOUBLE PRECISION,INTENT(OUT),OPTIONAl :: pha0,dis0,dsun0,elo0,gallat0 ! phase, 
             ! distance to Earth, distance to Sun, elongation, galactic latitude
! =============END INTERFACE=========================================   
  double precision xea(6),xobs(6) ! cartesian coordinates of the Earth, of the observer   
  double precision xast(6) ! asteroid cartesian coordinates
  double precision dadx(3),dddx(3) ! first derivatives of alpha, delta, w.r. to coord
  DOUBLE PRECISION :: adot,ddot,pha,dis,dsun,elo,gallat ! proper motion, data to compute magnitude
  double precision dxde(6,6) ! first derivatives of cartesian coordinates with respect to ele
  integer j ! loop variables j=1,6;
  double precision prscal ! double precision functions
! **********************************************************************
!****************                                                       
!   static memory not required                                          
!****************                                                       
! Orbit propagation:                                                    
! full n-body numerical integration                                     
  call propag (el,tauj,xast,xea,ider,dxde,twobo)
  xobs(1:3)=xea(1:3)+pos
  xobs(4:6)=xea(4:6)+vel
!  WRITE(11,899)tauj,xast,xobs
!899 FORMAT(f15.9,3(1x,f15.12),3(1x,f15.12),3(1x,f15.12),3(1x,f15.12))
!  IF(kill_propag) STOP '**** alph_del: kill propag *****'
  IF(kill_propag) RETURN
! Computation of observations                                           
  call oss_dif(xast,xea,tauj,iocj,pos,vel,alj,dej,ider,dadx,dddx,&
     &       adot,ddot,pha,dis,dsun,elo,gallat)
! store true phase, etc.
  IF(PRESENT(adot0))adot0=adot
  IF(PRESENT(ddot0))ddot0=ddot
  IF(PRESENT(pha0))pha0=pha 
  IF(PRESENT(dis0))dis0=dis 
  IF(PRESENT(dsun0))dsun0=dsun 
  IF(PRESENT(elo0))elo0=elo 
  IF(PRESENT(gallat0))gallat0=gallat 
  if(ider.lt.1)return                                               
! derivatives with respect to equinoctal elements                       
  DO  j=1,6
     dade(j)=prscal(dadx,dxde(1,j)) 
     ddde(j)=prscal(dddx,dxde(1,j)) 
  ENDDO
END SUBROUTINE alph_del
! ===================================================================== 
! OSS_DIF vers. 3.0 
! vers. 3.0 A. Milani, November 2002 (fortran 90)
! vers. 1.3 A. Milani, June 14, 1997 (equatorial coordinates for observations)
! ===================================================================== 
! Corrections to observations                                           
! ===================================================================== 
! Input                                                                 
! xast asteroid cartesian coordinates at time tauj                      
! xea Earth cartesian coordinates at time tauj                          
!     both in the ecliptic system                                       
! tauj observation time                                                 
! idst station code, pos, vel geocentric position and velocity
! ider flag for derivatives options:                                    
!       =0 no derivatives                                               
!       =1 first deriv. of $\alpha,\delta$ w.r. to positions 
! Output                                                                
! alj,dej,alpha,delta computed at time tauj (in the equatorial system)  
! (if required)                                                         
! dadx,dddx matrices of first derivatives w. r. to ecliptic positions   
! (if required)                                                         
! ddadx,ddddx matrices of second deriv., 2-b approximation              
! ===================================================================== 
 SUBROUTINE oss_dif(xast,xea,tauj,idst,pos,vel,alj,dej,ider,dadx,dddx, &
     &       adot,ddot,pha,dis,dsun,elo,gallat)
   USE fund_const
   USE force_model 
! ===================================================================== 
! INPUT
! control derivatives                        
   integer, intent(in) :: ider
! cartesian coordinates, ecliptic, of asteroid, of earth, 
! WARNING: even if the velocities are not always propagated, they are 
! available at the time of the observations                                       
   DOUBLE PRECISION, INTENT(IN) :: xast(6),xea(6)
! time tdtd of the observation
   DOUBLE PRECISION, INTENT(IN) :: tauj
! observatory code, geocentric position, velocity 
   INTEGER, INTENT(in) :: idst
   DOUBLE PRECISION, INTENT(IN) :: pos(3),vel(3)
! OUTPUT
! observations alpha, delta (equatorial, radians)              
   DOUBLE PRECISION, INTENT(OUT) :: alj,dej
! partials of equatorial alpha, delta w.r. to cartesian ecliptic coord. 
   DOUBLE PRECISION, INTENT(OUT) :: dadx(3),dddx(3)
! proper motion, data to compute magnitude
   DOUBLE PRECISION, INTENT(OUT) :: adot,ddot,pha,dis,dsun,elo,gallat
! END INTERFACE
! difference vector
   DOUBLE PRECISION d(6)
! topocentric position of the observatory ***** temporary for check
   double precision xo(3),vo(3) 
! phase and elongation cosine                                           
   double precision cospha,coselo,sinelo,vvv(3) 
! auxiliary var.                                          
   double precision dz 
! real functions                                                        
   double precision vsize,prscal 
! integer loop index                                                    
   integer i,ii,ij,j 
! rotation matrices from reference_systems.mod, rotated vectors
   double precision deq(6),tmp(3) 
! ===================================================================== 
! Difference vector 
   d=xast-xea ! difference in ecliptic coordinates gives a geocentric vector
! ===================================================================== 
! Displacement of the station with respect to the center of the Earth   
   if(istat.gt.0)then 
      IF(idst.ne.500)THEN 
! check *********************
!           call pvobs(tauj,idst,xo,vo) 
!           write(*,*)'pos -xo ', pos-xo
!           write(*,*)'vel -vo ', vel-vo
! end check *********************
!            call vdiff(d,xo,d) 
!            call vdiff(d(4),vo,d(4)) 
         d(1:3)=d(1:3)-pos
         d(4:6)=d(4:6)-vel
      ENDIF
   endif
! ===================================================================== 
! Aberration (only time delay)                                          
   if(iaber.gt.0)then 
      call aber1(d,xast(4),d) 
   endif
! ===================================================================== 
! Computation of solar distance, earth distance, phase, elongation      
   dsun=vsize(xast) 
   dis=vsize(d) 
   cospha=prscal(d,xast)/(dis*dsun) 
   pha=acos(cospha) 
   coselo=-prscal(d,xea)/(dis*vsize(xea)) 
   elo=acos(coselo) 
   CALL prvec(d,xea,vvv)
   sinelo=-vvv(3)
   If(sinelo.lt.0.d0)elo=-elo
! ===================================================================== 
! rotation to the equatorial reference system                           
   call prodmv(deq,roteceq,d) 
   call prodmv(deq(4),roteceq,d(4)) 
! trick to change as little as possible from vers. 1.2 to 1.3           
   d(1:6)=deq(1:6) 
! ===================================================================== 
! galactic latitude                                                     
   gallat=pig/2d0-acos((d(1)*gax+d(2)*gay+d(3)*gaz)/dis) 
! ===================================================================== 
! Computation of observation: right ascension (radians)                 
   dz=d(1)**2+d(2)**2 
   if (dz.le.100.d0*epsilon(1.d0)) then 
      alj=0.d0 
   else 
      alj=atan2(d(2),d(1)) 
      if (alj.lt.0.d0) then 
         alj=alj+dpig 
      endif
   endif
! Computation of observation: declination (radians)                     
   dej=asin(d(3)/dis) 
! ===================================================================== 
! Computation of first derivatives of $\alpha$ and $\delta$ w.r. to posi
! (if required): we derived eq. (2.20)                                  
   dadx(1)=-d(2)/dz 
   dadx(2)=d(1)/dz 
   dadx(3)=0.d0 
   dddx(1)=-d(3)*(d(1)/(sqrt(dz)*dis**2)) 
   dddx(2)=-d(3)*(d(2)/(sqrt(dz)*dis**2)) 
   dddx(3)=sqrt(dz)/dis**2 
! ===================================================================== 
! Apparent motion:                                                      
   adot=prscal(dadx,d(4)) 
   ddot=prscal(dddx,d(4)) 
   if(ider.eq.0)return 
! ===================================================================== 
! rotation to the equatorial reference system 
   dadx=matmul(roteqec,dadx)
   dddx=matmul(roteqec,dddx) 
 END SUBROUTINE oss_dif
! ==================================
!
! R _ R D O T
!
! radar observations
SUBROUTINE r_rdot (el,tr,ioc,tech,posr,post,r,v,drde,dvde,ider) 
  USE reference_systems 
  USE propag_state
  USE astrometric_observations, ONLY: radius ! to set default value
  USE orbit_elements
  USE fund_const
  USE ever_pitkin
! =============INPUT====================                                
  TYPE(orbit_elem),intent(IN) :: el ! asteroid equinoctal elements 
! observation technology (used for surface bounce correction)                 
  CHARACTER*1, INTENT(IN):: tech 
! observation time (MJD)       
  double precision, intent(in) :: tr 
! observatory code (integer) with encoded transmitter and receiver
  integer, intent(in) :: ioc 
! flag to control computation of derivatives                            
  integer, intent(in):: ider 
  DOUBLE PRECISION, INTENT(IN) :: posr(3), post(3) 
! normally position and velocity of observer
! but in this special case, body fixed geocentric position 
! of receiver and of transmitterer
! ============OUTPUT====================                                
! observations: range and range rate in AU, AU/day                      
  double precision, intent(out) :: r,v 
! partial derivatives of range and range w.r. to asteroid coordinates   
  double precision, intent(out) ::  drde(6),dvde(6) 
! =============END INTERFACE=========================================   
! cartesian coordinates of the Earth, of the asteroid, id. at receive time; 
! barycentric cood. of Sun
  double precision xea(6),xast(6),xastr(6),xsun(6) 
! first partial derivatives of r, rdot, w.r. to ast. coordinates, vel.  
  double precision drdx(6),dvdx(6) 
! asteroid radius (for surface bounce correction)                       
  double precision rb 
! station codes                                                         
  INTEGER iotr,iore 
! station positions: geocentric, heliocentric                           
  double precision xre(3),yre(3),xtr(3),ytr(3)
  double precision xre1(3),yre1(3),xtr1(3),ytr1(3) 
  double precision rre(3),vre(3),rtr(3),vtr(3) 
! difference vector, distance, velocity difference, size                
  double precision rhorv(3),rhor,rhordv(3),rhord 
  double precision rhotv(3),rhot,rhotdv(3),rhotd 
  double precision vsize,prscal !functions 
! solar system barycentric velocities (for 1/c^2 Doppler level)         
  double precision vressb(3),vastssb(3),vtrssb(3) 
  double precision vressb2,vtrssb2 
! down leg time, bounce time                                            
  double precision taud,taudold,tb,tbold 
! up leg time, transmit time                                            
  double precision tauu,tt,ttold,tauuold 
! delta tau correction function is part of module 
! speed of light from fund-const.mod
! radius of asteroid from astrometric_observations.mod 
! iteration index, control and flag for the position/velocity subroutine 
  INTEGER i,ifla 
  INTEGER, PARAMETER :: itmax=10 
! control on convergence set at 0.05 microseconds (Yeomans et al. AJ 199
  DOUBLE PRECISION ept 
!     PARAMETER (ept=6.d-13)                                            
! correction 9/1/2002: control on taud, not on time in MJD              
  PARAMETER (ept=1.d-16) 
  DOUBLE PRECISION tbdif(itmax),taudif(itmax) 
! time scale                                                            
  double precision tretdb,temp,tdiffr,tdifft 
! 2-body step                                                           
  DOUBLE PRECISION eqast(6),enne 
! first derivatives of cartesian coordinates with respect to elements   
  double precision dxde(6,6) 
! auxiliary scalars for Doppler GR & tropospheric corrections           
  double precision rgeo,rgeod,deldoptr,deldopre,rast,rastd 
  double precision rsta,rstad,scal1,scal2,scal3,levelc 
  double precision dtaud,dtauu
! =================================================                     
! deal with surface vs. mass center return:                             
  if(tech.eq.'s')then 
     rb=radius 
     if(rb.le.0)then 
        write(*,*)  '**** rrdot: internal error', radius 
        stop 
     endif
  elseif(tech.eq.'c')then 
! no correction is needed                                               
     rb=0 
  else
     WRITE(*,*)' r_rdot; error in bs%tech', tech,' at time ', tr 
     WRITE(ierrou,*)' r_rdot; error in bs%tech', tech,' at time ', tr 
     numerr=numerr+1
     rb=0        
  endif
! ======================================                                
! Displacement of the stations with respect to the center of the Earth  
! find codes of two observatories                                       
  iotr=ioc/10000 
  iore=ioc-iotr*10000 
  ifla=1 
! compute position of receiver at time tr                               
  call pvobs3(tr,posr,xre,yre)
!     WRITE(*,*)' rec pos ',xre-xre1
!     WRITE(*,*)' rec vel ',yre-yre1
! =======================================                               
! get receive time in TDB                                               
  call times(tr+2400000.5,temp,tdiffr) 
  tretdb=tr+tdiffr/86400.d0 
! Compute down leg time                                                 
! Orbit propagation at time tr                                          
  call propag(el,tretdb,xastr,xea,ider,dxde) 
! ... new initial conditions are xastr                                  
! receive station position at time tr                                   
  rre=xea(1:3) + xre  ! CALL vsumg(3,xea,xre,rre) 
  vre=xea(4:6) + yre  ! CALL vsumg(3,xea(4),yre,vre) 
! initial guess is bounce time(position) equal to receive time(position)
  tb=tretdb 
  taud=0.d0 
  xast=xastr 
 ! Loop on down leg                                                      
  DO i=1,itmax 
     rhorv=xast(1:3)-rre ! CALL vdiff(xast,rre,rhorv) 
     rhor=vsize(rhorv) 
     dtaud=deltau(xast,xre,rhorv,rre) 
! correction                                                            
     taudold=taud 
     taud=(rhor-rb)/vlight + dtaud 
     tbold=tb 
     tb=tretdb-taud 
! convergence control                                                   
     tbdif(i)=tb-tbold 
     taudif(i)=taud-taudold 
!         WRITE(*,*)tb,tb-tbold,taud,taud-taudold                       
     IF(abs(taud-taudold).lt.ept) GOTO 9 
! Orbit propagation at time tb    
!     CALL coocha(xastr,'CAR',gms,eqast,'EQU',enne) 
     CALL fser_propag(xastr(1:3),xastr(4:6),tretdb,tb,gms,xast(1:3),xast(4:6))
!     CALL prop2b(tretdb,eqast,tb,xast,gms,0,dxde,ddxde) 
  ENDDO
! too many iterations                                                   
  WRITE(*,*)' slow conv. on down leg time ',(tbdif(i),i=1,itmax),   &
     &  (taudif(i),i=1,itmax)                                           
! compute relative velocity between asteroid and receiver               
9 CONTINUE 
  rhordv=xast(4:6)-vre  ! CALL vdiff(xast(4),vre,rhordv) 
  rhord=prscal(rhorv,rhordv)/rhor 
! solar velocity at the receive and bounce epochs to get the asteroid   
! barycentric velocity (vastssb) and the receiver barycentric velocity  
! (vressb)                                                              
  ifla=2 
! - receive                                                             
  CALL earcar(tretdb,xsun,ifla) 
  vressb=xsun(4:6)+vre ! CALL vsumg(3,vre,xsun(4),vressb) 
  vressb2=vressb(1)*vressb(1)+vressb(2)*vressb(2)+vressb(3)*vressb(3)
! - bounce                                                              
  CALL earcar(tb,xsun,ifla) 
  vastssb=xast(4:6)+xsun(4:6)  ! CALL vsumg(3,xast(4),xsun(4),vastssb) 
! =======================================                               
! compute upleg time                                                    
! fist guess is up leg time = down leg time                             
  ifla=1 
  tt=tb-taud 
  tauu=taud 
! Loop on upleg time                                                    
  DO i=1,itmax 
! compute transmitter position at estimated transmit time tt            
! call the Earth-rotation model in TDT                                  
     call times(tt+2400000.5,temp,tdifft) 
     CALL earcar(tt-(tdifft/86400.d0),xea,ifla) 
     CALL pvobs3(tt,post,xtr,ytr)
!        WRITE(*,*) ' tra pos ',xtr-xtr1
!        WRITE(*,*) ' tra vel ',ytr-ytr1
     rtr=xea(1:3) + xtr   !         CALL vsumg(3,xea,xtr,rtr) 
! upleg time                                                            
     rhotv=xast(1:3)-rtr ! CALL vdiff(xast,rtr,rhotv) 
     rhot=vsize(rhotv) 
     dtauu=deltau(xast,xtr,rhotv,rtr) 
     tauuold=tauu 
     tauu=(rhot-rb)/vlight + dtauu 
     ttold=tt 
     tt=tb-tauu 
! convergence control                                                   
!        WRITE(*,*)tt,tt-ttold,tauu, tauu-tauuold                       
     IF(abs(tauu-tauuold).lt.ept) GOTO 19 
  ENDDO
! too many iterations                                                   
  WRITE(*,*)' slow convergence up leg time ',tt-ttold, tauu-tauuold 
19 CONTINUE 
!      WRITE(*,*)' tauu at convergence ', tauu
!      WRITE(*,*)' rtr ', rtr
! compute relative velocity between asteroid and transmitter            
  vtr=xea(4:6)+ytr    ! CALL vsumg(3,xea(4),ytr,vtr) 
  rhotdv=xast(4:6)-vtr   !  CALL vdiff(xast(4),vtr,rhotdv) 
  rhotd=prscal(rhotv,rhotdv)/rhot 
! ==========================================================            
! compute distance                                                      
  r=0.5d0*(tauu+taud+(tdifft-tdiffr)/86400.d0)*vlight 
! compute relative frequency shift (up to 1/c^3 level);                 
! solar velocity at the bounce epochs and the asteroid barycentric      
! velocity (vastssb)                                                    
  ifla=2 
  CALL earcar(tt,xsun,ifla) 
  vtrssb=xsun(4:6) + vtr  ! CALL vsumg(3,vtr,xsun(4),vtrssb) 
  vtrssb2=vtrssb(1)*vtrssb(1)+vtrssb(2)*vtrssb(2)+ vtrssb(3)*vtrssb(3)
  levelc=rhotd+rhord 
  scal1=prscal(rhotv,vtrssb)/rhot/vlight 
  scal2=prscal(rhorv,vastssb)/rhor/vlight 
  scal3=0.5d0*(vtrssb2-vressb2)+                                    &
     &      gms*((1.d0/vsize(rtr))-(1.d0/vsize(rre)))                   
  v=0.5d0*(levelc+rhotd*scal1*(1.d0+scal1)                          &
     &               -rhord*scal2*(1.d0-scal2)                          &
     &               -(rhotd*rhord*(1.d0+scal1-scal2)                   &
     &               -scal3*(1.d0-(levelc/vlight)))/vlight)             
! - get GR and tropospheric corrections to Doppler                      
! -- GR stuff                                                           
  rast=vsize(xast) 
  rastd=prscal(xast,xast(4))/rast 
! a) upleg piece                                                        
  rsta=vsize(rtr) 
  rstad=prscal(rtr,vtr)/rsta 
  call deldop1(rast,rastd,rhot,rhotd,rsta,rstad,deldoptr) 
! b) downleg piece                                                      
  rsta=vsize(rre) 
  rstad=prscal(rre,vre)/rsta 
  call deldop1(rast,rastd,rhor,rhord,rsta,rstad,deldopre) 
  v=v+0.5d0*(deldoptr+deldopre) 
! -- troposheric stuff                                                  
! a) at transmit passage -->                                            
  rgeo=vsize(xtr) 
  rgeod=prscal(xtr,ytr)/rgeo 
  call deldop2(xast,xtr,ytr,rgeo,rgeod,rhotv,rhotdv,rhot,rhotd,deldoptr)
! b) at receive passage <--                                             
  rgeo=vsize(xre) 
  rgeod=prscal(xre,yre)/rgeo 
  call deldop2(xast,xre,yre,rgeo,rgeod,rhorv,rhordv,rhor,rhord,deldopre)
  v=v+0.5d0*(deldoptr+deldopre) 
! rem interplanetary environment effects (e^- plasma) neglected         
! ==========================================================            
! Derivatives                                                           
  IF(ider.eq.0)RETURN 
! derivs of r,rdot wrt cartesian                                        
  do i=1,3 
! d(range)/d(r)                                                         
     drdx(i)=(rhotv(i)/rhot +rhorv(i)/rhor)/2.d0 
! d(range)/d(v) = 0                                                     
     drdx(i+3)=0 
! d(range-rate)/d(r)                                                    
     dvdx(i)=-((rhotd*rhotv(i)/rhot-rhotdv(i))/rhot +               &
     &             (rhord*rhorv(i)/rhor-rhordv(i))/rhor)/2.d0           
! d(range-rate)/d(v)= c * d(range)/d(r)                                 
     dvdx(i+3)=drdx(i) 
  enddo
! derivs of cartesian with respect to elements               
  do i=1,6 
     drde(i)=DOT_PRODUCT(drdx,dxde(1:6,i))   ! drde(i)=prscag(6,drdx,dxde(1,i)) 
     dvde(i)=DOT_PRODUCT(dvdx,dxde(1:6,i))   ! dvde(i)=prscag(6,dvdx,dxde(1,i)) 
  enddo
  RETURN 
CONTAINS
! ======================================================================
! DELTAU - "small" corrections to radar time of flight                  
! ======================================================================
!     xast - asteroid position, heliocentric                            
!     xsta - station position, relative to Earth center                 
!     rho - asteroid position, relative to station                      
!     r - station position, heliocentric                                
  DOUBLE PRECISION FUNCTION deltau(xast,xsta,rho,r) 
    double precision xast(3),xsta(3),rho(3),r(3) 
    double precision vsize 
    double precision rsta,e,p,q,cosz,cotz,deltau1,deltau2 
!     double precision sinha,sin2ha,fghz,ampli,finte1,fun1,fun2         
!     double precision aprim,bprim,deltau3,omeg1,alpha                  
! ================================                                      
! Relativistic delay                                                    
    e=vsize(r) 
    p=vsize(xast) 
    q=vsize(rho) 
    deltau1=2d0*gms*log(abs((e+p+q)/(e+p-q)))/vlight**3 
! Earth ionospheric/tropospheric delay                                  
! ref. EM Standish, A&A 233, 252 (1990)                                 
    rsta=vsize(xsta) 
    if(rsta.lt.1d-12)then 
       write(*,*)'deltau: radar station at geocenter!' 
       cosz=1d0-1d-8 
    else 
       cosz=prscal(xsta,rho)/rsta/q ! cosz=prscag(3,xsta,rho)/rsta/q 
    endif
    if(cosz.eq.1d0)cosz=cosz-1d-8 
    cotz=cosz/sqrt(1d0-cosz**2) 
    deltau2=(7d-9/86400d0)/(cosz+1.4d-3/(4.5d-2+cotz)) 
! Interplanetary medium propagation effect                              
! ref. EM Standish, A&A 233, 252 (1990) with constants of DE118         
! rem. a more precise model might be needed here; e.g.                  
!      Muhleman & Anderson, ApJ 247, 1093 (1981) or newer               
!      alpha=q/e                                                        
!      cosha=(r(1)*rho(1)+r(2)*rho(2)+r(3)*rho(3))/e/q                  
!      sin2ha=1.d0-cosha*cosha                                          
!      sinha=dsqrt(sin2ha)                                              
! X- or S-band;                                                         
!c !!! Information about the frequency is not passed here at the        
!c     moment; it should be decided manually !!!                        
!c      fghz=2.38d0                                                     
!      fghz=8.51d0                                                      
!      ampli=(2.01094d-8)/fghz/fghz/e/86400.d0                          
!      aprim=1.237265d-6                                                
!      bprim=9.524021d0                                                 
!      finte1=(datan((alpha+cosha)/sinha)-datan(cosha/sinha))/sinha     
!      fun1=bprim*finte1                                                
!      omeg1=1.d0+alpha*(2.d0*cosh+alpha)                               
!      fun2=aprim*(0.25d0*((alpha+cosha)*((1.d0/omeg1)+(1.5d0/sin2ha))  
!     .           /omeg1-cosha*(1.d0+(1.5d0/sin2ha)))/sin2h             
!     .           +0.375d0*finte1/sin2ha/sin2ha)/(e**4)                 
!      deltau3=ampli*(fun1+fun2)                                        
! Add 'em up                                                            
!c      deltau=deltau1+deltau2+deltau3                                  
    deltau=deltau1+deltau2 
    return 
  END FUNCTION deltau
! ======================================================================
! DELDOP1 - "small" corrections to radar-rate measurements              
! ======================================================================
  SUBROUTINE deldop1(p,pdot,q,qdot,e,edot,deldop) 
    double precision p,pdot,q,qdot,e,edot,deldop 
    double precision brac1,brac2 
! ================================                                      
! relativistic range-rate correction                                    
    brac1=-q*(edot+pdot)-qdot*(e+p) 
    brac2=((e+p)**2)-(q**2) 
    deldop=4.d0*gms*brac1/brac2/(vlight**2) 
  END SUBROUTINE deldop1
! ======================================================================
! DELDOP2 - "small" corrections to radar-rate measurements              
! ======================================================================
!     xast(6) - asteroid r & v heliocentric                             
!     xtr(3),ytr(3) - station r & v relative to Earth center            
!     rsta,drsta - |xtr| & d|xtr|/dt                                    
!     rhov(3),drhov(3) - asteroid r & v relative to station             
!     rho,drho - |rhov| & d|rhov|/dt                                    
  SUBROUTINE deldop2(xast,xsta,vsta,rsta,drsta,rhov,drhov,rho,drho,deldop)
    double precision xast(6),xsta(3),vsta(3),rhov(3),drhov(3) 
    double precision rsta,drsta,rho,drho,deldop 
    double precision cosz,sinz,cotz,dcoszdt,phiz 
    double precision brac,brac1,brac2,scal1,scal2 
! ================================                                      
! rate of change of the Earth ionospheric/tropospheric ~ Doppler shift  
    if(rsta.lt.1d-12)then 
       write(*,*)'deltau: radar station at geocenter!' 
       cosz=1d0-1d-8 
    else 
       cosz=(rhov(1)*xsta(1)+rhov(2)*xsta(2)+rhov(3)*xsta(3))/rsta/rho 
    endif
    sinz=sqrt(1.d0-cosz**2) 
    cotz=cosz/sinz 
    brac=0.045d0+cotz 
    brac1=1.d0-(0.0014d0/brac/brac/(sinz**3)) 
    brac2=cosz+(0.0014d0/brac) 
    phiz=-vlight*(7.d-9/86400.d0)*brac1/(brac2**2) 
    scal1=drhov(1)*xsta(1)+drhov(2)*xsta(2)+drhov(3)*xsta(3) 
    scal2=rhov(1)*vsta(1)+rhov(2)*vsta(2)+rhov(3)*vsta(3) 
    dcoszdt=((scal1+scal2)/rho/rsta)-cosz*((drho/rho)+(drsta/rsta)) 
    deldop=phiz*dcoszdt 
  END SUBROUTINE deldop2

END SUBROUTINE r_rdot

END MODULE pred_obs  

! Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: February 24, 1997                                            
!                                                                       
!  ***************************************************************      
!  *                                                             *      
!  *                          A B E R 1                          *      
!  *                                                             *      
!  *          Correzione approssimata per aberrazione            *      
!  *                  stellare e/o planetaria                    *      
!  *                                                             *      
!  ***************************************************************      
!                                                                       
!                                                                       
! INPUT:    XREL(3)   -  Posizione relativa vera del corpo osservato (UA
!           VREL(3)   -  Velocita` relativa (UA/d)                      
!                                                                       
! OUTPUT:   XCOR(3)   -  Posizione relativa apparente (tenendo conto del
!                        aberrazione)                                   
!                                                                       
! NOTA: in generale per ottenere la correzione completa (comprendente le
!       cosiddette aberrazioni "stellare" + "planetaria") bisogna che VR
!       sia la velocita` relativa del corpo osservato rispetto all'osser
!       tore:                                                           
!                 VREL  =   V(pianeta) - V(osservatore)                 
!                                                                       
!       Se si vuole ottenere solo la correzione per l'aberrazione "stell
!       bisogna porre VREL =  - V(osservatore).                         
! 
SUBROUTINE aber1(xrel,vrel,xcor) 
  USE fund_const
  IMPLICIT NONE
  DOUBLE PRECISION :: xrel(3),vrel(3)
  DOUBLE PRECISION :: xcor(3) ! removed INTENT to allow xrel in the same location as xcor
  DOUBLE PRECISION ro,vsize,dt
!  distanza
  ro=vsize(xrel) 
!  effetto di ritardo                                                   
  dt=ro/vlight
  xcor=xrel-dt*vrel 
END SUBROUTINE aber1
! ===================================================================== 
! OUTOBC                                                                
! ===================================================================== 
!  output of predicted observation, possibly with confidence ellipse    
!   input: iun   = output unit                                          
!          type  = observation type                                     
!          ids = station code                                           
!          t1 = time of observation (UTC)                               
!          alpha, delta, hmagn = observation                            
!          adot,ddot = proper motion                                    
!          elo,dis = elongation, distance from Earth                    
!                                                                       
!          icov  = 1 for observations only, 2 to add confidence ellipse 
!          gamad,sig,axes = covariance matrix, sigmas along axes        
!                      (only for icov=2, otherwise dummy)               
! ===================================================================== 
      SUBROUTINE outobc(iun,type,ids,t1,alpha,delta,hmagn,adot,ddot,    &
     &     elo,dis,icov,gamad,sig,axes)  
      USE fund_const                               
      implicit none 
! needs AU value in km: from fund_cons                                                  
! output unit, station code, obs. type                                  
      integer iun,ids
      CHARACTER*(1) type
! observations                                                          
      double precision t1,alpha,delta,hmagn,adot,ddot,elo,dis 
! covariance                                                            
      integer icov 
      double precision gamad(2,2),axes(2,2),sig(2) 
! ================end interface===============================          
      double precision princ 
      integer i,j 
! time variables                                                        
      integer ideg,iday,imonth,iyear,ihour,imin,isec,ln,truncat 
      double precision hour,minu,sec 
      CHARACTER*22 timstr 
      CHARACTER*19 tmpstr 
      CHARACTER*12 rastri,rdstri 
      CHARACTER*1 signo 
! convert time                                                          
      CALL mjddat(t1,iday,imonth,iyear,hour) 
! convert hour to 12:12:12                                              
      ihour=truncat(hour,1d-7) 
      minu=(hour-ihour)*60.d0 
      imin=truncat(minu,1d-5) 
      sec=(minu-imin)*60.d0 
      isec=truncat(sec,1d-3) 
      WRITE(timstr,192) iyear,imonth,iday,ihour,imin,isec,sec-isec 
  192 FORMAT(I4,'/',I2.2,'/',I2.2,1x,I2.2,':',I2.2,':',I2.2,f3.2) 
! =================== select by observation type ===================    
      IF(type.eq.'O'.or.type.eq.'S')THEN 
! %%%%%%%%%%%% ASTROMETRY %%%%%%%%%%%%%%%%                              
! convert RA                                                            
         alpha=princ(alpha) 
         CALL sessag(alpha*degrad/15.d0,signo,ihour,imin,sec) 
         IF(signo.eq.'-')STOP 'wrirms error: negative right ascension.' 
! prepare RA string                                                     
         WRITE(tmpstr,FMT='(F6.3)') sec 
         CALL rmsp(tmpstr,ln) 
         IF(ln.lt.6)tmpstr='0'//tmpstr 
         WRITE(rastri,130)ihour,imin,tmpstr 
  130    FORMAT(I2.2,':',I2.2,':',a6) 
! convert DEC                                                           
         CALL sessag(delta*degrad,signo,ideg,imin,sec) 
! prepare DEC string                                                    
         WRITE(tmpstr,FMT='(F5.2)') sec 
         CALL rmsp(tmpstr,ln) 
         IF(ln.lt.5)tmpstr='0'//tmpstr 
         WRITE(rdstri,170)signo,ideg,imin,tmpstr 
  170    FORMAT(A1,I2.2,1x,I2.2,1x,a5) 
                                                                        
         write(iun,101)timstr,t1,ids,                                   &
     &        rastri,alpha*degrad,                                      &
     &        rdstri,delta*degrad,                                      &
     &        secrad*adot/24.d0,secrad*ddot/24.d0,                      &
     &        dis,elo*degrad,hmagn                                      
         write(*,101)timstr,t1,ids,                                     &
     &        rastri,alpha*degrad,                                      &
     &        rdstri,delta*degrad,                                      &
     &        secrad*adot/24.d0,secrad*ddot/24.d0,                      &
     &        dis,elo*degrad,hmagn                                      
  101    format('Astrometric Observation Prediction'/                   &
     &        'For ',a19,' (UTC); ',f12.5,'(MJD)'/                      &
     &        'Observatory code= ',i4.4/                                &
     &        'RA= ',a12,' (HH:MM:SS); ',f11.5,' (deg)'/                &
     &        'DEC= ',a12,' (deg min sec); ',f11.5,' (deg)'/            &
     &        'RA/DEC Apparent motion=',2(2x,f9.2),' (arcsec/hour)'/    &
     &        'Earth distance= ',f8.4,' (AU)'/                          &
     &        'Solar elongation= ',f7.2,' (deg)'/                       &
     &        'Apparent magnitude= ',f5.2)                              
         IF(icov.eq.1)RETURN 
! rescaling in arcsec                                                   
         do  i=1,2 
            sig(i)=sig(i)*secrad 
            do  j=1,2 
               gamad(i,j)=gamad(i,j)*secrad**2 
            enddo 
         enddo 
         write(iun,201)(sig(j),(axes(i,j),i=1,2),j=1,2) 
         write(*,201)(sig(j),(axes(i,j),i=1,2),j=1,2) 
  201    format(                                                        &
     &'Size and orientation of 1-sigma uncertainty ellipse'/            &
     &'Short axis : Size= ',1p,g12.6 ,' (arcsec); Direction= ',         &
     & 0p,2(1x,f8.5)/                                                   &
     &'Long axis : Size= ',1p,g12.6 ,' (arcsec); Direction= ',          &
     & 0p,2(1x,f8.5))                                                   
      ELSEIF(type.eq.'R'.or.type.eq.'V')THEN 
! %%%%%%%%%%%% RADAR %%%%%%%%%%%%%%%%                                   
         write(iun,102)t1,ids,alpha*aukm,delta*aukm 
         write(*,102)t1,ids,alpha*aukm,delta*aukm 
  102    format('time, MJD=',f13.6,'  station=',i4/                     &
     &       ' range (KM)         = ',f16.5/                            &
     &       ' range rate (KM/DAY)=  ',f15.5)                           
         IF(icov.eq.1)RETURN 
! rescaling in km, km/day                                               
         do  i=1,2 
            sig(i)=sig(i)*aukm 
            do  j=1,2 
               gamad(i,j)=gamad(i,j)*aukm**2 
            enddo 
         enddo 
         write(iun,202)(sig(j),(axes(i,j),i=1,2),j=1,2) 
         write(*,202)(sig(j),(axes(i,j),i=1,2),j=1,2) 
  202    format(' in the range (KM), range-rate (KM/DAY) plane'/        &
     &          ' sigma1 = ',1p,g14.7 ,' axis1= ',2(1x,g12.5)/          &
     &          ' sigma2 = ',1p,g14.7 ,' axis2= ',2(1x,g12.5))          
      ELSEIF(type.eq.'P')THEN 
! %%%%%%%%%%%% PROPER MOTION %%%%%%%%%%%%%%%%                           
         write(iun,109)t1,ids,alpha*secrad/24.d0,delta*secrad/24.d0 
         write(*,109)t1,ids,alpha*secrad/24.d0,delta*secrad/24.d0 
  109    format('time, MJD=',f13.6,'  station=',i4/                     &
     &       ' RA motion (arcsec/hour)     = ',f9.2/                    &
     &       ' DEC motion (arcsec/hour)    = ',f9.2)                    
         IF(icov.eq.1)RETURN 
! rescaling in arcsec/hour                                              
         do  i=1,2 
            sig(i)=sig(i)*secrad/24.d0 
            do  j=1,2 
               gamad(i,j)=gamad(i,j)*(secrad/24.d0)**2 
            enddo 
         enddo 
         write(iun,209)(sig(j),(axes(i,j),i=1,2),j=1,2) 
         write(*,209)(sig(j),(axes(i,j),i=1,2),j=1,2) 
  209    format('sigma1 (arcsec/hr)= ',1p,g14.7 ,' axis1= ',2(1x,g12.5)/&
     &          'sigma2 (arcsec/hr)= ',1p,g14.7 ,' axis2= ',2(1x,g12.5))
      ELSE 
         WRITE(*,*)'outobc: type=',type,' not understood' 
      ENDIF 
      return 
      END SUBROUTINE outobc



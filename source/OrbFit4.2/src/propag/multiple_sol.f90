! ========MODULE multiple_sol================================
!
MODULE multiple_sol  
USE output_control
USE fund_const
USE astrometric_observations
USE orbit_elements
USE least_squares
IMPLICIT NONE
PRIVATE

! LIST OF PUBLIC ENTITIES
! PUBLIC ROUTINEs
PUBLIC :: f_multi, nmulti, mult_input
PUBLIC :: fmuobs, fmupro, fmuplo, outmul, prop_sig, step_fit2
! former lov_int.mod
PUBLIC :: lovinit,lovobs,lovmagn,lovinterp
! PUBLIC data
! maximum number of multiple solutions
  INTEGER, PARAMETER::  mulx=9999
! first, last, reference solution                            
  INTEGER imip,imim,imi0
! multiple solution arrays                                                
  TYPE(orbit_elem), DIMENSION(mulx) :: elm
  DOUBLE PRECISION csinom(mulx),delnom(mulx)
  TYPE(orb_uncert), DIMENSION(mulx) :: unm
  DOUBLE PRECISION tdt_cat ! common epoch time
  DOUBLE PRECISION delta_sigma ! step in sigma 
  DOUBLE PRECISION v_inf(mulx)
  PUBLIC mulx, imip,imim,imi0,elm,unm,csinom,delnom,v_inf,tdt_cat,delta_sigma
  DOUBLE PRECISION moid_m(mulx), dnp_m(mulx), dnm_m(mulx) ! moid and nodal dist
  DOUBLE PRECISION sigq(mulx) ! sigma_Q,
  INTEGER, DIMENSION(mulx) :: imul ! obsolete copy of indexes, imul(j)=j or 0 
  PUBLIC moid_m, dnp_m ,dnm_m, sigq,imul
! common data, not public so far                            


CONTAINS                                                                     
! PUBLIC ROUTINES                                                       
!               f_multi
!               mult_input
!               fmuobs
!               fmupro
!               fmuplo
!               outmul                                                  
!                                                                       
! MODULE CONTAINS                                                       
!  ROUTINES                                                             
!                prop_sig                                                
!                int_step                                                  
! OUT OF MODULE                                              
!                graha_1   
!                                                                       
!                                                                       
!  HEADERS      
!        multiple_sol.o: \
!	../include/parobx.h90 \ only to dimension sorted data in difvin
!	../include/sysdep.h90 \
!	../suit/FUND_CONST.mod \
!	../suit/astrometric_observations.mod \
!	least_squares.o  \
!       pred_obs.o 
!
! =======================================                               
!  F_MULTI                                                               
! =======================================                               
! interface for multiple solutions               
! version with adaptive stepsize                                        
! ===============INTERFACE========================                      
SUBROUTINE f_multi(batch,obsc,inic,ok,covc,         &
     &     el0,uncert,csinor,delnor,                  &
     &     mc, obs, obsw,sigma,imult,sig1,sig2) 
! ================INPUT===========================                      
  LOGICAL batch ! batch/interactive
  INTEGER mc ! number of observations
! new data types
  TYPE(ast_obs),DIMENSION(mc) :: obs
  TYPE(ast_wbsr),DIMENSION(mc) :: obsw
! initial conditions: epoch, elements, abs. magnitude, gmag             
  TYPE(orbit_elem), INTENT(IN) :: el0 ! 
! normal and covariance matrices, norms of residuals, of last correction
  TYPE(orb_uncert), INTENT(IN) :: uncert ! 
  DOUBLE PRECISION, INTENT(IN) :: csinor,delnor 
  LOGICAL obsc,inic,covc ! logical flags  
! ========INPUT, but can also be assigned inside ===================
! number of alternate solutions
  INTEGER, INTENT(INOUT) ::  imult 
! max of sigma along the line                                           
  DOUBLE PRECISION, INTENT(INOUT) ::  sigma 
! ===============OUTPUT=========================                
! success flag                                                          
  LOGICAL, INTENT(OUT) ::  ok 
  DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: sig1,sig2 ! extreme values 
!               of sigma_Q at imim, imip
! ==============END INTERFACE===================== 
  DOUBLE PRECISION rescov  ! function to compute sigma_q
! temporary orbit elements
  TYPE(orbit_elem) eltmp
! weak direction, rms along it                                          
  DOUBLE PRECISION wdir(6),wdir0(6),sdir 
! no obs, used                                                          
  INTEGER nused 
! ======== differential correction flags and controls ======            
! scalar temporaries                                                    
  DOUBLE PRECISION dn,sigmam 
! loop indexes                                                          
  INTEGER i, j, imi,jj 
! success flag for diff. correction                                     
  LOGICAL succ 
! rms magnitudes                                             
  DOUBLE PRECISION rmsh  
  INTEGER iconv(mulx) 
! hack for catalog                                                      
  INTEGER iunctc,le 
  CHARACTER*20 astna0 
! control of integration method                                         
  INTEGER imint,k 
  DOUBLE PRECISION hh,ratio,ecc
  TYPE(orbit_elem) :: eqf, elc 
  LOGICAL fail 
! control of bizarre orbits                                             
  LOGICAl bizarre 
  INTEGER iun,iundump,iunint 
  DOUBLE PRECISION, DIMENSION(6) :: units ! for scaling
  iun=abs(iun_log) 
  IF(verb_mul.ge.30)THEN 
     iundump=abs(iun_log) 
  ELSE 
     iundump=-abs(iun_log) 
  ENDIF
  IF(verb_mul.ge.20)THEN 
     iunint=abs(iun_log) 
  ELSE 
     iunint=-abs(iun_log) 
  ENDIF
! ===================================================================== 
! check availability of observations and initial condition              
  IF(.not.obsc)THEN 
     WRITE(*,*)'f_multi: NO OBSERVATIONS' 
     ok=.false. 
     RETURN 
  ENDIF
  CALL chereq(2,inic,covc,el0%t,iunint,ok) 
  IF(.not.ok)THEN 
     WRITE(*,*)' f_multi: no initial data ',inic,covc,el0%t 
     RETURN 
  ENDIF
! ===================================================================== 
! check availability of JPL ephemerides and ET-UT table                 
  CALL chetim(obs(1)%time_tdt,obs(mc)%time_tdt,ok) 
  IF(.not.ok) THEN 
     WRITE(*,*)' f_multi: no JPL ephem/ET-UT table ',obs(1)%time_tdt,obs(mc)%time_tdt
     RETURN 
  ENDIF
! ===================================================================== 
! compute line of variations                                            
  CALL weak_dir(uncert%g,wdir,sdir,iunint,el0%coo,el0%coord,units)
  wdir0=wdir 
! input parameters of segment on the variations line                    
  IF(batch)THEN 
! sigma and imult have to be passed in the call                         
     IF(2*imult+1.gt.mulx)THEN 
        WRITE(*,*)' too many; max is ',mulx 
        STOP ' mmulti: too many mutsol' 
     ENDIF
  ELSE 
     WRITE(*,*)' how many sigma?' 
     READ(*,*) sigma 
259  WRITE(*,*)' how many steps on each side?' 
     READ(*,*) imult 
     WRITE(iun,*)' sigma=',sigma,'  in ',imult,' steps' 
     IF(2*imult+1.gt.mulx)THEN 
        WRITE(*,*)' too many; max is ',mulx 
        GOTO 259 
     ENDIF
  ENDIF
! ===================================================================== 
! nominal stepsize (in the sigma space)                                 
  dn=1.d0/float(imult) 
! store information on stepsize
  delta_sigma=dn*sigma
! ===================================================================== 
! nominal solution at the center of the list                            
  imi0=imult+1 
  elm(imi0)=el0 
  unm(imi0)=uncert
  csinom(imi0)=csinor 
  delnom(imi0)=delnor 
  sigq(imi0)=0.d0 
! orbital distance                                                      
  CALL nomoid(el0%t,elm(imi0),moid_m(imi0),dnp_m(imi0),dnm_m(imi0))
  iconv(imi0)=0
! ===================================================================== 
! main loop on number of steps (positive side)                          
! ===================================================================== 
  DO 5 i=1,imult 
     imi=imult+i+1 
     IF(.not.batch)WRITE(*,*)' alternate solution no. ',imi 
     IF(.not.batch)WRITE(iun,*)' alternate solution no. ',imi 
     CALL prop_sig(batch,elm(imi-1),elc,dn,sigma,mc,obs,obsw,wdir,sdir,units,fail)
! check for hyperbolic                                                  
     IF(fail)THEN 
        IF(.not.batch)WRITE(*,*)'step ',imi,' failed' 
!        IF(.not.batch)WRITE(*,*)elc 
        WRITE(iun,*)'fail in prop_sig, imi= ',imi,dn
        imi=imi-1
        GOTO 6 
     ELSEIF(bizarre(elc,ecc))THEN 
        IF(.not.batch)WRITE(*,*)'step ',imi,' bizarre', 'ecc=', ecc 
        WRITE(iun,*)'bizarre out of  prop_sig, imi= ',imi,ecc,dn
        imi=imi-1 
        GOTO 6 
     ELSE 
! constrained differential corrections:
         CALL constr_fit(mc,obs,obsw,elc,wdir,elm(imi),unm(imi),     &
     &       csinom(imi),delnom(imi),rmsh,nused,succ)
! exit if not convergent                                                
        IF(.not.succ) THEN 
           WRITE(iun,*)'fail in constr_fit, imi= ',imi,dn
           imi=imi-1 
           GOTO 6 
        ENDIF
! orbital distance                                                      
        CALL nomoid(elm(imi)%t,elm(imi),moid_m(imi),dnp_m(imi),dnm_m(imi))
        iconv(imi)=0
! check for sigQ.le.sigma                                               
        sigq(imi)=sqrt(abs(csinom(imi)**2-csinom(imi0)**2)*nused)
        IF(csinom(imi).lt.csinom(imi0))sigq(imi)=-sigq(imi)
        IF(batch.and.sigq(imi).gt.sigma)THEN 
           WRITE(iun,*)'too large sigma_q ',imi,sigq(imi)
           GOTO 6 
        ENDIF
     ENDIF
5 END DO
! ===================================================================== 
! line of variations, negative sigma
6 imip=imi 
  wdir=wdir0
  DO 7 i=1,imult 
! ===================================================================== 
! main loop on number of steps (negative side)                          
! ===================================================================== 
     imi=imult-i+1 
     IF(.not.batch)WRITE(*,*)' alternate solution no. ',imi 
     IF(.not.batch)WRITE(iun,*)' alternate solution no. ',imi 
     sigmam=-sigma 
     CALL prop_sig(batch,elm(imi+1),elc,dn,sigmam,mc,obs,obsw,wdir,sdir,units,fail) 
     IF(fail)THEN ! check for divergence
        IF(.not.batch)WRITE(*,*)'step ',imi,' failed' 
        WRITE(iun,*)'fail in prop_sig, imi= ',imi,dn
        imi=imi+1 
        GOTO 8 
     ELSEIF(bizarre(elc,ecc))THEN ! check for hyperbolic
        IF(.not.batch)WRITE(*,*)'step ',imi,' bizarre' 
        IF(.not.batch)WRITE(*,*)elc
        WRITE(iun,*)'bizarre out of prop_sig, imi= ',imi,ecc,dn
        imi=imi+1 
        GOTO 8 
     ELSE 
! differential corrections:                                             
! constrained differential corrections:
         CALL constr_fit(mc,obs,obsw,elc,wdir,elm(imi),unm(imi),     &
     &       csinom(imi),delnom(imi),rmsh,nused,succ)           
! exit if not convergent                                                
        IF(.not.succ)THEN 
           WRITE(iun,*)'fail in constr_fit, imi= ',imi,dn
           imi=imi+1 
           GOTO 8 
        ENDIF
! orbital distance                                                      
        CALL nomoid(elm(imi)%t,elm(imi),moid_m(imi),dnp_m(imi),dnm_m(imi))
        iconv(imi)=0
        sigq(imi)=sqrt(abs(csinom(imi)**2-csinom(imi0)**2)*nused)
        IF(csinom(imi).lt.csinom(imi0))sigq(imi)=-sigq(imi)
        IF(batch.and.sigq(imi).gt.sigma)THEN 
           WRITE(iun,*)'too large sigma_q ',imi,sigq(imi)
           GOTO 8 
        ENDIF
     ENDIF
7 ENDDO
! ===================================================================== 
! summary table                                                         
! ===================================================================== 
8 imim=imi 
! extreme values of sigma_Q
  IF(PRESENT(sig1))sig1=sigq(imim)
  IF(PRESENT(sig2))sig2=sigq(imip)
! catalog reference time
  tdt_cat=el0%t
  IF(batch)RETURN 
  CALL tee(iun,'SUMMARY OF MULTIPLE SOLUTIONS=')
  IF(el0%coo.eq.'EQU')THEN 
     CALL tee(iun,'no       a      h      k      p      q      lambda=')
  ELSEIF(el0%coo.eq.'CAR')THEN 
     CALL tee(iun,'no       x      y      z    xdot    ydot    zdot=')
  ELSEIF(el0%coo.eq.'KEP')THEN 
     CALL tee(iun,'no       a      e      I    Omeg    omeg    mean.an=')
  ELSEIF(el0%coo.eq.'COM')THEN
    CALL tee(iun,'no        q      e      I    Omeg    omeg    t.peri=')
  ELSEIF(el0%coo.eq.'ATT')THEN
    CALL tee(iun,'no      alpha  delta   adot  ddot     r      rdot=')
  ENDIF
  CALL filopn(iunctc,'mult.ctc','unknown') 
  CALL wromlh(iunctc,'ECLM','J2000') 
  DO i=imim,imip 
     WRITE(*,144)i,elm(i)%coord 
     WRITE(iun,144)i,elm(i)%coord
144  FORMAT(i5,5f12.8,f15.8) 
  ENDDO
  CALL tee(iun,'no  RMS ,lastcor,  magn,  MOID ,nod+,nod-, sigQ=') 
  DO i=imim,imip 
     WRITE(*,145)i,csinom(i),delnom(i),elm(i)%h_mag,moid_m(i),dnp_m(i),dnm_m(i),iconv(i),sigq(i)       
     WRITE(iun,145)i,csinom(i),delnom(i),elm(i)%h_mag,moid_m(i),dnp_m(i),dnm_m(i),iconv(i),sigq(i)
145  FORMAT(i5,1x,1p,e13.5,e11.3,2x,0p,f5.2,1x,f8.5,1x,f8.5,1x,f8.5,1x,i2,1x,f7.3)                 
     WRITE(astna0,111)i 
111  FORMAT('mult_',i4) 
     CALL rmsp(astna0,le) 
     CALL write_elems(elm(i),astna0(1:le),'ML',' ',iunctc,unm(i))
  ENDDO
  CALL filclo(iunctc,' ') 
END  SUBROUTINE f_multi
! =======================================                               
!  NMULTI interface to (distribution) f_multi                            
! =======================================                               
! multiple solutions for a differential correction problem              
! ===============INTERFACE========================                      
SUBROUTINE nmulti(nam0,elc,uncert,csinor,delnor,            &
     &     mc,obs,obsw,imult,sigma,moid_min) 
! ================INPUT=========================== 
  USE name_rules, ONLY: name_len
! number of observations
  INTEGER, INTENT(IN) ::  mc
! new data types
  TYPE(ast_obs),DIMENSION(mc), INTENT(IN) :: obs
  TYPE(ast_wbsr),DIMENSION(mc), INTENT(IN) :: obsw
! ====initial state at epoch=================                           
! asteroid name                                                         
  CHARACTER*(name_len), INTENT(IN) ::  nam0 
! initial conditions: epoch, elements, abs. magnitude, gmag             
  TYPE(orbit_elem), INTENT(IN) :: elc ! includes tc,eqc(6),hc,gmag 
! normal and covariance matrices
  TYPE(orb_uncert), INTENT(IN) :: uncert ! includes gc(6,6),cc(6,6)
! norms of residuals, of last correction
  DOUBLE PRECISION, INTENT(IN) ::  csinor,delnor 
! options for generation of multiple solutions                          
! number of multiple sol. is 2*imult+1                                  
  INTEGER, INTENT(IN) :: imult 
  INTEGER :: imult1 ! to avoid intent problems
! interval on LOV from -sigma to +sigma                                 
  DOUBLE PRECISION, INTENT(INOUT) ::  sigma 
! ==========OUTPUT==============                                        
  DOUBLE PRECISION, INTENT(OUT) ::  moid_min 
! first, last of the multiple solution indexes                          
!  INTEGER, INTENT(OUT) ::  imimout,imipout 
! renormalization factor                                                
  DOUBLE PRECISION rescov 
! reference solution                                                    
  INTEGER imi0 
! logical variables used for call to fmulti                             
  LOGICAL obsc,inic,covc,ok 
! ==============END INTERFACE=====================                      
! loop indexes                                                          
  INTEGER i 
  LOGICAL batch ! batch control
! catalog obj name (with underscore)                                    
  CHARACTER*20 astna0 
  INTEGER le
! output files                                                          
  CHARACTER*60 file 
  INTEGER iunctc,iunrep,iuncat 
  TYPE(orbit_elem) :: elk
  INTEGER fail_flag
! =============================================                         
! set flags to true to avoid checks                                     
  obsc=.true. 
  inic=.true. 
  covc=.true. 
! batch mode                                                            
  batch=.true. 
  imult1=imult
! ******************                                                    
! call to f_multi, distribution version, but with batch mode             
  CALL f_multi(batch,obsc,inic,ok,covc,elc,uncert,csinor,delnor,mc,obs,obsw,sigma,imult1)
! ================================================                      
! output                                                                
! multiple solutions catalog: multiline format                          
  CALL filnam('multsol',nam0,'ctc',file,le) 
  CALL filopn(iunctc,file(1:le),'unknown') 
  CALL wromlh(iunctc,'ECLM','J2000') 
! single line format KEP                                                    
  CALL filnam('multsol',nam0,'cat',file,le) 
  CALL filopn(iuncat,file(1:le),'unknown') 
  CALL wro1lh(iuncat,'ECLM','J2000','KEP') 
! report file (with MOID, RMS, etc.)                                    
  CALL filnam('multsol',nam0,'mrep',file,le) 
  CALL filopn(iunrep,file(1:le),'unknown') 
! information to be passed to resret2
  WRITE(iunrep,199)imult,sigma 
199 FORMAT('MULTIPLE SOLUTIONS SUMMARY, imult=',i4,' maxsigma=',f5.2) 
  WRITE(iunrep,299)scaling_lov,second_lov
299 FORMAT('Scaling LOV=',L1/'Second  LOV=',L1)
! header
  WRITE(iunrep,*)'no  RMS ,lastcor,  magn,  MOID ,nod+,nod-, sigQ' 
! loop on multiple solutions effectively computed                       
  moid_min=1.d3 
  DO i=imim,imip 
! orbital distance                                                      
     moid_min=min(moid_m(i),moid_min) 
     WRITE(iunrep,145)i,csinom(i),delnom(i),elm(i)%h_mag,moid_m(i),dnp_m(i),dnm_m(i),0,sigq(i)
145  FORMAT(i5,1x,1p,e13.5,e11.3,2x,0p,f5.2,1x,f8.5,1x,f10.5,1x,f10.5,1x,i2,1x,f7.3)
     astna0=' '
     WRITE(astna0,111)i 
111  FORMAT('mult_',i4) 
     CALL rmsp(astna0,le) 
! output multiline catalog                                              
     CALL write_elems(elm(i),astna0(1:le),'ML',' ',iunctc,unm(i))
! output one line catalog, if possible
     CALL coo_cha(elm(i),'KEP',elk,fail_flag)
     IF(fail_flag.ge.4)THEN
        WRITE(iun_log,*)' hyperbolic orbit in output ', fail_flag, i     
     ELSE
        CALL write_elems(elk,astna0(1:le),'1L',' ',iuncat)
     ENDIf
  ENDDO
  CALL filclo(iunctc,' ') 
  CALL filclo(iunrep,' ') 
  CALL filclo(iuncat,' ') 
!=============================================                          
END SUBROUTINE nmulti
!
! Copyright 2006, The Orbfit Consortium                                 
! ===================================================================   
! STEP_FIT2 differential corrector slippping along the LOV by steps  
! ===================================================================   
! version 3.3.2, 10 January 2006; with prop_sig for smooth slipping
!                                        
! Input: m observations number                                          
!        obs observations
!        obsw weights/bias
! Input/Output corrected:
!        el0 asteroid elements (replaced by corrected if successful)
!        unc0 includes covariance matrix and inverse
!        csino0 residuals norm
!        delno0 differential corrections norm
! Output 
!        nused  number of non-discarded obs
!        succ success flag; the solution is nominal, or a better one along LOV, 
!           if succ=.true.   
! ============= REMARK ===============================================  
! The weights are constant in this routine                              
! =============INTERFACE===== ========================================= 
SUBROUTINE step_fit2(m,obs,obsw,el0,unc0,csino0,delno0,nused,succ)
! ===================================================================== 
  USE pred_obs
  USE close_app, ONLY: kill_propag
! ================input data==========================
  INTEGER, INTENT(IN) ::  m ! number of observations
  TYPE(ast_obs), DIMENSION(m), INTENT(IN)  :: obs
  TYPE(ast_wbsr), DIMENSION(m), INTENT(INOUT) :: obsw 
! ================input/output ==========================
  TYPE(orbit_elem), INTENT(INOUT) :: el0  ! epoch time, initial elements       
                     ! assumed to be on the LOV, anyway check
  TYPE(orb_uncert), INTENT(INOUT)  :: unc0 ! normal and covar. matrix
  DOUBLE PRECISION, INTENT(INOUT) :: delno0,csino0 ! corr,res norm 
! ================output ==========================
  INTEGER, INTENT(OUT) ::  nused ! no obs. used 
  LOGICAL, INTENT(OUT) ::  succ ! success flag 
! =============END INTERFACE============================================
  DOUBLE PRECISION :: deq6_lov(6),deq6_const(6),ecc ! diff. corrections, eccentricity 
  DOUBLE PRECISION :: peq(6),deq6(6),deq(6),sdir, units(6) ! LOV vector, correction, step
  TYPE(orbit_elem) :: elc  ! corrected elements
  TYPE(orb_uncert) :: uncert ! normal and covar. matrix
  DOUBLE PRECISION :: delnor,csinor ! corr,res norm 
  TYPE(orbit_elem) elb ! temporary elements
! input data sorted 
  INCLUDE 'parobx.h90' 
  TYPE(ast_obs),DIMENSION(nobx) :: obs_s
  TYPE(ast_wbsr),DIMENSION(nobx) :: obsw_s        
  INTEGER iposs(nobx), nsolv
  INTEGER icor6(6) ! variables to be solved
  DOUBLE PRECISION csino1, delnor_lov, delnor_const, snorm, rescov, mu, dn
  LOGICAL lov_step, bizarre, fail
  INTEGER it,itst,itc,i,j, ir ! loop indexes
  INTEGER, PARAMETER :: itcmax=10 ! max no of constrained iterations
  INTEGER, PARAMETER :: itstmax=20 ! max no of iterations along LOV
! ==================================================================== 
! number of solve-for variables                                         
  icor6=1 ! for final matrices
! definition of blocks if necessary
  CALL blockset(m,obs,obsw)
! sort of times and reordering of obs. and weights                      
  CALL sort_obs(el0%t,obs,obsw,m,iposs,obs_s,obsw_s)                         
! Initialisation with starting value for elements                       
  elc=el0;  csino1=csino0; succ=.false.
  if(verb_dif.gt.9)write(iun_log,*)'starting values ', el0 
  if(verb_dif.gt.19) write(*,*)'starting values ',  el0 
! ================== main loop on step iterations=========================
  it=0 ! counter of all iterations
  itst=0 ! counter of LOV steps only (limited by itstmax)
  itc=0  ! constrained correction iterations (limited by itcmax)
60 CONTINUE 
  CALL sin_cor(m,obs_s,obsw_s,elc,icor6,iun_log,.true.,delnor,csinor,uncert,deq6)
  IF(.not.uncert%succ) GOTO 9
  it=it+1 
  CALL weak_dir(uncert%g,peq,sdir,-1,elc%coo,elc%coord,units)
  deq6=deq6/units ! diff. corr. in scaled coords
  IF(DOT_PRODUCT(peq,deq6).lt.0.d0) peq=-peq ! sign rule to decrease Q 
  deq6_lov=DOT_PRODUCT(peq,deq6)*peq ! component along the LOV
  deq6_const=deq6-deq6_lov ! component orthogonal to the LOV
  deq6_lov=deq6_lov*units ! in unscaled coordinates
  deq6_const=deq6_const*units ! in unscaled coordinates
  delnor_lov=snorm(deq6_lov,uncert%c,6,6) ! how big the diff.cor. along the LOV
  delnor_const=snorm(deq6_const,uncert%c,6,6) ! how big the diff.cor. normal to the LOV
! Check if we need another constrained iteration
  IF(delnor_const.ge.del_constr)THEN
! yes, need to get back to the LOV
     itc=itc+1
     IF(verb_dif.ge.19)write(*,*)' constr.corr., delnor_const=',delnor_const
     IF(verb_dif.ge.9)write(iun_log,*)' constr.corr., delnor_const=',delnor_const
     lov_step=.false.
! use constrained corrections  
     deq=deq6_const
  ELSEIF(delnor_lov.lt.delcr)THEN
! no, and also no need to move along the LOV
     succ=.true.
     IF(verb_dif.ge.19)write(*,*)' success, nominal solution, RMS, norms ',delnor_const,delnor_lov
     write(iun_log,*)' success, nominal solution, RMS, norms ',delnor_const,delnor_lov
     goto 70
  ELSE
! no, is the current result better than the stored one?
     IF(csinor.lt.csino1)THEN ! potential LOV output
        el0=elc; unc0=uncert; csino0=csinor; delno0=delnor; csino1=csinor
        succ=.true.
        IF(verb_dif.ge.19)write(*,*)' success, LOV solution, RMS, norms ',csinor,delnor_const,delnor_lov
        write(iun_log,*)' success, LOV solution, RMS, norms ',csinor,delnor_const,delnor_lov
     ENDIF
! but it is useful to try and move along the LOV
     itc=0 !restart count of constr. steps for next time
     itst=itst+1 ! LOV step
! use unconstrained corrections, of limited length in the weak direction
     IF(verb_dif.ge.19)write(*,*)' full corr., delnor_lov=',delnor_lov
     IF(verb_dif.ge.9)write(iun_log,*)' full corr., delnor_lov=',delnor_lov
     lov_step=.true.
     IF(delnor_lov.gt.step_sig)THEN
!        dn=step_sig
        deq=deq6_const+deq6_lov*(step_sig/delnor_lov)
     ELSE
!        dn=delnor_lov
        deq=deq6_const+deq6_lov
     ENDIF
! ALTERNATIVE METHOD
!     elb=elc
!     CALL prop_sig(.true.,elb,elc,dn,sdir,m,obs,obsw,peq,sdir,units,fail)
! END ALTERNATIVE METHOD
  ENDIF
! relaxation loop
  elb=elc
  DO ir=1,4
! Update solution 
     elb%coord=elc%coord+deq 
     IF(bizarre(elb,ecc))THEN
        IF(verb_dif.ge.19)WRITE(*,*)' constr_fit: short step to avoid ecc=',ecc
        IF(verb_dif.ge.9)WRITE(iun_log,*)' constr_fit: short step to avoid ecc=',ecc
        deq=deq/2.d0
        elb%coord=elc%coord+deq 
     ELSE
        EXIT
     ENDIF
  ENDDO
  elc=elb
  delnor=snorm(deq,uncert%c,6,6) ! full norm of the actual correction used
  IF(verb_dif.ge.9)write(iun_log,200)it,itc,csinor,delnor,elc%coord
  IF(verb_dif.ge.19)write(*,200)it,itc,csinor,delnor,elc%coord 
200 FORMAT(' *** iteration ',i3,2x,I3,' RMS residuals =',1p,d12.4,    &
   &        '   norm corr =',d12.4,'  new elem values:'/0p,6f13.7/)
! control against hyperbolic and bizarre orbits                         
  IF(bizarre(elc,ecc))THEN 
     IF(verb_dif.ge.19)write(*,*)' step_fit: iter. ',it,' bizarre; e=',ecc
     IF(verb_dif.ge.9)write(iun_log,*)' step_fit: iter. ',it,' bizarre; e=',ecc
     RETURN
  ENDIF
! check if iterations can go on
  IF(itc.gt.itcmax.or.itst.gt.itstmax)THEN
     GOTO 70 !interrupt attempt
  ELSE
     GOTO 60 ! keep trying
  ENDIF
70 CONTINUE
  IF(.not.succ)THEN
     IF(verb_dif.ge.19)THEN
        WRITE(*,*)' non convergent after ',it,itst,itc,' iterations ' 
     ELSEIF(verb_dif.gt.2)THEN
        WRITE(iun_log,*)' non convergent after ',it,itst,itc,' iterations '
     ENDIF
  ENDIF
! =========covariance (and normal matrix) rescaling==============
  nsolv=6 
  nused=0 
  DO i=1,m 
     IF(obsw_s(i)%sel_coord.gt.0)THEN
        IF(obs_s(i)%type.eq.'O'.or.obs_s(i)%type.eq.'S')THEN
           nused=nused+2
        ELSEIF(obs_s(i)%type.eq.'R'.or.obs_s(i)%type.eq.'V')THEN
           nused=nused+1
        ENDIF
     ENDIF
  ENDDO
  mu=rescov(nsolv,nused,csinor) 
! apply rescaling to both covariance and normal matrix 
  uncert%g=uncert%g*mu**2 
  uncert%c=uncert%c/mu**2 
! Output final result, with norm of the residuals                            
  IF(verb_dif.ge.9.and.verb_dif.lt.19)write(*,201)it,csinor,delnor 
201  format(' done constr. iter. ',i3,'  RMS=',1p,d12.4,' last corr. =',d12.4)
!  IF(verb_dif.lt.9.and.verb_dif.gt.2)write(iun,201)it,csinor,delnor 
! reordering the residuals for output                                   
  CALL unsort_obs(iposs,m,obsw_s,obsw)
  RETURN ! normal ending
9 CONTINUE ! impossible normal matrix inversion: give up
  IF(verb_mul.ge.9)THEN
     WRITE(*,*)' inversion failed '
     WRITE(iun_log,*)' inversion failed '
  ENDIF
  succ=.false.
  csinor=csino0
END SUBROUTINE step_fit2


! =======================================                               
! mult_input initializes storage of multiple solution data              
! =======================================                               
SUBROUTINE mult_input(catname,ok)  
  USE name_rules
! ------------INPUT------------------                                   
! file with catalog of multiple solutions                               
  CHARACTER*160, INTENT(IN) :: catname 
! ------------OUTPUT-----------------                                   

! succesful input                                                       
  LOGICAL, INTENT(OUT) :: ok 
! ------------END INTERFACE--------- 
! index range of multiple solution
  INTEGER :: m1,m2,m0 
  DOUBLE PRECISION v_infty 
! -----for call to rdorb----                                            
! names                                                                 
  CHARACTER*(idname_len) name0 
  CHARACTER*(name_len) name1 
  INTEGER le
 CHARACTER*160 catname1
! magnitude slope parameter, mass, epoch                                
  double precision sl,mass,t 
! type of orbital elements (KEP/EQU/CAR),                               
  character*3 eltype 
! reference system type (EQUM/EQUT/ECLM),                               
  character*4 rsys 
! epoch specification (J2000/OFDATE)                                    
  character*6 epoch 
! record number                                                         
  integer no 
! avalaibility of covariance and normal matrices                        
  logical defcov,defnor 
! end of file                                                           
  LOGICAL eof 
! indexes of multiple solutions                                         
  INTEGER imul(mulx),norb,j 
! velocity w.r.to Earth for each orbit                                  
  DOUBLE PRECISION v_max,v_min 
! temporary input arrays                                                
  DOUBLE PRECISION eq(6) 
  DOUBLE PRECISION g(6,6),c(6,6),h 
! loop indexes                                                          
  INTEGER i 
! system dependencies                                                   
  INCLUDE 'sysdep.h90' 
! -------------------------------------------------                     
  v_min=100.d0 
  v_max=0.d0 
! opening and reading multiple solution catalog                         
  catname1=catname
  CALL rmsp(catname1,le) 
  INQUIRE(file=catname1(1:le),exist=ok) 
  IF(.not.ok)THEN 
     WRITE(*,*)' file ',catname1(1:le),' not found' 
     RETURN 
  ENDIF
  CALL oporbf(catname1(1:le),0) 
  DO i=1,mulx 
     CALL rdorb(name0,eq,eltype,t,g,defcov,c,defnor,h,sl,mass,rsys,epoch,no,eof)
     IF(eof)THEN 
        norb=i-1 
        GOTO 2 
     ENDIF
! control on time                                                       
     IF(i.eq.1)THEN 
        tdt_cat=t 
     ELSE 
        IF(t.ne.tdt_cat)THEN 
           WRITE(*,*)'mult_input: time discrepancy from tcat=',tdt_cat 
           WRITE(*,*)'mult_input: at record ',i,' t=',t 
           ok=.false. 
           RETURN 
        ENDIF
     ENDIF
! availability of covariance                                            
     IF (.NOT.(defnor.AND.defcov)) THEN 
        WRITE(*,*)'mult_input',  name0, ' matrices not avalaible' 
        ok=.false. 
        RETURN 
     ENDIF
! handling of name; but note that the multiple solutions are assumed to 
! and with consecutive indexes!!! that is imul(i)=imul(1)+i-1           
     CALL splinam(name0,name1,imul(i)) 
     IF(imul(i).ne.imul(1)+i-1)THEN 
        WRITE(*,*)'mult_input: indexes not in order: ',             &
     &           imul(i),' at record ',i,' should be ',imul(1)+i-1      
!            RETURN                                                     
     ENDIF
     j=imul(i) 
! copy in output arrays
     elm(j)=undefined_orbit_elem
     IF(rhs.EQ.2)THEN
        elm(j)%center=3
     END IF
     elm(j)%coord= eq
     elm(j)%t=t
     elm(j)%coo=eltype 
     IF(h.gt.0.d0)elm%mag_set=.true.
     elm(j)%h_mag=h 
     unm(j)%c=c
     unm(j)%g=g
! v_infinity computation                                                
     v_inf(i)=v_infty(elm(j)) 
     v_min=min(v_min,v_inf(i)) 
     v_max=max(v_max,v_inf(i)) 
  ENDDO
! increase mulx                                                        
  WRITE(*,*)'mult_input: increase mulx, file is longer than',mulx 
2 CONTINUE 
  norb=i-1 
  ok=norb.gt.0 
  imim=imul(1) 
  imip=imul(norb) 
  WRITE(*,*)' mult_input: input of ',norb,' multiple solutions' 
  WRITE(*,*)' indexes between ', imim, ' and ', imip 
  WRITE(*,*)' max of v_infty ',v_max,' min ',v_min 
3 WRITE(*,*)' which one is the nominal?' 
  READ(*,*)imi0 
  IF(imi0.lt.imim.or.imi0.gt.imip)THEN 
     WRITE(*,*)' OUT OF RANGE ',imim,imip 
     GOTO 3 
  ELSE 
     WRITE(*,*)' nominal is ',imi0 
  ENDIF
END SUBROUTINE mult_input
!
! =======================================                               
!  FMUOBS                                                               
! =======================================                               
! output of multiple observations                                       
! ===============INTERFACE========================                      
SUBROUTINE fmuobs(type,ids,t1,tut1,sigma,         &
     &     aobs,dobs,iff,titnam,filnam,iun20) 
  USE pred_obs   
! sigma value                   
  DOUBLE PRECISION sigma 
! actual observations                                                   
  DOUBLE PRECISION aobs,dobs 
! strings with asteroid names, output unit                              
  CHARACTER*80 titnam 
  CHARACTER*60 filnam 
  INTEGER, INTENT(IN) ::  iun20 
! station code, flag for use of act.obs., observation type              
  INTEGER, INTENT(IN) ::  ids,iff
  CHARACTER*(1) type
! target time (TDT, UTC)                                                
  DOUBLE PRECISION, INTENT(IN) :: t1,tut1 
! =================END INTERFACE====================                    
  DOUBLE PRECISION alm(mulx),dem(mulx)
!       ,adotm(mulx),ddotm(mulx),disv(mulx) 
  INTEGER ng,i,npop, icov, inl
  DOUBLE PRECISION eqm1(6,mulx),amagn(mulx),eq1(6),alpha,delta 
!      DOUBLE PRECISION amagn1,alpha1,delta1 ! for test comparison
  DOUBLE PRECISION :: adot,ddot ! proper motion
  DOUBLE PRECISION :: pha,dis,dsun,elo,gallat ! phase, distance to Earth, distance to Sun
  DOUBLE PRECISION :: gamad(2,2),axes(2,2),sig(2) ! covariance on sky plane
! ====================================================                  
  type='O' ! and does not handle radar observations
  inl=1 !
! first compute nominal prediction                                      
  CALL  predic_obs(elm(imi0),ids,t1,type,       &
     &        alpha,delta,amagn(imi0),inl,        &
     &        GAMAD=gamad,SIG=sig,AXES=axes,    &
     &        ADOT0=adot,DDOT0=ddot,DIS0=dis,ELO0=elo)
  WRITE(*,*)' nominal solution ' 
  icov=1
  CALL outobc(iun20,type,ids,tut1,alpha,delta,amagn(imi0),adot,ddot,&
     &     elo,dis,icov,gamad,sig,axes)                                 
  alm(imi0)=0.d0 
  dem(imi0)=0.d0 
! initialize revolution counter                                         
  ng=0 
! loop on existing multiple solutions: first forward, then backward     
  DO  i=imi0+1,imip 
     CALL  predic_obs(elm(i),ids,t1,type,alm(i),dem(i),amagn(i),inl,   &
     &        ADOT0=adot,DDOT0=ddot,DIS0=dis,ELO0=elo)
     IF(alm(i).gt.pig)alm(i)=alm(i)-dpig 
     WRITE(*,*)' alternate obs.no. ',i
     icov=1 
     CALL outobc(iun20,type,ids,tut1,alm(i),dem(i),amagn(i),        &
     &     adot,ddot,elo,dis,icov,gamad,sig,axes)                       
     alm(i)=alm(i)-alpha 
     dem(i)=dem(i)-delta 
! update revolution counter                                             
     CALL angupd(alm(i),alm(i-1),ng) 
  ENDDO
  DO  i=imi0-1,imim,-1 
     CALL  predic_obs(elm(i),ids,t1,type,alm(i),dem(i),amagn(i),inl,  &
     &        ADOT0=adot,DDOT0=ddot,DIS0=dis,ELO0=elo)
     WRITE(*,*)' alternate obs.no. ',i 
     icov=1
     CALL outobc(iun20,type,ids,tut1,alm(i),dem(i),amagn(i),        &
     &     adot,ddot,elo,dis,icov,gamad,sig,axes)                       
     alm(i)=alm(i)-alpha 
     dem(i)=dem(i)-delta 
! keep count of lost revolutions in alpha                               
     CALL angupd(alm(i),alm(i+1),ng) 
  END DO
! ===============================================                       
! output multiple prediction of observations                            
  CALL outmul(titnam,filnam,tut1,sigma,alpha,delta,                 &
     &        alm,dem,amagn,imim,imip,imi0,iff,aobs,dobs,type)            
END SUBROUTINE fmuobs 

!                                                                       
! ====================================================                  
! FMUPRO multiple state propagation for FITOBS                          
! ====================================================                  
SUBROUTINE fmupro(iun20,tr) 
  USE propag_state
! =================INPUT========================================= 
  INTEGER, INTENT(IN) :: iun20! output units
  DOUBLE PRECISION, INTENT(IN) ::tr ! target time
! ================END INTERFACE==========================               
! ======== output moid =====================                            
  INTEGER iconv(mulx)  ! obsolete, for compatibility with old nomoid
  INTEGER j,i ! loop indexes        
! ===================================================================== 
! main loop                                                             
! propagation to time tr                                                
  DO j=imim,imip 
     WRITE(*,*)' orbit ',j 
     CALL pro_ele(elm(j),tr,elm(j),unm(j),unm(j)) 
! orbital distance                                                      
     CALL nomoid(tr,elm(j),moid_m(j),dnp_m(j),dnm_m(j))
     iconv(j)=0
  ENDDO
! ===================================================================== 
! summary table                                                         
! ===================================================================== 
  CALL tee(iun20,'SUMMARY OF MULTIPLE SOLUTIONS=') 
  WRITE(iun20,223) tr
  WRITE(*,223) tr 
223 FORMAT(' elements at time ',f8.1,' (MJD):') 
  IF(elm(1)%coo.eq.'EQU')THEN 
     CALL tee(iun20,'no       a      h      k      p      q      lambda=')
  ELSEIF(elm(1)%coo.eq.'CAR')THEN 
     CALL tee(iun20,'no       x      y      z    xdot    ydot    zdot=')
  ELSEIF(elm(1)%coo.eq.'KEP')THEN 
     CALL tee(iun20,'no       a      e      I    Omeg    omeg    mean.an=')
  ELSEIF(elm(1)%coo.eq.'COM')THEN
    CALL tee(iun20,'no        q      e      I    Omeg    omeg    t.peri=')
  ELSEIF(elm(1)%coo.eq.'ATT')THEN
    CALL tee(iun20,'no      alpha  delta   adot  ddot     r      rdot=')
  ENDIF
  DO i=imim,imip 
     WRITE(*,144)i,elm(i)%coord 
144  FORMAT(i5,6f12.8) 
     WRITE(iun20,144)i,elm(i)%coord
  ENDDO

  CALL tee(iun20,'no.,  magn,  MOID ,  nod+  ,  nod-=') 
  DO i=imim,imip 
     WRITE(*,145)i,(i),moid_m(i),dnp_m(i),dnm_m(i),iconv(i) 
     WRITE(iun20,145)i,elm(i)%h_mag,moid_m(i),dnp_m(i),dnm_m(i),iconv(i) 
145  FORMAT(i4,2x,f5.2,1x,f8.5,1x,f8.5,1x,f8.5,1x,i2) 
  ENDDO
END SUBROUTINE fmupro
!                                                                       
! =======================================                               
!  FMUPLO                                                               
! =======================================                               
! graphic output of multiple orbits                                     
! ===============INTERFACE========================                      
SUBROUTINE fmuplo(titnam,sigma) 
  USE util_suit
  CHARACTER*80, INTENT(IN) :: titnam 
  DOUBLE PRECISION, INTENT(IN) :: sigma
! ============END INTERFACE==================                           
  DOUBLE PRECISION a(mulx),e(mulx),aa,ee ,q,qg,enne
  INTEGER i,j,k, fail_flag
  TYPE(orbit_elem) :: elcom
  CHARACTER*3 coo
  INTEGER icoo,numb
  CHARACTER*60 xlab,ylab
  CHARACTER*20 menunam ! characters for menu
! ===========================================
  menunam='plotcoo'  
  CALL menu(icoo,menunam,3,' which coordinates to plot?=',            &
     &      'keplerian a,e=',                                          &
     &      'cometary q,e=',                                           &
     &      'attributable r, rdot=')   
  IF(icoo.eq.0)RETURN
  numb=imip-imim+1
  coo=elm(1)%coo
  IF(icoo.eq.1)THEN      
     xlab='Semimajor axis (AU)'
     ylab='Eccentricity'             
! a-e plot of multiple solutions
     IF(coo.eq.'EQU')THEN
        DO i=1,numb 
           j=i+imim-1
           a(i)=elm(j)%coord(1) 
           e(i)=sqrt(elm(j)%coord(2)**2+elm(j)%coord(3)**2) 
        ENDDO
        aa=elm(imi0)%coord(1)
        ee=sqrt(elm(imi0)%coord(2)**2+elm(imi0)%coord(3)**2) 
     ELSEIF(coo.eq.'KEP')THEN
        DO i=1,numb 
           j=i+imim-1
           a(i)=elm(j)%coord(1) 
           e(i)=elm(j)%coord(2)
        ENDDO
        aa=elm(imi0)%coord(1)
        ee=elm(imi0)%coord(2)
     ELSE
        DO i=1,numb 
           j=i+imim-1
           CALL coo_cha(elm(j),'COM',elcom,fail_flag)
           IF(fail_flag.le.10.and.fail_flag.ne.4)THEN
              e(i)=elcom%coord(2)
              a(i)=abs(elcom%coord(1)/(1.d0-e(i)))
           ELSE
              WRITE(*,*)' fmuplo: conversion failure ', elcom%coo 
              RETURN
           ENDIF
        ENDDO
        CALL coo_cha(elm(imi0),'COM',elcom,fail_flag)
        IF(fail_flag.le.10.and.fail_flag.ne.4)THEN
           ee=elcom%coord(2)
           aa=abs(elcom%coord(1)/(1.d0-ee))
        ELSE
           WRITE(*,*)' fmuplo: conversion failure ', elcom%coo 
           RETURN
        ENDIF
     ENDIF
  ELSEIF(icoo.eq.2)THEN
     xlab='Perihelion distance (AU)'
     ylab='Eccentricity'     
     DO i=1,numb 
        j=i+imim-1
        CALL coo_cha(elm(j),'COM',elcom,fail_flag)
        IF(fail_flag.le.10.and.fail_flag.ne.4)THEN
           e(i)=elcom%coord(2)
           a(i)=elcom%coord(1)
        ELSE
           WRITE(*,*)' fmuplo: conversion failure ', fail_flag
           RETURN
        ENDIF
     ENDDO
     CALL coo_cha(elm(imi0),'COM',elcom,fail_flag)
     aa=elcom%coord(1)
     ee=elcom%coord(2)
  ELSEIF(icoo.eq.3)THEN
     xlab='Geocentric distance (AU)'
     ylab='Geocentric range rate (AU/d)'     
     DO k=1,numb 
        j=k+imim-1
        CALL coo_cha(elm(j),'ATT',elcom,fail_flag)
        IF(fail_flag.le.10.and.fail_flag.ne.4)THEN
           e(k)=elcom%coord(6)
           a(k)=elcom%coord(5)
        ELSE
           WRITE(*,*)' fmuplo: conversion failure ', fail_flag
           RETURN
        ENDIF
     ENDDO
     CALL coo_cha(elm(imi0),'ATT',elcom,fail_flag)
     ee=elcom%coord(6)
     aa=elcom%coord(5)
  ENDIF
  CALL ploae(elm(imi0)%t,a,e,aa,ee,sigma,numb,titnam,xlab,ylab)
   
END SUBROUTINE fmuplo 
! ===================================================================== 
!  OUTMUL                                                               
! ===================================================================== 
! output multiple observations                                          
! =============INTERFACE=============================================== 
SUBROUTINE outmul(titnam,filnam,t1,sigma,alpha,delta,             &
     &     alm,dem,hmagn,imloc,iploc,i0loc,iff,aobs,dobs,type)
! =============INPUT=================================================== 
! file name                                                             
  CHARACTER*80,INTENT(IN) :: titnam 
  CHARACTER*60,INTENT(IN) ::  filnam 
! first, last and central index control for closed curve,obs type             
  INTEGER,INTENT(IN) ::  imloc,iploc,i0loc, iff
  CHARACTER*(1), INTENT(OUT) :: type
! observation time MJD, sigma value, nominal prediction                 
  DOUBLE PRECISION,INTENT(IN) :: t1,sigma,alpha,delta 
! observation: predicted  value alpha, delta, magnitude, actual         
  DOUBLE PRECISION,INTENT(IN) ::  alm(iploc),dem(iploc),hmagn(iploc),aobs,dobs 
! =============END INTERFACE=========================================== 
! conversion to sessagesimal                                            
  DOUBLE PRECISION seca,secd 
  INTEGER inta,mina,intd,mind 
  CHARACTER*1 siga,sigd 
! conversion of time                                                    
  INTEGER iy,imo,iday 
  DOUBLE PRECISION hour 
! scalar temporaries, differences                                       
  DOUBLE PRECISION dee,daa,ado,ddo 
! file name                                                             
  INTEGER le 
  CHARACTER*90 file, filnam1 
! loop indexes, units                                                   
  INTEGER n,iun7 
! ======================================================================
! open output file
  filnam1=filnam                                                      
  CALL rmsp(filnam1,le) 
  file=filnam1(1:le)//'.cbd' 
  CALL filopn(iun7,file,'unknown') 
! date and sigma value                                                  
  CALL mjddat(t1,iday,imo,iy,hour) 
  WRITE(iun7,297)iday,imo,iy,hour,sigma 
297 FORMAT(i3,i3,i5,f8.4,f5.2) 
! line of variations                                                    
  DO n=imloc,iploc 
     daa=alpha+alm(n) 
     daa=mod(daa,dpig) 
     IF(daa.lt.0.d0)daa=daa+dpig 
     IF(daa.gt.dpig)daa=daa-dpig 
     daa=daa*degrad/15 
     IF(daa.lt.0.d0.or.daa.gt.24.d0)THEN 
        WRITE(*,*)' outmul: daa out of range ', daa 
     ENDIF
     CALL sessag(daa,siga,inta,mina,seca) 
     dee=(delta+dem(n))*degrad 
     CALL sessag(dee,sigd,intd,mind,secd) 
! proper motion in arcsec/hour                                          
!        ado=adotm(n)*secrad/24.d0 
!        ddo=ddotm(n)*secrad/24.d0 
! output                                                                
     IF(siga.eq.'+')siga=' ' 
     WRITE(iun7,396)n,siga,inta,mina,seca,sigd,intd,mind,secd,hmagn(n)
!     &       disv(n),ado,ddo,hmagn(n)                                   
396  FORMAT(i3,1x,a1,i2,1x,i2,1x,f4.1,2x,a1,i2,1x,i2,1x,f4.1,1x,f5.2)
!     &       1x,f8.5,1x,f8.2,1x,f8.2,                           
  ENDDO
  CALL filclo(iun7,' ') 
! graphics output                                                       
  IF(iff.eq.1)THEN 
     CALL plocbd(titnam,alpha,delta,sigma,t1,                       &
          &         alm(imloc),dem(imloc),iploc-imloc+1,type)
  ELSEIF(iff.eq.2)THEN 
     CALL ploobs(titnam,alpha,delta,sigma,t1,                       &
     &         alm(imloc),dem(imloc),iploc-imloc+1, aobs,dobs)
  ENDIF
END SUBROUTINE outmul
! =========================================                             
! PROP_SIG                                                               
! propagator in sigma space                                             
SUBROUTINE prop_sig(batch,el1,el2,dn,sigma,mc,obs,obsw,wdir,sdir,units,fail)
! ==============input================= 
  LOGICAL batch ! batch control
! current elements and epoch; stepsize factor, target sigma 
  TYPE(orbit_elem), INTENT(IN) :: el1            
  DOUBLE PRECISION, INTENT(IN) :: dn,sigma 
! ======observations====                                                
  INTEGER, INTENT(IN) :: mc ! number of obs
! new data types
  TYPE(ast_obs),DIMENSION(mc), INTENT(IN) :: obs
  TYPE(ast_wbsr),DIMENSION(mc), INTENT(IN) :: obsw
! weak direction, rms along it                                          
  DOUBLE PRECISION, DIMENSION(6), INTENT(INOUT) :: wdir, units
  DOUBLE PRECISION, INTENT(INOUT) :: sdir
! ==============output================   
  TYPE(orbit_elem), INTENT(OUT) :: el2 ! elements
  LOGICAl, INTENT(OUT) :: fail ! failure flag
!  ===========end interface============                                  
  DOUBLE PRECISION hh,ch,h,ridmax 
  INTEGER iun,imint 
! intermediate point: state, covariance, normal matr., RMS, norm corr   
  TYPE(orbit_elem) :: elt,eltp
  TYPE(orb_uncert) :: uncert 
  DOUBLE PRECISION :: csinor,delnor, wdir0(6) 
! obs,residuals control, flag (for difcor)         
  INTEGER inew, icor(6),ncor,itmaxold,j 
! success control is dummy                                              
  LOGICAL succ 
! use weak direction to find initial conditions                         
  wdir0=wdir
  hh=dn*sigma 
  IF(.not.batch.and.verb_mul.ge.9)THEN 
     iun=abs(iun_log) 
  ELSE 
     iun=-abs(iun_log) 
  ENDIF
! 1=Euler 2= RK2                                                        
  imint=2 
! iteration on stepsize reductions                                      
!     ridmax=30.d0                                                      
  ridmax=130.d0 
  h=hh 
  ch=0.d0 
  elt=el1 
! try to do the step at once...                                         
1 CONTINUE 
!  write(*,*)'prop_sig 0: h,a, sdir',h,elt%coord(1),sdir 
!  WRITE(*,*)'mmulti: ch,h, eqt ', ch,h,(eqt(j),j=1,3) 
!  WRITE(*,*)' prop-sig: about to enter int_step',tc, obs(mc)%time_tdt
  CALL int_step(elt,el2,h,imint,mc,obs,obsw,iun,wdir,sdir,units,fail) 
!  write(*,*)'prop_sig 1: a, sdir',elt%coord(1),sdir 
! if failed, half step                                                  
  IF(fail)THEN 
     IF(abs(h).ge.abs(hh)/ridmax)THEN 
        h=h/2.d0 
        IF(verb_mul.ge.20.and.iun.gt.0)THEN 
           WRITE(iun,*)'prop_sig: halving', h,hh 
        ENDIF
        GOTO 1 
     ELSE 
        IF(verb_mul.ge.5.and.iun.gt.0)THEN 
           WRITE(iun,*)'prop_sig: too many halvings ',h,hh,ridmax 
        ENDIF
        WRITE(*,*)'prop_sig: too many halvings ',h,hh,ridmax
        RETURN 
     ENDIF
  ELSE 
     ch=ch+h 
  ENDIF
  IF(abs(ch).lt.abs(hh))THEN 
     elt=el2 
! compute line of variations                                            
! compute vectorfield at intermediate point                             
     itmaxold=itmax 
     itmax=0 
     CALL whicor(0,icor,ncor,inew) 
! compute covariance and residual norm     
!  WRITE(*,*)' prop-sig: about to enter diff-cor',tc, obs(mc)%time_tdt 
     CALL diff_cor(mc,obs,obsw,elt,icor,iun,eltp,uncert,csinor,delnor,succ)
     itmax=itmaxold 
! compute weak direction and length                                     
     CALL weak_dir(uncert%g,wdir,sdir,iun,eltp%coo,eltp%coord,units) 
!     write(*,*)'prop_sig 2: a, sdir',elt%coord(1),sdir 
     IF(DOT_PRODUCT(wdir,wdir0).lt.0.d0)wdir=-wdir
     GOTO 1 
  ENDIF
  IF(iun.gt.0.and.verb_mul.gt.20)WRITE(iun,*)' prop_sig: stepsize ',h 
END SUBROUTINE prop_sig
!=================================================                      
! INT_STEP integration step for propagation along LOV                    
SUBROUTINE int_step(el0,el1,hh,imint,mc,obs,obsw,iun,wdir,sdir,units,fail)
! stepsize, integration method                                          
  DOUBLE PRECISION, INTENT(IN) ::  hh 
  INTEGER, INTENT(IN) :: imint 
! weak direction (in input computed at eq0)             
  DOUBLE PRECISION, INTENT(INOUT) :: wdir(6),units(6),sdir
! ======observations====                                                
! number of obs.
  INTEGER, INTENT(IN) ::  mc
 ! new data types
  TYPE(ast_obs),DIMENSION(mc), INTENT(IN) :: obs
  TYPE(ast_wbsr),DIMENSION(mc), INTENT(IN) :: obsw
! initial conditions: epoch, elements, new elements after step 
  TYPE(orbit_elem), INTENT(IN) :: el0
  TYPE(orbit_elem), INTENT(OUT) :: el1
! output unit                                                           
  INTEGER iun 
! failure flag                                                          
  LOGICAl fail
! =========end interface================                                
! intermediate point: state, covariance, normal matr., RMS, norm corr   
  TYPE(orbit_elem) ::  el12,el12p
  TYPE(orb_uncert) :: un12
  DOUBLE PRECISION :: csino12,delno12 
! weak direction (in input computed at eq0), length of axis             
  DOUBLE PRECISION wdir12(6),sdir12,wdirst(6),sdirst 
! obs,residuals control, flag (for difcor)         
  INTEGER inew, icor(6),ncor,itmaxold 
! success control is dummy,bizarre is impor                             
  LOGICAL succ,bizarre 
! loop indexes                                                          
  INTEGER j 
! controls of fixed point iteration                                     
  INTEGER itx,it 
  DOUBLE PRECISION eps,cosa,ang,dsi,e1,e12,direc,q,qg,enne 
! ================================================                      
  IF(imint.eq.1)THEN 
! Euler method
     el1=el0                                                          
     el1%coord=el0%coord+units*wdir*sdir*hh 
     fail=.false. 
  ELSEIF(imint.eq.2)THEN 
! Runge-Kutta-Gauss of order 2:                                         
     itx=5 
     eps=1.d-4 
! first store vectorfield at el0                            
     wdirst=wdir
     sdirst=sdir 
! compute intermediate point (first guess)                              
     el12=el0 
     el12%coord=el0%coord+units*wdir*sdir*hh*0.5d0 
! before running differential corrections, tests that orbit is still good
     IF(bizarre(el12,e12))THEN 
        fail=.true. 
        IF(iun.gt.0)WRITE(iun,*)'int_step: bizarre ',e12 
        el1=el12
        RETURN 
     ENDIF
! fixed point iteration                                                 
     DO it=1,itx 
! compute vectorfield at intermediate point                             
        itmaxold=itmax 
        itmax=0 
        CALL whicor(0,icor,ncor,inew) 
! compute covariance and residual norm                                  
!       write(*,*)mc,iobs(1),ioco(1)                                 
        CALL diff_cor(mc,obs,obsw,el12,icor,iun,el12p,un12,csino12,delno12,succ)
        itmax=itmaxold 
! compute weak direction and length                                     
        CALL weak_dir(un12%g,wdir12,sdir12,iun,el12p%coo,el12p%coord,units) 
        direc=DOT_PRODUCT(wdir12,wdir)  
        IF(direc.lt.0.d0) wdir12=-wdir12
! store vectorfield at current iteration                               
        wdirst=wdir12
        sdirst=sdir12 
! recompute intermediate point  
        el12%coord=el0%coord+units*wdir12*sdir12*hh*0.5d0   
!    write(*,*)'weak sigma, iter, eq12',sdir12,it,(eq12(j),j=1,3) 
! before a new iteration, tests that orbit is still elliptic            
        IF(bizarre(el12,e12))THEN 
           IF(iun.gt.0)WRITE(iun,*)'int_step: bizarre ',e12 
           fail=.true. 
           el1=el12
           RETURN 
        ENDIF
! convergence control                                                   
        cosa=DOT_PRODUCT(wdir12,wdirst) ! cosa=prscag(6,wdir12,wdirst) 
        IF(cosa.le.1.d0)THEN 
           ang=acos(cosa)
        ELSE 
           ang=0.d0 
        ENDIF
        dsi=abs(sdirst-sdir12)/sdirst 
!    WRITE(*,*)' int_step: it, ang, dsi, sdir12,sdirst ',          
!    +          it,ang,dsi,sdir12,sdirst                                
        IF(ang.lt.eps.and.dsi.lt.eps)THEN 
!             WRITE(*,*)' int_step: it, ang, dsir ',it,ang,dsi           
! convergence achieved                                                  
           GOTO 2 
        ENDIF
        wdirst=wdir12
        sdirst=sdir12 
     ENDDO
! convergence failed                                                    
! WARNING: very verbose
     WRITE(*,*)' int_step: failed convergence, it, ang, dsig'         
     WRITE(*,*)it-1,ang,dsi 
! WARNING                      
     fail=.true. 
     IF(iun.gt.0)WRITE(iun,*)'int_step: non convergent ',el12 
     el1=el12 
     RETURN 
! convergence succeeded                                                 
2    CONTINUE 
! compute final point 
     el1=el0                                                          
     el1%coord=el0%coord+units*wdir12*sdir12*hh                              
! tests that orbit is still acceptable 
     IF(bizarre(el1,e1))THEN 
        fail=.true. 
        IF(iun.gt.0)WRITE(iun,*)'int_step: bizarre ',e1 
     ELSE 
        fail=.false. 
     ENDIF
! but we go on                                                          
  ELSE 
     WRITE(*,*)' int_step: not invented yet, imint=',imint 
     STOP 
  ENDIF
!     write(*,*)'int_step: at exit ',fail, (eq1(j),j=1,3)                
END SUBROUTINE int_step

! =================lov_int.mod===================================
! Created by Giacomo Tommei, 25 November 2002    
! Last changes, 17 December 2002                       
! ==================================================================
! LIST OF PUBLIC ENTITIES:
!
!  SUBROUTINES:                                                          
!                 lovinit                                                    
!                 lovobs                                                     
!                 lovmagn                                                    
!                 lovinterp: 
! ==================================================================
! ======================================================
!  lovinit 
! initializes storage of multiple solution data                 
! ======================================================
SUBROUTINE lovinit(astname,mulsodir,obsdir,progna,iunout,        &
       &     succ,vel_inf,imult,sigmax)
  USE obssto
  USE least_squares 
  USE tp_trace, ONLY: imul0                      
! ==================INPUT==============================  
  CHARACTER(LEN=9), INTENT(IN) :: astname           ! asteroid name 
  CHARACTER(LEN=80), INTENT(IN) ::  obsdir,mulsodir ! obs and mult. sol. 
                                                      ! directory  name 
  INTEGER, INTENT(IN) ::  iunout                    ! output file 
  CHARACTER(LEN=*), INTENT(IN) ::  progna           ! program name
!==================OUTPUT=============================
  LOGICAL, INTENT(OUT) :: succ             ! success flag 
  DOUBLE PRECISION, INTENT(OUT) :: vel_inf ! velocity at infinite 
                                             ! with respect to Earth (circular approx)
  INTEGER, INTENT(OUT) :: imult            ! multiple solution specifications       
  DOUBLE PRECISION, INTENT(OUT) :: sigmax 
! ===============OUTPUT THROUGH MODULE MULTIPLE_SOL===============
! time of initial conditions 
! index range of multiple solution, of nominal solution  
! =============for call to read_elems=============================
  CHARACTER(LEN=160) :: catname,repname,file   ! names 
  INTEGER iuncat 
  CHARACTER(LEN=19) :: name0 
  CHARACTER(LEN=9) :: name1 
  INTEGER :: le,iunrep, imtmp, iconv, norb
  DOUBLE PRECISION hmagn,tcat
  TYPE(orbit_elem) :: eltmp
  TYPE(orb_uncert) :: unctmp
  INTEGER :: j                            ! VA number 
  LOGICAL :: eof                          ! end of file  
  DOUBLE PRECISION :: v_infty,v_max,v_min ! velocity w.r.to Earth max and min
  CHARACTER(LEN=60) :: rwofi0             ! residuals and weights file name   
  LOGICAL :: obs0!,precob,change           ! logical flag et al
  CHARACTER(LEN=20) :: error_model ! error model file 
  INTEGER :: iun20                        ! unit 
  INTEGER :: i,ii                       ! loop indexes  
  INCLUDE 'sysdep.h90'                    ! system dependencies      
!===========================================================
  v_min=100.d0 
  v_max=0.d0 
  imul=0 ! flag to find if the elements have been read
! opening and reading multiple solution catalog                         
  succ=.false. 
  catname=mulsodir//dircha//astname//'.ctc' 
  CALL rmsp(catname,le)
  CALL filopn(iuncat,catname(1:le),'old')
  file=' '
  CALL oporbf(file,iuncat) 
! multiple solution report file                                         
  repname=mulsodir//dircha//astname//'.mrep' 
  CALL rmsp(repname,le) 
  CALL filopn(iunrep,repname(1:le),'old') 
  READ(iunrep,199)imult,sigmax 
199 FORMAT(34x,i4,10x,f5.2) 
  WRITE(iunout,200)imult, sigmax
200 FORMAT(' imult, sigma x ',i4,10x,f5.2) 
  READ(iunrep,299)scaling_lov,second_lov
299 FORMAT(12x,L1/12x,L1)
  WRITE(iunout,300)scaling_lov,second_lov
300 FORMAT(' scaling_lov ',L1,' second_lov ',L1)
  READ(iunrep,*) 
  WRITE(*,*) ' imult, sigmax from neomult ', imult,sigmax 
  imi0=imult+1 
  imul0=imi0 ! align tp_trace with multiple_sol
! read orbit one by one                                                 
  DO i=1,mulx 
     CALL read_elems(eltmp,name0,eof,file,iuncat,unctmp)
     IF(eof) GOTO 2 
! check consistency of epoch times
     IF(i.eq.1)THEN 
        tcat=eltmp%t 
     ELSE 
        IF(eltmp%t.ne.tcat)THEN 
           WRITE(*,*)'lovinit: time discrepancy from tcat=',tcat 
           WRITE(*,*)'lovinit: at record ',i,' t=',eltmp%t 
           STOP 
        ENDIF
     ENDIF
! control on elements type NO!!! 
! availability of covariance                                            
     IF (.NOT.(unctmp%succ)) THEN 
        WRITE(*,*)'lovinit:',  name0, ' matrices not avalaible' 
        STOP
     ENDIF
! handling of name; but note that the multiple solutions are assumed to 
! and with consecutive indexes!!! that is imul(i)=imul(1)+i-1           
     CALL splinam(name0,name1,j) 
     IF(j.gt.0.and.j.le.mulx)THEN
        imul(i)=j
     ELSEIF(j.lt.0)THEN
        WRITE(*,*)'lovinit: VA with negative index', j
        STOP
     ELSE
        WRITE(*,*)'lovinit: VA with index too large=',j,' increase mulx=',mulx
        STOP
     ENDIF
! note that the multiple solutions are assumed to be written
! with consecutive indexes!!! 
     IF(i.gt.1)THEN 
        IF(imul(i-1).eq.0)THEN
           WRITE(*,*)'lovinit: indexes not in order: VA ',j,' at record ',i,' missing ',j-1
        ENDIF
     ENDIF
! ok, copy in multiple_sol.mod
     elm(j)=eltmp
     unm(j)=unctmp
     READ(iunrep,145)ii,csinom(j),delnom(j),hmagn                   &
            &                    ,moid_m(j),dnp_m(j),dnm_m(j),iconv,sigq(j)
145  FORMAT(i5,1x,1p,e13.5,e11.3,2x,0p,f5.2,1x,                     &
            &              f8.5,1x,f10.5,1x,f10.5,1x,i2,1x,f6.3)
     IF(ii.ne.j)WRITE(*,*)'lovinit: discrepancy between mrep and ctc files'
! v_infinity computation                                                
     v_inf(j)=v_infty(elm(j)) 
     v_min=min(v_min,v_inf(j)) 
     v_max=max(v_max,v_inf(j)) 
! nominal solution                                                      
!     IF(j.eq.imi0)csinor=csinom(j) 
  ENDDO
! increase mulx                                                        
  WRITE(*,*)'lovinit: increase nmax, file is longer than',mulx 
2   CONTINUE 
  CALL clorbf 
  CALL filclo(iunrep,' ') 
  norb=i-1 
  imim=imul(1) 
  imip=imul(norb) 
  WRITE(iunout,*)' lovinit: input of ',norb,' multiple solutions' 
  WRITE(iunout,*)' lovinit: max of v_infty ',v_max,' minimum ',v_min 
  vel_inf=v_min 
  tdt_cat=tcat
! input observation data                                                
!    precob=.false. 
  iun20=iunout 
! to get error_model from rwo file; does not work with .obs file
  CALL retinobs(obsdir,astname,obs0,error_model,rms_m,rmsmag_m)
  CALL errmod_set(error_model)
! forces the error_model, but can use .obs file
  IF(m_m.gt.0.and.obs0)succ=.true. 
END SUBROUTINE lovinit
! =====================================================                  
! LOVMAGN                                                               
! provides magnitude and v_inty for the VA nearest to the               
! real index x; also velocity at infinity                       
! =====================================================                  
SUBROUTINE lovmagn(x,v_i,h) 
! ================INPUT===============================
  DOUBLE PRECISION, INTENT(IN) :: x       ! fractional index 
! ================OUTPUT==============================           
  DOUBLE PRECISION, INTENT(OUT) :: v_i    ! velocity at 
                                            ! infinity (U in au/day)
  DOUBLE PRECISION, INTENT(OUT) :: h      ! H magnitude 
! =============END INTERFACE===========================                 
  INTEGER :: j                            ! loops index
! ===================================================
  DO j=imim,imip 
     IF(x.lt.j)THEN 
        v_i=v_inf(j) 
        h=elm(j)%h_mag 
        GOTO 2 
     ENDIF
  ENDDO
  h=elm(imip)%h_mag 
  v_i=v_inf(imip) 
2 CONTINUE 
END SUBROUTINE lovmagn

! ===================================================                    
! LOVINTERP                                                             
! provides an interpolated orbit along the LOV                          
! (with non-integer index rindex)                                       
! ===================================================   
                 
SUBROUTINE lovinterp(rindex,deltasig,el0,unc0,succ)
  USE obssto
  USE least_squares
! ====================INPUT=============================      
  DOUBLE PRECISION, INTENT(IN) :: deltasig      ! interval in sigma 
                                                  ! between solutions
  DOUBLE PRECISION, INTENT(IN) ::  rindex       ! real index
! ====================OUTPUT============================ 
  TYPE(orbit_elem), INTENT(OUT) :: el0          ! elements     
  TYPE(orb_uncert), INTENT(OUT) :: unc0  ! normal and covariance matrices
                                                  ! corresponding to rindex
  LOGICAL, INTENT(OUT) :: succ                  ! success flag  
! ==================END INTERFACE=============================
  LOGICAL :: batch                              ! batch control
! ================LOV INTERPOLATION===========================
  DOUBLE PRECISION :: wdir(6),units(6),sdir,h0,diff(6) 
  LOGICAL :: fail, bizarre 
!  INTEGER :: iun20 
  DOUBLE PRECISION :: csinew,delnew,rmshnew    ! for constr_fit   
  INTEGER :: nused 
  TYPE(orbit_elem) :: elc                ! corrected 
! =============LOOP INDEXES, LOV INDEXES=======================
  INTEGER :: j,i,imu,ir 
  DOUBLE PRECISION :: s, ecc 
! =============================================================
! failure case ready                                                    
  succ=.false. 
! interpolate elements                                                  
  ir=NINT(rindex)
! extrapolation?                                                        
! WARNING: will result in failure if extrapolation is by a              
! long step beyond the end of the LOV string already computed           
  IF(ir.lt.imim)THEN 
     WRITE(*,*)'lovinterp: extrapolation rindex=',rindex,           &
            &        ' beyond ',imim                                        
     imu=imim 
  ELSEIF(ir.gt.imip)THEN 
     WRITE(*,*)'lovinterp: extrapolation rindex=',rindex,           &
            &        ' beyond ',imip                                     
     imu=imip 
  ELSE 
     imu=ir 
  ENDIF
! linear interpolation between the two consecutive ones                 
  s=rindex-imu 
!  iun20=-1  
  CALL weak_dir(unm(imu)%g,wdir,sdir,-1,elm(imu)%coo,elm(imu)%coord,units) 
! anti reversal check
  IF(imu.lt.imip)THEN
     diff=elm(imu+1)%coord-elm(imu)%coord
  ELSE
     diff=elm(imu)%coord-elm(imu-1)%coord
  ENDIF
  IF(DOT_PRODUCT(diff,wdir).lt.0.d0)wdir=-wdir
! nonlinear interpolation   
  CALL prop_sig(batch,elm(imu),el0,s,deltasig,m_m,obs_m,obsw_m,wdir,sdir,units,fail)
! check for hyperbolic                                                  
  IF(fail)THEN 
     WRITE(*,*)'step ',rindex,' hyperbolic' 
     WRITE(*,*)el0%coo,el0%coord 
     RETURN 
  ELSEIF(bizarre(el0,ecc))THEN 
     WRITE(*,*)'step ',rindex,' byzarre' 
     WRITE(*,*) el0%coo,el0%coord, ecc
     RETURN 
  ENDIF
! constrained corrections: 
  CALL constr_fit(m_m,obs_m,obsw_m,el0,wdir,elc,unc0,csinew,delnew,rmshnew,nused,succ)
! exit if not convergent                                                
  IF(.not.succ) THEN 
     WRITE(*,*)'lovinterp: const_fit failed for ',rindex 
     RETURN 
  ELSE
     el0=elc 
     succ=.true. 
  ENDIF
END SUBROUTINE lovinterp

! ====================================================                  
! LOVOBS                                                                
! ====================================================                  
! provides number and arc of observational data                         
SUBROUTINE lovobs(m0,nrej,calend1,calend2) 
  USE obssto
! ======================OUTPUT==========================
  INTEGER, INTENT(OUT) :: m0,nrej   ! number of obs., number rejected
  CHARACTER(LEN=14), INTENT(OUT) :: calend1,calend2 ! calendar dates of first 
                                                      ! and last
! =====================================================
!  INCLUDE 'parobx.h90'
  INTEGER :: j, ind(nobx)  ! loop index, sort index      
! =====================================================   
! calendar date
  CALL heapsort(obs_m(1:m_m)%time_tdt,m_m,ind)
  CALL calendwri(obs_m(ind(1))%time_tdt,calend1) 
  CALL calendwri(obs_m(ind(m_m))%time_tdt,calend2) 
  m0=m_m 
  nrej=0 
  DO j=1,m_m 
     IF(obsw_m(j)%sel_coord.eq.0)THEN
         nrej=nrej+1
     ENDIF
  ENDDO
END SUBROUTINE lovobs

END MODULE multiple_sol
! ====================================================================  
! Graham- Schmidt procedure to generate an orthonormal basis v          
! starting from 1  n-vector a                                           
! The new basis must be such that the first  vector is a                
SUBROUTINE graha_1(a,n,v) 
  implicit none 
  integer, intent(in) ::  n 
  double precision, intent(in) ::  a(n)
  double precision, intent(out) :: v(n,n)
! end interface 
  integer j,jok,jj,i 
  double precision cc,cc1,epsi,vl 
  integer,parameter :: nx=10  ! workspace
  double precision ws(nx) 
  logical ize 
! dimension check                                                       
  if(n.gt.nx)then 
     write(*,*)'n =',n,' larger than nx=',nx,' in graha' 
     stop 
  endif
! selection of the control for "zero" vectors                           
  cc1=sqrt(DOT_PRODUCT(a(1:n),a(1:n)))  ! cc1=sqrt(prscag(n,a,a)) 
  epsi=1.d-12*cc1 
  if(epsi.eq.0.d0)then 
     write(*,*)' a has rank zero' 
!        stop                                                           
  endif                                                                  
! V1 is the versor of A                                                 
  call versor(n,a,epsi,v(1,1),vl,ize) 
  if(ize)then 
     write(*,*)' first vector of a is too small' 
!        stop                                                           
  endif
! we now use the vectors of the canonic basis to supplement the span of 
  jok=0 
  do 1 j=1,n 
! remove the components along span(A), that is along V1                 
     cc1=-v(j,1) 
     do  i=1,n 
        ws(i)=cc1*v(i,1) 
     enddo
     ws(j)=ws(j)+1.d0 
     call versor(n,ws,epsi,v(1,2+jok),vl,ize) 
     if(.not.ize)then 
! now V(3+jok) is orthogonal to span(A); remove the components along    
! the previous ones (unless it is the first)                            
        if(jok.gt.0)then 
           do  jj=1,jok 
              cc=-DOT_PRODUCT(v(1:n,2+jok),v(1:n,1+jj)) 
                   ! cc=-prscag(n,v(1,2+jok),v(1,1+jj)) 
              v(1:n,2+jok)=v(1:n,2+jok)+cc*v(1:n,1+jj) 
                    ! call lincog(n,v(1,2+jok),1.d0,v(1,1+jj),cc,v(1,2+jok)) 
           enddo
           call versor(n,v(1,2+jok),epsi,v(1,2+jok),vl,ize) 
           if(ize)then 
              goto 1 
           endif
        endif
! the new versor is a good one                                          
        jok=jok+1 
        if(jok.eq.n-1)then 
           goto 2 
        endif
     endif
1 continue 
2 continue 
  if(jok.lt.n-1)then 
     write(*,*)' graha_1: something went wrong, jok=',jok 
  endif
END SUBROUTINE graha_1

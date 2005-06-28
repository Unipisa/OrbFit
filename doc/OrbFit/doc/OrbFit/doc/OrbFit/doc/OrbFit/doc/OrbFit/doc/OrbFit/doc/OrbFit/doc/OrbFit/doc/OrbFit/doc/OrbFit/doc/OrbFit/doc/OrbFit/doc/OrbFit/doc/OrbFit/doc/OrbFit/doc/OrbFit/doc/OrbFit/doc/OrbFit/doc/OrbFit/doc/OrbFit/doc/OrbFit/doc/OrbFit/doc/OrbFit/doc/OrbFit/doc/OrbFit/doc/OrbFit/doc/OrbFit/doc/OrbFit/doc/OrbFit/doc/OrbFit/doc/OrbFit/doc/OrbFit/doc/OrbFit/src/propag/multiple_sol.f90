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
PUBLIC :: fmuobs, fmupro, fmuplo, outmul, prop_sig
! former lov_int.mod
PUBLIC :: lovinit,lovobs,lovmagn,lovinterp
! PUBLIC data
! maximum number of multiple solutions
  INTEGER, PARAMETER::  mulx=4001
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
  DOUBLE PRECISION, DIMENSION(6) :: scales
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
  CALL weak_dir(uncert%g,wdir,sdir,iunint,el0%coo,el0%coord)
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
  IF(batch)THEN 
! preselect a  stepsize which could not result in hyperbolic orbit      
! hypothetical final point                                              
     ratio=sqrt(2.d0) 
     IF(scaling_lov)THEN
        CALL scale_coef(el0%coo,el0%coord,scales)
        scales=1.d0/scales
     ELSE
        scales=1.d0
     ENDIF
     DO k=1,10 
        eqf=el0
        DO j=1,6 
           eqf%coord(j)=el0%coord(j)+wdir(j)*scales(j)*sdir*dn*sigma*imult 
        ENDDO
        IF(bizarre(eqf,ecc))THEN 
           dn=dn/ratio 
           IF(verb_mul.gt.9)WRITE(iun,*)' f_multi: shortening to avoid ecc= ',ecc,k,dn
           IF(verb_mul.gt.19) WRITE(iun,*)eqf               
        ELSE 
           GOTO 9 
        ENDIF
     ENDDO
9    CONTINUE 
  ENDIF
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
     CALL prop_sig(batch,elm(imi-1),elc,dn,sigma,mc,obs,obsw,wdir,sdir,fail)
! check for hyperbolic                                                  
     IF(fail)THEN 
        IF(.not.batch)WRITE(*,*)'step ',imi,' failed' 
        IF(.not.batch)WRITE(*,*)elc 
        IF(verb_mul.gt.9)WRITE(iun,*)'fail in prop_sig, imi= ',imi,dn
        imi=imi-1
        GOTO 6 
     ELSEIF(bizarre(elc,ecc))THEN 
        IF(.not.batch)WRITE(*,*)'step ',imi,' bizarre', 'ecc=', ecc 
        IF(.not.batch)WRITE(*,*)elc 
        IF(verb_mul.gt.9)WRITE(iun,*)'bizarre out of  prop_sig, imi= ',imi,ecc,dn
        imi=imi-1 
        GOTO 6 
     ELSE 
! constrained differential corrections:
         CALL constr_fit(mc,obs,obsw,elc,wdir,elm(imi),unm(imi),     &
     &       csinom(imi),delnom(imi),rmsh,nused,succ)
! exit if not convergent                                                
        IF(.not.succ) THEN 
           IF(verb_mul.gt.19)WRITE(iun,*)'fail in constr_fit, imi= ',imi,dn
           imi=imi-1 
           GOTO 6 
        ENDIF
! orbital distance                                                      
        CALL nomoid(elm(imi)%t,elm(imi),moid_m(imi),dnp_m(imi),dnm_m(imi))
        iconv(imi)=0
! check for sigQ.le.sigma                                               
        sigq(imi)=sqrt(abs(csinom(imi)**2-csinom(imi0)**2)*nused)/   &
     &             rescov(6,nused,csinor)  
        IF(csinom(imi).lt.csinom(imi0))sigq(imi)=-sigq(imi)
        IF(batch.and.sigq(imi).gt.sigma)THEN 
           IF(verb_mul.gt.19)WRITE(iun,*)'too large sigma_q ',imi,sigq(imi)
           IF(sigq(imi).gt.sigma*1.1d0) imi=imi-1
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
     CALL prop_sig(batch,elm(imi+1),elc,dn,sigmam,mc,obs,obsw,wdir,sdir,fail) 
! check for hyperbolic                                                  
     IF(fail)THEN 
        IF(.not.batch)WRITE(*,*)'step ',imi,' failed' 
        IF(.not.batch)WRITE(*,*)elc 
        IF(verb_mul.gt.19)WRITE(iun,*)'fail in prop_sig, imi= ',imi,dn
        imi=imi+1 
        GOTO 8 
     ELSEIF(bizarre(elc,ecc))THEN 
        IF(.not.batch)WRITE(*,*)'step ',imi,' bizarre' 
        IF(.not.batch)WRITE(*,*)elc
        IF(verb_mul.gt.19)WRITE(iun,*)'bizarre out of prop_sig, imi= ',imi,ecc,dn
        imi=imi+1 
        GOTO 8 
     ELSE 
! differential corrections:                                             
! constrained differential corrections:
         CALL constr_fit(mc,obs,obsw,elc,wdir,elm(imi),unm(imi),     &
     &       csinom(imi),delnom(imi),rmsh,nused,succ)           
! exit if not convergent                                                
        IF(.not.succ)THEN 
           IF(verb_mul.gt.19)WRITE(iun,*)'fail in constr_fit, imi= ',imi,dn
           imi=imi+1 
           GOTO 8 
        ENDIF
! orbital distance                                                      
        CALL nomoid(elm(imi)%t,elm(imi),moid_m(imi),dnp_m(imi),dnm_m(imi))
        iconv(imi)=0
        sigq(imi)=sqrt(abs(csinom(imi)**2-csinom(imi0)**2)*nused)/   &
     &             rescov(6,nused,csinor)  
        IF(csinom(imi).lt.csinom(imi0))sigq(imi)=-sigq(imi)
        IF(batch.and.sigq(imi).gt.sigma)THEN 
           IF(verb_mul.gt.19)WRITE(iun,*)'too large sigma_q ',imi,sigq(imi)
           IF(sigq(imi).gt.sigma*1.1d0) imi=imi+1
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
  CALL tee(iun,'no       a      h      k      p      q      lambda=')          
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
  INTEGER iconv(mulx) ! obsolete moid flag 
! catalog obj name (with underscore)                                    
  CHARACTER*20 astna0 
  INTEGER le,nscal,j
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
     CALL nomoid(elc%t,elm(i),moid_m(i),dnp_m(i),dnm_m(i)) 
     iconv(i)=0
     moid_min=min(moid_m(i),moid_min) 
! warning: assumption that scalar obs=2*no.obs may be invalid for radar
     nscal=0
     DO j=1,mc
        IF(obs(j)%type.eq.'O'.or.obs(j)%type.eq.'S'.and.obsw(j)%sel_coord.gt.0)THEN
           nscal=nscal+2
        ELSEIF(obs(j)%type.eq.'R'.or.obs(j)%type.eq.'V'.and.obsw(j)%sel_coord.gt.0)THEN
           nscal=nscal+1
        ENDIF
     ENDDO
     IF(csinom(i).gt.csinor)THEN 
        sigq(i)=sqrt((csinom(i)**2-csinor**2)*mc*2)                  &
     &          /rescov(6,nscal,csinor)                                  
     ELSE 
        sigq(i)=0.d0 
     ENDIF
     WRITE(iunrep,145)i,csinom(i),delnom(i),elm(i)%h_mag             &
     &                    ,moid_m(i),dnp_m(i),dnm_m(i),iconv(i),sigq(i)       
145  FORMAT(i5,1x,1p,e13.5,e11.3,2x,0p,f5.2,1x,                      &
     &              f8.5,1x,f10.5,1x,f10.5,1x,i2,1x,f7.3)                 
     astna0=' '
     WRITE(astna0,111)i 
111  FORMAT('mult_',i4) 
     CALL rmsp(astna0,le) 
! output multiline catalog                                              
     CALL write_elems(elm(i),astna0(1:le),'ML',' ',iunctc,unm(i))
! output one line catalog, if possible
     CALL coo_cha(elm(i),'KEP',elk,fail_flag)
     IF(fail_flag.ge.4)THEN
        WRITE(*,*)' hyperbolic orbit in output ', fail_flag, i     
     ELSE
        CALL write_elems(elk,astna0(1:le),'1L',' ',iuncat)
     ENDIf
  ENDDO
  CALL filclo(iunctc,' ') 
  CALL filclo(iunrep,' ') 
  CALL filclo(iuncat,' ') 
!=============================================                          
END SUBROUTINE nmulti
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
  CALL rmsp(catname,le) 
  INQUIRE(file=catname(1:le),exist=ok) 
  IF(.not.ok)THEN 
     WRITE(*,*)' file ',catname(1:le),' not found' 
     RETURN 
  ENDIF
  CALL oporbf(catname(1:le),0) 
  DO i=1,mulx 
     CALL rdorb(name0,eq,eltype,t,g,defcov,                         &
     &        c,defnor,h,sl,mass,rsys,epoch,no,eof)                     
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
! control on elements type                                              
     IF(eltype.ne.'EQU')THEN 
        WRITE(*,*)'mult_input: non equinoctal, but of type ',eltype 
!        RETURN 
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
  CALL tee(iun20,                                                   &
     &  'no.,     a      h      k      p      q      lambda=')          
  DO i=imim,imip 
     WRITE(*,144)i,elm(i)%coord 
144  FORMAT(i3,6f12.8) 
     WRITE(iun20,144)i,elm(i)%coord
  ENDDO
  CALL tee(iun20,'no.,  magn,  MOID ,  nod+  ,  nod-=') 
  DO i=imim,imip 
     WRITE(*,145)i,(i),moid_m(i),dnp_m(i),dnm_m(i),iconv(i) 
     WRITE(iun20,145)i,elm(i)%h_mag,moid_m(i),dnp_m(i),dnm_m(i),iconv(i) 
145  FORMAT(i3,2x,f5.2,1x,f8.5,1x,f8.5,1x,f8.5,1x,i2) 
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
SUBROUTINE prop_sig(batch,el1,el2,dn,sigma,mc,obs,obsw,wdir,sdir,fail)
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
  DOUBLE PRECISION, DIMENSION(6), INTENT(INOUT) :: wdir
  DOUBLE PRECISION, INTENT(IN) :: sdir
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
!  WRITE(*,*)'mmulti: ch,h, eqt ', ch,h,(eqt(j),j=1,3) 
!  WRITE(*,*)' prop-sig: about to enter int_step',tc, obs(mc)%time_tdt
  CALL int_step(elt,el2,h,imint,                                     &
     &        mc,obs,obsw,iun,wdir,sdir,fail) 
! write(*,*)fail,eq2                                                
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
     CALL diff_cor(mc,obs,obsw,elt,icor,iun,              &
     &        eltp,uncert,csinor,delnor,succ)                        
     itmax=itmaxold 
! compute weak direction and length                                     
     CALL weak_dir(uncert%g,wdir,sdir,iun,eltp%coo,eltp%coord) 
     IF(DOT_PRODUCT(wdir,wdir0).lt.0.d0)wdir=-wdir
     GOTO 1 
  ENDIF
  IF(iun.gt.0.and.verb_mul.gt.20)WRITE(iun,*)' prop_sig: stepsize ',h 
END SUBROUTINE prop_sig
!=================================================                      
! INT_STEP integration step for propagation along LOV                    
SUBROUTINE int_step(el0,el1,hh,imint,                              &
     &       mc,obs,obsw,iun,wdir,sdir,fail)                            
! stepsize, integration method                                          
  DOUBLE PRECISION, INTENT(IN) ::  hh 
  INTEGER, INTENT(IN) :: imint 
! weak direction (in input computed at eq0)             
  DOUBLE PRECISION, INTENT(IN) :: wdir(6)
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
  DOUBLE PRECISION, DIMENSION(6)   ::  scales
! =========end interface================                                
! intermediate point: state, covariance, normal matr., RMS, norm corr   
  TYPE(orbit_elem) ::  el12,el12p
  TYPE(orb_uncert) :: un12
  DOUBLE PRECISION :: csino12,delno12 
! weak direction (in input computed at eq0), length of axis             
  DOUBLE PRECISION sdir,wdir12(6),sdir12,wdirst(6),sdirst 
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
  IF(scaling_lov)THEN
     CALL scale_coef(el0%coo,el0%coord,scales)
     scales=1.d0/scales
  ELSE
     scales=1.d0
  ENDIF 
  IF(imint.eq.1)THEN 
! Euler method
     el1=el0                                                          
     el1%coord=el0%coord+scales*wdir*sdir*hh 
     fail=.false. 
  ELSEIF(imint.eq.2)THEN 
! Runge-Kutta-Gauss of order 2:                                         
     itx=5 
     eps=1.d-4 
! first store vectorfield as it is on eq0                               
     wdirst=wdir
     sdirst=sdir 
! compute intermediate point (first guess)                              
     el12=el0 
     el12%coord=el0%coord+scales*wdir*sdir*hh*0.5d0 
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
        CALL diff_cor(mc,obs,obsw,el12,icor,iun,   &
     &     el12p,un12,csino12,delno12,succ)                    
        itmax=itmaxold 
! compute weak direction and length                                     
        CALL weak_dir(un12%g,wdir12,sdir12,iun,el12p%coo,el12p%coord) 
        direc=DOT_PRODUCT(wdir12,wdir)  
        IF(direc.lt.0.d0)   wdir12=-wdir12
! recompute intermediate point  
        el12%coord=el0%coord+scales*wdir*sdir*hh*0.5d0   
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
        wdirst=wdir12 ! CALL vcopy(6,wdir12,wdirst) 
        sdirst=sdir12 
     ENDDO
! convergence failed                                                    
!       WRITE(*,*)' int_step: failed convergence, it, ang, dsig'         
!       WRITE(*,*)it-1,ang,(sdir12-sdirst)/sdirst                       
     fail=.true. 
     IF(iun.gt.0)WRITE(iun,*)'int_step: non convergent ',el12 
     el1=el12 
     RETURN 
! convergence succeeded                                                 
2    CONTINUE 
! compute final point 
     el1=el0                                                          
     el1%coord=el0%coord+scales*wdir12*sdir*hh                              
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
        imul(j)=j
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
        IF(imul(j-1).eq.0)THEN
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
  DOUBLE PRECISION :: wdir(6),sdir,h0,diff(6) 
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
  ir=rindex
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
  CALL weak_dir(unm(imu)%g,wdir,sdir,-1,elm(imu)%coo,elm(imu)%coord) 
! anti reversal check
  IF(imu.lt.imip)THEN
     diff=elm(imu+1)%coord-elm(imu)%coord
  ELSE
     diff=elm(imu)%coord-elm(imu-1)%coord
  ENDIF
  IF(DOT_PRODUCT(diff,wdir).lt.0.d0)wdir=-wdir
! nonlinear interpolation   
  CALL prop_sig(batch,elm(imu),el0,s,deltasig,m_m,obs_m,obsw_m,wdir,sdir,fail)        
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
     WRITE(*,*)'lovinterp: diff_vin failed for ',rindex 
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

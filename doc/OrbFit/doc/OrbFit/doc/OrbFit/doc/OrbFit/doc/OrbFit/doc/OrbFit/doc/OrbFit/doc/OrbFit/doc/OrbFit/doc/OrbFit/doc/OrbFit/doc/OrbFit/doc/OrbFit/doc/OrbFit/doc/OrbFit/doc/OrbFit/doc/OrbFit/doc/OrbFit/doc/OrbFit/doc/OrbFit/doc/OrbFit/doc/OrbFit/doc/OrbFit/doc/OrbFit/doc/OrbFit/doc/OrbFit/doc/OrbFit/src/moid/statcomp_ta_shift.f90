
! ************************************************************      
! ***********  C O M P U T I N G   M O I D   for   ***********
! ***********  C O M E T A R Y     O R B I T S     ***********
! ************************************************************      
! ****** written by GIOVANNI F. GRONCHI (October 2004) *******
! ************************************************************      
! last modified 24/06/2005, GFG
! ============================================================      
  PROGRAM statcomp_ta_shift
    USE critical_points_shift
    USE output_control
    USE fund_const                                    
    USE orbit_elements
    USE planet_masses
!    USE dist_grad
    IMPLICIT NONE 
    TYPE(orbit_elem) :: eccom   ! cometary elements of the comet
    TYPE(orbit_elem) :: eqp     ! equinoctal elements of the Earth
    TYPE(orbit_elem) :: ecpl    ! cometary elements of the Earth
! cometary elements: (q=a(1-e), e, I, Omnod, omega, t_0)
    REAL(KIND=8) :: t1 ! time for earth call
    REAL(KIND=8),DIMENSION(6) :: eqptmp ! auxiliary
!    INTEGER,PARAMETER :: poldeg=16
    REAL(KIND=8) :: fpl(poldeg),fcom(poldeg) ! true anomalies 
    REAL(KIND=8) :: D2(poldeg) ! squared distance function 
                               ! at critical points 
    INTEGER :: D2_ord(poldeg) ! sorted indexes
    INTEGER nstat ! number of stationary points (maximum = 16)   
    INTEGER nummin,nummax ! number of relative minima/maxima found
! --------- error flags ---------
    INTEGER :: answer(poldeg) ! type of singular point : 
!                             answer(j) =  1    MAXIMUM
!                             answer(j) = -1    MINIMUM                
!                             answer(j) =  0    SADDLE                 
!                             answer(j) =  2    CANNOT DECIDE          
    CHARACTER*7 :: singtype 

    LOGICAL :: morse ! check with Morse theory 
!                        morse = true      OK   
!                        morse = false     ERROR  
    LOGICAL :: weier ! check with Weierstrass theory:
!                        weier = true      OK 
!                        weier = false     ERROR
    LOGICAL :: warnflag(3) ! program warning flags:
!      warnflag(1) = true    OK
!      warnflag(1) = false   leading coefficient of resultant is very small
!      warnflag(2) = true    OK                                          
!      warnflag(2) = false   higher degree terms in resultant are not small
!      warnflag(3) = true    OK                                           
!      warnflag(3) = false   low precision for resultant coefficients     
    LOGICAL :: sflag(6) ! solving system messages:
!         sflag(1) = true    OK
!         sflag(1) = false   there are two good solutions             
!         sflag(2) = true    OK                                               
!         sflag(2) = false   neither of the two evaluation is close to 0
!         sflag(3) = true    OK                                               
!         sflag(3) = false   the s-component of the solution is complex
!         sflag(4) = true    OK    
!         sflag(4) = false   leading coeff. of the 2nd degree pol. is small
!         sflag(5) = true    OK    
!         sflag(5) = false   1st and 2nd coeffs except of 2nd degree poly
!                            are small                       
!         sflag(6) = true    OK    
!         sflag(6) = false   2nd degree poly have all coefficients small
    LOGICAL :: hzflag ! hzflag = .true.  OK!
                      ! hzflag = .false. abs(root)>10^5
                      ! (large root of the resultant poly)
    LOGICAL :: hwflag ! hwflag = .true.  OK!
                      ! hwflag = .false. abs(root)>10^5
                      ! (large root of 2nd deg poly)
    LOGICAL :: multfl ! multfl = .true.  OK!
                      ! multfl = .false. 0 has multiplicity > 4
    LOGICAL :: hevalflag ! hevalflag = .true. OK!
                         ! hevalflag = .false. unsuccessful hessian evaluation
! --------------------------------------------------------------------
    LOGICAL :: wri_mat ! wri_mat=.true.   write file for matlab
                       ! wri_mat=.false.  do not write file for matlab
    LOGICAL :: chk_ecc_anomalies ! chk_ecc_anomalies=.true.  perform the check
                                 ! chk_ecc_anomalies=.false. do nothing
!   auxiliary                                               
    REAL(KIND=8) :: ffpl,ffcom,chfpl,chfcom
    REAL(KIND=8) :: DD2 
!   auxiliary                                               
    REAL(KIND=8) :: cosupl,sinupl,cosu,sinu,uupl,uu,chupl,chu
! ------------------------------------------
    INTEGER iundat,iunsta
    INTEGER j  ! loop index     
!
    INCLUDE 'jplhdr.h90' ! to compute mu of the Earth
!
    INTEGER fail_flag
    character*6 progna 
    character*80 run 
! ==================================================================
 
    progna='fitobs'
    run='moid_comets'
! read options
    CALL trivopt(progna,run,iun_log)
    CALL rmodel
    CALL trange

! flags initialization                                              
    weier = .true.
!    weier = .false. ! dummy
    morse = .true.
    sflag(1) = .true. 
    sflag(2) = .true. 
    sflag(3) = .true.
    sflag(4) = .true.
    sflag(5) = .true.
    sflag(6) = .true.
    warnflag(1) = .true.
    warnflag(2) = .true.
    warnflag(3) = .true.
    hzflag = .true.
    hwflag = .true.
    multfl = .true.
! verbosity control
    verb_moid=1 ! no verbosity
!    verb_moid = 21 ! high verbosity
! write file for matlab
    wri_mat =.true. ! write file
!    wri_mat =.false. ! do not write
! check eccentric anomalies (only for elliptic orbits)
    chk_ecc_anomalies = .false.
!    chk_ecc_anomalies = .true.

! Earth elements
    t1 = 53450.d0
    CALL earth(t1,eqptmp(1:6))
    eqp=undefined_orbit_elem
    eqp%coord(1:6)=eqptmp(1:6)
    eqp%coo='EQU'
    eqp%t=t1
    CALL coo_cha(eqp,'COM',ecpl,fail_flag) !conversion into cometary elems
    write(*,*)'Earth elems (cometary):',ecpl%coord(1:2),ecpl%coord(3:5)*degrad

! comet elements
    CALL filopn(iundat,'comelems.dat','old') 
    READ(iundat,*)
    eccom=undefined_orbit_elem
    READ(iundat,*) eccom%coord(1:6)
    CALL filclo(iundat,' ') 
    eccom%coo='COM'
    eccom%t=eccom%coord(6) ! it is not used later
    eccom%coord(3:5) = radeg*eccom%coord(3:5)
    write(*,*)'comet elems (cometary):',eccom%coord(1:2), &
         & eccom%coord(3:5)*degrad

! dummy
!    ecpl%coord(2) = 0.d0

    CALL compute_critical_points_shift(ecpl%coord(1:5),eccom%coord(1:5), &
         & fpl,fcom,nstat,nummin,nummax,answer,warnflag,sflag,morse, &
         & weier,hzflag,hwflag,multfl,hevalflag)
!    write(*,*)'nstat',nstat
!    write(*,*)'critical points',fcom(1:nstat),fpl(1:nstat)
!   check
    if(.not.sflag(3)) then
       write(*,*)'stat_comp: negative discriminant!'
    endif

! loop on number of stationary points                               
    DO j = 1,nstat 
! conversion in radians for D2eval                                  
       ffpl=fpl(j)*radeg 
       ffcom=fcom(j)*radeg 
! compute d^2 at critical (fpl,fcom)
       CALL d2eval_ta(ffpl,ffcom,DD2) 
       D2(j)=DD2 
! angles between 0 and 360 degrees                                  
       CALL choosedeg(fpl(j),chfpl) 
       CALL choosedeg(fcom(j),chfcom)        
       fpl(j) = chfpl 
       fcom(j) = chfcom 
    ENDDO

! sorting 3-uples (fcom,fpl,D2) according to value of D2               
    CALL heapsort(D2,nstat,D2_ord)

! stat_data file for matlab 
    if(wri_mat) then
       CALL filopn(iunsta,'stat_data','unknown') 
    endif

! ******************** WRITING OUTPUT ON SCREEN ********************
    WRITE(*,*)'#######################################################&
         &########'                                                         
    WRITE(*,*)'######### S T A T I O N A R Y   P O I N T S ###########&
         &########'                                                         
    WRITE(*,*)'#######################################################&
         &########'                                                         
    WRITE(*,*)'       V              v             DIST          TYP&
         &E'                                                                
    WRITE(*,*)'=======================================================&
         &========'                                                         
! loop on number of stationary points                               
    DO j = 1,nstat 
       write(*,*) 
       IF (answer(D2_ord(j)).eq.1) THEN 
          singtype='MAXIMUM' 
       ELSEIF (answer(D2_ord(j)).eq.-1) THEN 
          singtype='MINIMUM' 
       ELSEIF (answer(D2_ord(j)).eq.0) THEN 
          singtype='SADDLE' 
       ELSEIF (answer(D2_ord(j)).eq.-2) THEN 
          singtype='ERROR' 
       ENDIF
       WRITE(*,105)fpl(D2_ord(j)),fcom(D2_ord(j)),dsqrt(D2(D2_ord(j))), &
            & singtype 

! write on stat_data file for matlab 
       if(wri_mat) then
          WRITE(iunsta,106)nstat,nummin,answer(D2_ord(j)), &
               & fcom(D2_ord(j)),fpl(D2_ord(j)) 
       endif

    ENDDO

105 FORMAT(2x,f13.8,3x,f13.8,3x,f13.8,5x,a7) 
106 FORMAT(1x,i3,2x,i3,2x,i3,2x,f11.6,2x,f11.6) 

    IF(chk_ecc_anomalies) THEN
! ************************************************************************
! write the corresponding eccentric anomalies
       WRITE(*,*)'#######################################################&
            &########'                                                         
       WRITE(*,*)'######### S T A T I O N A R Y   P O I N T S ###########&
            &########'                                                         
       WRITE(*,*)'#######################################################&
            &########'                                                         
       WRITE(*,*)'       upl              u             DIST          TYP&
            &E'                                                                
       WRITE(*,*)'=======================================================&
            &========'                                                         
       DO j = 1,nstat 
          cosupl = (ecpl%coord(2) + cos(radeg*fpl(D2_ord(j))))/ &
               & (1.d0+ecpl%coord(2)*cos(radeg*fpl(D2_ord(j))))
          sinupl = sin(radeg*fpl(D2_ord(j)))*sqrt(1.d0+ecpl%coord(2))/&
               &(1.d0+ecpl%coord(2)*cos(radeg*fpl(D2_ord(j))))
          uupl = datan2(sinupl,cosupl)
          CALL choosedeg(degrad*uupl,chupl) 
          uupl=chupl
          cosu = (eccom%coord(2) + cos(radeg*fcom(D2_ord(j))))/ &
               & (1.d0+eccom%coord(2)*cos(radeg*fcom(D2_ord(j))))
          sinu = sin(radeg*fcom(D2_ord(j)))*sqrt(1.d0+eccom%coord(2))/&
               &(1.d0+eccom%coord(2)*cos(radeg*fcom(D2_ord(j))))
          uu = datan2(sinu,cosu)
          CALL choosedeg(degrad*uu,chu) 
          uu=chu
          WRITE(*,108)uupl,uu,dsqrt(D2(D2_ord(j))),answer(D2_ord(j))
       ENDDO
108 FORMAT(2x,f13.7,3x,f13.7,3x,f13.8,5x,i2) 

! ************************************************************************
    ELSE
! do nothing
    ENDIF

    if(wri_mat) then
       CALL filclo(iunsta,' ') 
    endif

  END PROGRAM statcomp_ta_shift

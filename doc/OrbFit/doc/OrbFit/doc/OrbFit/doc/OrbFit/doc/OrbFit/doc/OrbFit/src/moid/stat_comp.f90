! ************************************************************      
! ** C O M P U T I N G   S T A T I O N A R Y    P O I N T S **      
! ************************************************************      
! ********** written by GIOVANNI F. GRONCHI (2003) ***********      
! ******* Department of Mathematics, UNIVERSITY of PISA ******      
! ============================================================      
  PROGRAM stat_comp
    USE fund_const                                    
    IMPLICIT NONE 
!   elements of the Earth and of the asteroid                         
    DOUBLE PRECISION,DIMENSION(6) :: ekpl,elkep 
    DOUBLE PRECISION :: u(20),upl(20) ! eccentric anomalies 
    DOUBLE PRECISION :: D2(20) ! squared distance function 
                               ! at critical points 
    INTEGER nstat ! number of stationary points (maximum = 20)   
    INTEGER nummin,nummax ! number of relative minima/maxima found
    INTEGER :: answer(20) ! type of singular point : 
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
!   auxiliary                                               
    DOUBLE PRECISION uu,uupl,chu,chupl 
    DOUBLE PRECISION DD2 
!   mutual orbital elements 
    DOUBLE PRECISION :: mutI,mutom,mutompl 
    DOUBLE PRECISION :: a,e,i,om,apl,epl,ompl 
    DOUBLE PRECISION :: beta,betapl 
    COMMON/bbpl/beta,betapl 
!     ------------------------------------------
!     ordered 3-uples, for sorting
    DOUBLE PRECISION ordu(20),ordupl(20),ordD2(20) 
    INTEGER ordans(20) 
!     ------------------------------------------
    INTEGER iundat,iunmat,iunsta
    INTEGER j  ! loop index     
!     ==================================================================


!   flags initialization                                              
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

!   reading orbital elements
    CALL filopn(iundat,'orbelems.dat','old') 
    READ(iundat,*)
    READ(iundat,*) elkep
    READ(iundat,*)
    READ(iundat,*) ekpl
    CALL filclo(iundat,' ') 

    CALL filopn(iunmat,'orbit_data','unknown') 
    WRITE(iunmat,103) elkep(1:5)
    WRITE(iunmat,103) ekpl(1:5)

!   conversion into radians
    ekpl(3:5) = radeg*ekpl(3:5)
    elkep(3:5) = radeg*elkep(3:5)
!   switch to mutual elements
    CALL mutualrefcha(ekpl,elkep,mutI,mutom,mutompl)

    WRITE(iunmat,103) mutI*degrad,mutom*degrad,mutompl*degrad,0.d0,0.d0
    CALL filclo(iunmat,' ') 

103 FORMAT(f10.5,2x,f10.5,2x,f10.5,2x,f10.5,2x,f10.5)

    a    = elkep(1) 
    e    = elkep(2) 
    apl  = ekpl(1) 
    epl  = ekpl(2) 
    i    = mutI 
    om   = mutom 
    ompl = mutompl 
  
    beta=dsqrt(1-e**2) 
    betapl=dsqrt(1-epl**2) 

!   ======================================================
!   COMPUTING FUNCTIONS of the CONSTANT ORBITAL ELEMENTS              
!   ======================================================
    CALL aical(a,e,i,om,apl,epl,ompl) 

!     ==================================================================
!     COMPUTING STATIONARY POINTS FOR D2 AND DECIDING WHICH ARE         
!     RELATIVE MAXIMA/MINIMA BETWEEN THEM                               
    CALL comp_heart_rot(u,upl,nstat,nummin,nummax, &
         & answer,warnflag,sflag,morse,weier,hzflag,hwflag,multfl)
!     ==================================================================

!   check
    if(.not.sflag(3)) then
       write(*,*)'stat_comp: negative discriminant!  ecc=',e
    endif

!     loop on number of stationary points                               
    DO j = 1,nstat 

!     conversion in radians for D2eval                                  
       uu=u(j)*radeg 
       uupl=upl(j)*radeg 

!     ==================================================================
!     COMPUTE SQUARED DISTANCE in the points u,upl                      
       CALL D2eval(uu,uupl,DD2) 
!     ==================================================================

       D2(j)=DD2 

!     angles between 0 and 360 degrees                                  
       CALL choosedeg(u(j),chu) 
       CALL choosedeg(upl(j),chupl) 
       
       u(j) = chu 
       upl(j) = chupl 

!     END MAIN DO LOOP                                                  
    ENDDO

!     =====================================================             
!     sorting 3-uples (u,upl,D2) according to value of D2               
    CALL sortd2(nstat,u,upl,D2,answer) 
!     =====================================================             


!   stat_data file for matlab 
    CALL filopn(iunsta,'stat_data','unknown') 

!     ******************************************************************
!     ******************** WRITING OUTPUT ON SCREEN ********************
!     ******************************************************************
    WRITE(*,*)'#######################################################&
         &########'                                                         
    WRITE(*,*)'######### S T A T I O N A R Y   P O I N T S ###########&
         &########'                                                         
    WRITE(*,*)'#######################################################&
         &########'                                                         
    WRITE(*,*)'       u              upl             DIST          TYP&
         &E'                                                                
    WRITE(*,*)'=======================================================&
         &========'                                                         
!     loop on number of stationary points                               
    DO j = 1,nstat 
       write(*,*) 
       IF (answer(j).eq.1) THEN 
          singtype='MAXIMUM' 
       ELSEIF (answer(j).eq.-1) THEN 
          singtype='MINIMUM' 
       ELSEIF (answer(j).eq.0) THEN 
          singtype='SADDLE' 
       ELSEIF (answer(j).eq.-2) THEN 
          singtype='ERROR' 
       ENDIF
                                                                        
!     writing on screen                                                 
       WRITE(*,105)u(j),upl(j),dsqrt(D2(j)),singtype 

!   write on stat_data file for matlab 
       WRITE(iunsta,106)nstat,nummin,answer(j),u(j),upl(j) 
       
    ENDDO

105 FORMAT(2x,f13.8,3x,f13.8,3x,f13.8,5x,a7) 
106 FORMAT(1x,i3,2x,i3,2x,i3,2x,f11.6,2x,f11.6) 

    CALL filclo(iunsta,' ') 

    STOP
  END PROGRAM stat_comp

! ******************************************************************
! ********** COMPUTATION of COEFFICIENTS of SYLVESTER MATRIX *******
! ******************************************************************
! *********** written by GIOVANNI F. GRONCHI (2001) ****************
! ******** Department of Mathematics, UNIVERSITY of PISA ***********
! last modified June 2004 (GFG)
! ==================================================================
  SUBROUTINE aical(a,e,i,om,apl,epl,ompl) 
    IMPLICIT NONE 
! ============ orbital elements =================                   
    DOUBLE PRECISION,INTENT(IN) :: a,e,i,om,apl,epl,ompl 
! ------------------------------ end interface ---------------------
! coefficient of 20th degree (for check)                            
    DOUBLE PRECISION :: coe20 
! ============ functions of the elements ================           
    DOUBLE PRECISION Aa,Ba,Ca,Da,Ea,Fa 
    DOUBLE PRECISION Ap,Bp,Cp,Dp,Ep,Fp 
!                                                                       
    DOUBLE PRECISION A1,A2,A3,A4,A5,A6,A7 
    DOUBLE PRECISION A8,A9,A10,A11,A12,A13,A14,A15 
    COMMON/Aj1to15/ A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15 
! check coefficient of 20th degree                                  
    DOUBLE PRECISION tmpcoe 
!                                                                       
    DOUBLE PRECISION beta,betapl 
    COMMON/bbpl/beta,betapl 
! ==================================================================
! HINT: angles must be passed in radians                            
                                                                        
    Aa = a*cos(om) 
    Ba =-a*beta*sin(om) 
    Ca =-a*e*cos(om) 
    Da = a*sin(om) 
    Ea = a*beta*cos(om) 
    Fa =-a*e*sin(om) 
!                                                                       
    Ap = apl*cos(ompl) 
    Bp =-apl*betapl*sin(ompl) 
    Cp =-apl*epl*cos(ompl) 
    Dp = apl*sin(ompl) 
    Ep = apl*betapl*cos(ompl) 
    Fp =-apl*epl*sin(ompl) 
                                                                        
! ===================================================               
! A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14                    
! ===================================================               
    A1 = (a**2)*(beta**2) 
    A2 = 0 
    A3 = a**2 
    A4 = (apl**2)*(betapl**2) 
    A5 = 0 
    A6 = apl**2 
    A7 = -2.d0*a*apl*beta*betapl*(sin(om)*sin(ompl)+                  &
         &     cos(i)*cos(om)*cos(ompl))                                    
    A8 = -2.d0*a*apl*beta*(-cos(ompl)*sin(om)+                        &
         &     cos(i)*cos(om)*sin(ompl))                                    
    A9 = -2.d0*a*apl*betapl*(-cos(om)*sin(ompl)+                      &
         &     cos(i)*sin(om)*cos(ompl))                                    
    A10 = -2.d0*a*apl*(cos(om)*cos(ompl)+cos(i)*sin(om)*sin(ompl)) 
    A11 = -2.d0*a*apl*beta*epl*(sin(om)*cos(ompl)-                    &
         &     sin(ompl)*cos(om)*cos(i))                                    
    A12 = -2.d0*a*(a*e - apl*epl*(cos(om)*cos(ompl)+                  &
         &     sin(om)*sin(ompl)*cos(i)))                                   
    A13 = -2.d0*a*apl*betapl*e*(sin(ompl)*cos(om)-                    &
         &     cos(ompl)*sin(om)*cos(i))                                    
    A14 = 2.d0*apl*(a*e*cos(om)*cos(ompl)-apl*epl +                   &
         &     a*e*sin(ompl)*sin(om)*cos(i))                                
    A15 = ((Ca-Cp)**2)+(Fa**2)+(Fp**2)-2.d0*cos(i)*Fa*Fp 
    
! =======================================================           
! CHECK of LEADING COEFFICIENT                                      
    tmpcoe = 32.d0*A9*A8**3*A7*A14 - 32.d0*A10*A11**2*A14*A7**2 -           &
         &128.d0*A4*A8**2*A6*A7**2 + 64.d0*A4*A11**3*A14*A8 -               &
         & 64.d0*A4*A8**3*A14*A11 + 64.d0*A4*A8**3*A10*A11 + 32.d0          &
         &*A7**3*A9*A14*A8 - 32.d0*A10*A11**2*A14*A8**2 - 32.d0*A7**3*A13*  &
         & A14*A8 + 128.d0*A4*A8**2*A6*A11**2 + 64.d0*A6*A8**3*A14*A11 -    &
         & 64*A6*A11**3*A14*A8 + 64*A7**3*A9*A6*A11 - 64*A7**3*A13*A6*A11 - &
         & 32.d0*A9*A11**2*A13*A7**2 + 64.d0*A10*A11**3*A6*A8 +             &
         & 64.d0*A13*A11**3*A7*A6 + 64.d0*A9*A11**3*A7*A4 - 32.d0*A10*A8*   &
         & A7**3*A9 - 64.d0*A7**3*A9*A4*A11 + 64.d0*A7**3*A13*A4*A11 - 64.d0&
         &*A10*A8**3*A6*A11 + 32.d0*A13*A8**3*A7*A10 + 32.d0*A9*A8**2*A13*  &
         &A7**2 - 64.d0*A9*A11**3*A7*A6 - 32.d0*A13*A8**3*A7*A14 + 32.d0*A10&
         &*A8**2*A14*A7**2 - 32.d0*A9*A11**2*A13*A8**2 - 32.d0*A9*A8**3*    &
         &A7*A10 - 64.d0*A4*A11**3*A10*A8                                   
    
    coe20 = tmpcoe + 128.d0*A4*A11**2*A6*A7**2 + 32.d0*A10*A8*A7**3*        &
         &A13 - 16.d0*A7**4*A13**2 + 64.d0*A4**2*A11**4 - 16.d0*A10**2*     &
         &A8**4 - 16.d0*A9**2*A8**2*A7**2 - 64.d0*A4**2*A11**2*A7**2 +      &
         &16.d0*A10**2*A11**2*A7**2 - 64.d0*A6**2*A11**2*A7**2 + 16.d0*     &
         &A13**2*A11**2*A7**2 - 16.d0*A10**2*A8**2*A7**2 + 16.d0*A9**2*     &
         &A11**2*A8**2 - 64.d0*A6**2*A8**2*A11**2 + 16.d0*A9**2*A11**2*A7**2&
         &- 16.d0*A14**2*A8**2*A7**2 + 16.d0*A13**2*A11**2*A8**2 - 16.d0*   &
         &A13**2*A8**2*A7**2 - 16.d0*A14**2*A8**4 + 64.d0*A6**2*A11**4      &
         &- 16.d0*A7**4*A9**2 - 64.d0*A13*A11**3*A7*A4 + 16.d0*A14**2*      &
         &A11**2*A7**2 - 128.d0*A4*A11**4*A6 + 64.d0*A4**2*A8**2*A7**2 +    &
         &32.d0*A7**4*A13*A9 - 64.d0*A4**2*A11**2*A8**2 + 16.d0*A14**2*     &
         &A8**2*A11**2 + 32.d0*A10*A8**4*A14 + 16.d0*A10**2*A8**2*A11**2    &
         &+ 64.d0*A6**2*A8**2*A7**2                                         
    
    IF(abs(coe20).lt.1.d-15) THEN 
       WRITE(*,*)'aical: 20th COEFFICIENT is small',coe20 
    ENDIF
    
    RETURN 
  END SUBROUTINE aical

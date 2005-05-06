
! ******************************************************************
! ********** SELECT AMONG MINIMA, MAXIMA and SADLE POINTS **********
! ******************************************************************
! *********** written by GIOVANNI F. GRONCHI (1999) ****************
! ********** Department of Mathematics, UNIVERSITY of PISA *********
! last modified June 2004
! ==================================================================
  SUBROUTINE hess(u,upl,ans) 
    USE fund_const
    IMPLICIT NONE 
!   u,upl are passed in radians                                       
    DOUBLE PRECISION,INTENT(IN) :: u,upl 
    INTEGER,INTENT(OUT) :: ans !  ans = 1: maximum 
!                                 ans = 0: saddle 
!                                 ans= -1: minimum 
!                                 ans= -2: cannot decide
!   -------- end interface -----------------------------------------
    DOUBLE PRECISION :: H11,H12,H21,H22,trH,detH 
!   eigenvalues                                                       
    DOUBLE PRECISION :: lam1,lam2 
    DOUBLE PRECISION, PARAMETER :: eps=1.d-20 !tolerance parameter
!                                                                       
    DOUBLE PRECISION :: A1,A2,A3,A4,A5,A6,A7,A8,A9,A10 
    DOUBLE PRECISION :: A11,A12,A13,A14,A15 
    COMMON/Aj1to15/ A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15 
!   ==================================================================
                                                                        
!   HESSIAN MATRIX                                                    
    H11 = 2.d0*A1*cos(u)**2-2.d0*A1*sin(u)**2+2.d0*A3*sin(u)**2-      &
         & 2.d0*A3*cos(u)**2-A7*sin(u)*sin(upl)-A8*sin(u)*cos(upl)-     &
         & A9*cos(u)*sin(upl)-A10*cos(u)*cos(upl)-A11*sin(u)-A12*cos(u) 
    H12 = A7*cos(u)*cos(upl)-A8*cos(u)*sin(upl)-A9*sin(u)*cos(upl)+   &
         & A10*sin(u)*sin(upl)                                          
    H21 = H12 
    H22 = 2.d0*A4*cos(upl)**2-2.d0*A4*sin(upl)**2+2.d0*A6*sin(upl)**2-&
         & 2.d0*A6*cos(upl)**2-A7*sin(u)*sin(upl)-A8*sin(u)*cos(upl)-   &
         & A9*cos(u)*sin(upl)-A10*cos(u)*cos(upl)-A13*sin(upl)-         &
         & A14*cos(upl)                                                 
!   TRACE OF THE HESSIAN                                              
    trH = H11 + H22 
!   DETERMINANT OF THE HESSIAN                                        
    detH = H11*H22 - H12*H21 
!   check                                                             
    IF (abs(detH).le.0.d0) THEN 
!       WRITE(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!' 
       WRITE(*,*)'hess: Hessian determinant <= 0! det(H) =',detH 
!       WRITE(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!' 
       ans = -2
       GOTO 10 
    ENDIF
!   compute eigenvalues                                               
    IF ((trH**2-4.d0*detH).ge.0.d0) THEN 
       lam1 = (trH + dsqrt(trH**2 - 4.d0*detH))/2.d0 
       lam2 = (trH - dsqrt(trH**2 - 4.d0*detH))/2.d0 
       IF ((lam1.gt.eps).and.(lam2.gt.eps)) THEN 
          ans = -1 
       ELSEIF ((lam1.lt.-eps).and.(lam2.lt.-eps)) THEN 
          ans = 1 
       ELSEIF ((lam1.gt.eps).and.(lam2.lt.-eps)) THEN 
          ans = 0 
       ELSEIF ((lam1.lt.-eps).and.(lam2.gt.eps)) THEN 
          ans = 0 
       ELSE 
          ans = -2
       ENDIF
    ELSE 
       WRITE(*,*)'hess: negative discriminant! complex eigenvalues!' 
       ans = -2
    ENDIF
    
10  CONTINUE 
    
    RETURN 
  END SUBROUTINE hess


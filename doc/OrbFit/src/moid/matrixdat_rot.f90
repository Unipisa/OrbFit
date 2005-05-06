
!     ******************************************************************
!     ****** Given Ai, this routine gets alpha,beta,gamma,A,B,C,D,E ****
!     ******************************************************************
!     *********** written by GIOVANNI F. GRONCHI (2003) ****************
!     ********** Department of Mathematics, UNIVERSITY of PISA *********
!     ==================================================================
      SUBROUTINE matrixdat_rot(util,upltil,N,alp,bet,gam,A,B,C,D,E) 
      IMPLICIT NONE 
      DOUBLE PRECISION,INTENT(IN) :: util,upltil ! rotation angles 
!     number of evaluations                                             
      INTEGER,INTENT(IN) :: N 
!     Sylvester matrix elements                                         
      DOUBLE PRECISION,INTENT(OUT) :: alp(N),bet(N),gam(N) 
      DOUBLE PRECISION,INTENT(OUT) :: A(N),B(N),C(N),D(N),E(N) 
!     ========== end interface =========================================
      DOUBLE PRECISION :: sutil,cutil,supltil,cupltil !auxiliary
      DOUBLE PRECISION :: A1,A2,A3,A4,A5,A6,A7,A8,A9,A10 
      DOUBLE PRECISION :: A11,A12,A13,A14,A15 
      COMMON/Aj1to15/ A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15 
!     loop indexes                                                      
      INTEGER :: i,j,l 
!     ==================================================================
 
      sutil=dsin(util)
      cutil=dcos(util)
      supltil=dsin(upltil)
      cupltil=dcos(upltil)
                                                                       
      alp(1) = 2.d0*A1*cutil*sutil - 2.d0*A3*cutil*sutil + &
           & A11*cutil - A12*sutil - A8*cutil*cupltil - &
           & A7*cutil*supltil + A10*sutil*cupltil + A9*sutil*supltil

      alp(2) = 2.d0*A8*sutil*cupltil + 2.d0*A9*cutil*supltil + &
           & 2.d0*A7*sutil*supltil + 2.d0*A10*cutil*cupltil - &
           & 2.d0*A11*sutil - 2.d0*A12*cutil - 8.d0*A3*cutil**2 +&
           & 8.d0*A1*cutil**2 + 4.d0*A3 - 4.d0*A1

      alp(3) = -12.d0*cutil*sutil*(-A3+A1)

      alp(4) = -2.d0*A11*sutil - 2.d0*A12*cutil - 4.d0*A3 + 4.d0*A1 + &
           & 2.d0*A7*sutil*supltil + 2.d0*A10*cutil*cupltil + &
           & 8.d0*A3*cutil**2 - 8.d0*A1*cutil**2 + &
           & 2.d0*A8*sutil*cupltil + 2.d0*A9*cutil*supltil

      alp(5) = -A11*cutil - 2.d0*A3*cutil*sutil + 2.d0*A1*cutil*sutil - &
           & A10*sutil*cupltil - A9*sutil*supltil + A12*sutil + &
           & A8*cutil*cupltil + A7*cutil*supltil

!      alp(1) = -A8+A11 
!      alp(2) = 4.d0*A1-4.d0*A3+2.d0*A10-2.d0*A12 
!      alp(3) = 0.d0 
!      alp(4) = 2.d0*A10+4.d0*A3-4.d0*A1-2.d0*A12 
!      alp(5) = -A11+A8 

      bet(1) = 2.d0*A7*cutil*cupltil + 2.d0*A10*sutil*supltil - &
           & 2.d0*A8*cutil*supltil - 2.d0*A9*sutil*cupltil

      bet(2) = -4.d0*A7*sutil*cupltil + 4.d0*A8*sutil*supltil - &
           & 4.d0*A9*cutil*cupltil + 4.d0*A10*cutil*supltil

      bet(3) = 0.d0

      bet(4) = -4.d0*A7*sutil*cupltil + 4.d0*A8*sutil*supltil - &
           & 4.d0*A9*cutil*cupltil + 4.d0*A10*cutil*supltil

      bet(5) = 2.d0*A9*sutil*cupltil - 2.d0*A10*sutil*supltil + &
           & 2.d0*A8*cutil*supltil - 2.d0*A7*cutil*cupltil
                                                                        
!      bet(1) = 2.d0*A7 
!      bet(2) = -4.d0*A9 
!      bet(3) = 0.d0 
!      bet(4) = -4.d0*A9 
!      bet(5) = -2.d0*A7 
 
      gam(1) = 2.d0*A1*cutil*sutil - 2.d0*A3*cutil*sutil + &
           & A11*cutil - A12*sutil + A8*cutil*cupltil + &
           & A7*cutil*supltil - A10*sutil*cupltil - A9*sutil*supltil

      gam(2) = -2.d0*A8*sutil*cupltil - 2.d0*A9*cutil*supltil -&
           & 2.d0*A7*sutil*supltil - 2.d0*A10*cutil*cupltil - &
           & 2.d0*A11*sutil - 2.d0*A12*cutil - 8.d0*A3*cutil**2 + &
           & 8.d0*A1*cutil**2 + 4.d0*A3 - 4.d0*A1

      gam(3) = -12.d0*cutil*sutil*(-A3+A1)

      gam(4) = -2.d0*A11*sutil - 2.d0*A12*cutil - 4.d0*A3 + &
           & 4.d0*A1 - 2.d0*A7*sutil*supltil - 2.d0*A10*cutil*cupltil +&
           8.d0*A3*cutil**2 - 8.d0*A1*cutil**2 - 2.d0*A8*sutil*cupltil -&
           & 2.d0*A9*cutil*supltil

      gam(5) = -A11*cutil - 2.d0*A3*cutil*sutil + 2.d0*A1*cutil*sutil +&
           & A10*sutil*cupltil + A9*sutil*supltil + A12*sutil - &
           & A8*cutil*cupltil - A7*cutil*supltil

!      gam(1) = A8+A11 
!      gam(2) = 4.d0*A1-2.d0*A12-2.d0*A10-4.d0*A3 
!      gam(3) = 0.d0 
!      gam(4) = -2.d0*A12 - 4.d0*A1 + 4.d0*A3 - 2.d0*A10 
!      gam(5) = -A11-A8 
                                                                        
      A(1) = -A13*cupltil + A14*supltil + A10*cutil*supltil - &
           & A9*cutil*cupltil + A8*sutil*supltil - A7*sutil*cupltil + &
           & 2.d0*A4*cupltil*supltil - 2.d0*A6*cupltil*supltil

      A(2) = 2.d0*A9*sutil*cupltil - 2.d0*A10*sutil*supltil + &
           & 2.d0*A8*cutil*supltil - 2.d0*A7*cutil*cupltil

      A(3) = 2.d0*A4*cupltil*supltil - 2.d0*A6*cupltil*supltil - &
           & A13*cupltil + A7*sutil*cupltil - A8*sutil*supltil + &
           & A9*cutil*cupltil - A10*cutil*supltil + A14*supltil

!      A(1) = -A9-A13 
!      A(2) = -2.d0*A7 
!      A(3) = A9-A13 
                                                                        
      B(1) = 8.d0*A6*cupltil**2 - 8.d0*A4*cupltil**2 -2.d0*A14*cupltil -&
           & 2.d0*A13*supltil - 4.d0*A6 + 4.d0*A4 - 2.d0*A7*sutil*supltil - &
           & 2.d0*A10*cutil*cupltil - 2.d0*A9*cutil*supltil - &
           & 2.d0*A8*sutil*cupltil

      B(2) = -4.d0*A7*cutil*supltil + 4.d0*A9*sutil*supltil + &
           & 4.d0*A10*sutil*cupltil - 4.d0*A8*cutil*cupltil

      B(3) = -4.d0*A6 + 4.d0*A4 + 2.d0*A7*sutil*supltil + &
           & 2.d0*A9*cutil*supltil + 2.d0*A8*sutil*cupltil + &
           & 2.d0*A10*cutil*cupltil + 8.d0*A6*cupltil**2 - &
           & 8.d0*A4*cupltil**2 - 2.d0*A14*cupltil - 2.d0*A13*supltil

!      B(1) = -4.d0*A4+4.d0*A6-2.d0*A10-2.d0*A14 
!      B(2) = -4.d0*A8 
!      B(3) = 2.d0*A10-4.d0*A4+4.d0*A6-2.d0*A14 
                                                                        
      C(1) = 12.d0*cupltil*supltil*(-A4+A6)

      C(2) = 0.d0

      C(3) = 12.d0*cupltil*supltil*(-A4+A6)

!      C(1) = 0.d0 
!      C(2) = 0.d0 
!      C(3) = 0.d0 
                                                                        
      D(1) = -2.d0*A13*supltil - 8.d0*A6*cupltil**2 + 8.d0*A4*cupltil**2 - &
           & 2.d0*A14*cupltil + 4.d0*A6 - 4.d0*A4 - 2.d0*A7*sutil*supltil -&
           & 2.d0*A9*cutil*supltil - 2.d0*A8*sutil*cupltil - &
           & 2.d0*A10*cutil*cupltil

      D(2) = -4.d0*A7*cutil*supltil + 4.d0*A9*sutil*supltil + &
           & 4.d0*A10*sutil*cupltil - 4.d0*A8*cutil*cupltil

      D(3) = 4.d0*A6 - 4.d0*A4 + 2.d0*A7*sutil*supltil + &
           & 2.d0*A9*cutil*supltil + 2.d0*A8*sutil*cupltil + &
           & 2.d0*A10*cutil*cupltil - 8.d0*A6*cupltil**2 + &
           & 8.d0*A4*cupltil**2 - 2.d0*A14*cupltil - 2.d0*A13*supltil

!      D(1) = -2.d0*A14+4.d0*A4-4.d0*A6-2.d0*A10 
!      D(2) = -4.d0*A8 
!      D(3) = 2.d0*A10-2.d0*A14+4.d0*A4-4.d0*A6 
                                                                        
      E(1) = A7*sutil*cupltil - A8*sutil*supltil - A10*cutil*supltil - &
           & 2.d0*A6*cupltil*supltil + 2.d0*A4*cupltil*supltil + &
           & A13*cupltil - A14*supltil + A9*cutil*cupltil

      E(2) = 2.d0*A7*cutil*cupltil + 2.d0*A10*sutil*supltil - &
           & 2.d0*A8*cutil*supltil - 2.d0*A9*sutil*cupltil

      E(3) = -A14*supltil + A13*cupltil - A7*sutil*cupltil + &
           & A10*cutil*supltil - A9*cutil*cupltil + A8*sutil*supltil + &
           & 2.d0*A4*cupltil*supltil - 2.d0*A6*cupltil*supltil

!      E(1) = A9+A13 
!      E(2) = 2.d0*A7 
!      E(3) = -A9+A13 
                                                                        
!     settings of coefficients of higher degree to zero                 
      DO j = 6,N 
         alp(j) = 0.d0 
         bet(j) = 0.d0 
         gam(j) = 0.d0 
      ENDDO 
                                                                        
      DO l = 4,N 
         A(l) = 0.d0 
         B(l) = 0.d0 
         C(l) = 0.d0 
         D(l) = 0.d0 
         E(l) = 0.d0 
      ENDDO 
                                                                        
      RETURN 
    END SUBROUTINE matrixdat_rot

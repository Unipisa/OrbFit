c     ********************************************************************
c     ****** Given Ai, this routine gets alpha,beta,gamma,A,B,C,D,E ******
c     ********************************************************************
c     *********** written by GIOVANNI F. GRONCHI (2001) ******************
c     ********** Department of Mathematics, UNIVERSITY of PISA ***********
c     ====================================================================
      SUBROUTINE matrixdat(N,alp,bet,gam,A,B,C,D,E)
      IMPLICIT NONE
c     number of evaluations
      INTEGER N
c     Sylvester matrix elements
      DOUBLE PRECISION alp( * ),bet( * ),gam( * )
      DOUBLE PRECISION A( * ),B( * ),C( * ),D( * ),E( * )
c     --------------------------------------- end interface --------------
      DOUBLE PRECISION A1,A2,A3,A4,A5,A6,A7,A8,A9,A10
      DOUBLE PRECISION A11,A12,A13,A14,A15
      COMMON/Aj1to15/ A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15
c     loop indexes
      INTEGER i,j,l
c     ====================================================================

      alp(1) = -A8+A11
      alp(2) = 4.d0*A1-4.d0*A3+2.d0*A10-2.d0*A12
      alp(3) = 0.d0
      alp(4) = 2.d0*A10+4.d0*A3-4.d0*A1-2.d0*A12
      alp(5) = -A11+A8

      bet(1) = 2.d0*A7
      bet(2) = -4.d0*A9
      bet(3) = 0.d0
      bet(4) = -4.d0*A9
      bet(5) = -2.d0*A7

      gam(1) = A8+A11 
      gam(2) = 4.d0*A1-2.d0*A12-2.d0*A10-4.d0*A3
      gam(3) = 0.d0
      gam(4) = -2.d0*A12-4.d0*A1+4.d0*A3-2.d0*A10
      gam(5) = -A11-A8

      A(1) = -A9-A13
      A(2) = -2.d0*A7
      A(3) = A9-A13 

      B(1) = -4.d0*A4+4.d0*A6-2.d0*A10-2.d0*A14
      B(2) = -4.d0*A8
      B(3) = 2.d0*A10-4.d0*A4+4.d0*A6-2.d0*A14

      C(1) = 0.d0 
      C(2) = 0.d0
      C(3) = 0.d0

      D(1) = -2.d0*A14+4.d0*A4-4.d0*A6-2.d0*A10
      D(2) = -4.d0*A8
      D(3) = 2.d0*A10-2.d0*A14+4.d0*A4-4.d0*A6

      E(1) = A9+A13
      E(2) = 2.d0*A7
      E(3) = -A9+A13
 
c     settings of coefficients of higher degree to zero 
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
      END

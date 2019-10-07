c     ********************************************************************
c     ************* COMPUTE COEFFICIENTS OF SYLVESTER MATRIX *************
c     ********************************************************************
c     *********** written by GIOVANNI F. GRONCHI (2001) ******************
c     ********** Department of Mathematics, UNIVERSITY of PISA ***********
c     ====================================================================
      SUBROUTINE clevcompsylv(alpha,beta,gamma,A,B,C,D,E,SSYLV)
      IMPLICIT NONE
c     Sylvester matrix elements 
      COMPLEX*16 alpha,beta,gamma
      COMPLEX*16 A,B,C,D,E
      COMPLEX*16 SSYLV(6,6)
c     ---------------------------------- end interface -------------------
c     loop index
      INTEGER h,k
c     ====================================================================

      SSYLV(1,1) = alpha 
      SSYLV(1,2) = (0.d0,0.d0)
      SSYLV(1,3) = (0.d0,0.d0)
      SSYLV(1,4) = (0.d0,0.d0)
      SSYLV(1,5) = A
      SSYLV(1,6) = (0.d0,0.d0)
c
      SSYLV(2,1) = beta
      SSYLV(2,2) = alpha
      SSYLV(2,3) = (0.d0,0.d0)
      SSYLV(2,4) = (0.d0,0.d0)
      SSYLV(2,5) = B
      SSYLV(2,6) = A
c
      SSYLV(3,1) = gamma
      SSYLV(3,2) = beta
      SSYLV(3,3) = alpha
      SSYLV(3,4) = (0.d0,0.d0)
      SSYLV(3,5) = C
      SSYLV(3,6) = B
c
      SSYLV(4,1) = (0.d0,0.d0)
      SSYLV(4,2) = gamma
      SSYLV(4,3) = beta
      SSYLV(4,4) = alpha
      SSYLV(4,5) = D
      SSYLV(4,6) = C
c
      SSYLV(5,1) = (0.d0,0.d0)
      SSYLV(5,2) = (0.d0,0.d0)
      SSYLV(5,3) = gamma
      SSYLV(5,4) = beta
      SSYLV(5,5) = E
      SSYLV(5,6) = D
c
      SSYLV(6,1) = (0.d0,0.d0)
      SSYLV(6,2) = (0.d0,0.d0)
      SSYLV(6,3) = (0.d0,0.d0)
      SSYLV(6,4) = gamma
      SSYLV(6,5) = (0.d0,0.d0)
      SSYLV(6,6) = E

      RETURN
      END

c     ********************************************************************
c     ******** computes the determinant of matrix A(6,6) employing *******
c     ******** LU factorization (A(6,6) has complex coefficients) ********
c     ********************************************************************
c     *********** written by GIOVANNI F. GRONCHI (2001) ******************
c     ********** Department of Mathematics, UNIVERSITY of PISA ***********
c     ====================================================================
      SUBROUTINE cdetcomp(A,detA)
      IMPLICIT NONE
c     A is a matrix 6x6
      COMPLEX*16  A(6,6),detA
c     -------------------------------------- end interface ---------------
      INTEGER            INFO, LDA, M, N
c     IPIV must have the same dimension of N 
      INTEGER            IPIV(6)
      INTEGER signdetA
c     loop indexes
      INTEGER i,j
c     ====================================================================

      LDA= 6
c     number of columns of A
      N=6
c     number of rows of A 
      M=6
c     initialization
      detA=1.d0
      signdetA=1
c     compute LU factorization
      CALL ZGETF2( M, N, A, LDA, IPIV, INFO )

      DO i=1,N
         IF (IPIV(i).ne.i) THEN
            signdetA = -signdetA
         ENDIF   
      ENDDO
c     compute determinant of A
      DO i=1,N
         detA=detA*A(i,i)
      ENDDO      
      detA = signdetA*detA
      
      RETURN
      END

!     ******************************************************************
!     ******** computes the determinant of matrix A(6,6) employing *****
!     ******** LU factorization (A(6,6) has complex coefficients) ******
!     ******************************************************************
!     *********** written by GIOVANNI F. GRONCHI (2001) ****************
!     ********** Department of Mathematics, UNIVERSITY of PISA *********
!     ==================================================================
      SUBROUTINE cdetcomp(A,detA) 
      IMPLICIT NONE 
!     A is a matrix 6x6                                                 
      COMPLEX*16  A(6,6),detA 
!     -------------------------------------- end interface -------------
      INTEGER            INFO, LDA, M, N 
!     IPIV must have the same dimension of N                            
      INTEGER            IPIV(6) 
      INTEGER signdetA 
!     loop indexes                                                      
      INTEGER i,j 
!     ==================================================================
                                                                        
      LDA= 6 
!     number of columns of A                                            
      N=6 
!     number of rows of A                                               
      M=6 
!     initialization                                                    
      detA=1.d0 
      signdetA=1 
!     compute LU factorization                                          
      CALL ZGETF2( M, N, A, LDA, IPIV, INFO ) 
                                                                        
      DO i=1,N 
         IF (IPIV(i).ne.i) THEN 
            signdetA = -signdetA 
         ENDIF 
      ENDDO 
!     compute determinant of A                                          
      DO i=1,N 
         detA=detA*A(i,i) 
      ENDDO 
      detA = signdetA*detA 
                                                                        
      RETURN 
      END                                           

! *********************************************************
! ***** COMPUTE COEFFICIENTS OF SYLVESTER MATRIX **********
! *********************************************************
  SUBROUTINE compsylv48_QP(AA,BB,SYLV) 
    IMPLICIT NONE 
  INTEGER, PARAMETER :: qkind=kind(1.q0) !for quadruple precision
! Sylvester matrix elements                                         
    COMPLEX(KIND=qkind),INTENT(IN) :: AA(0:20)
    COMPLEX(KIND=qkind),INTENT(IN) :: BB(0:2)
    COMPLEX(KIND=qkind),INTENT(OUT) :: SYLV(1:22,1:22) 
! ======== end interface ==================================
    INTEGER :: i,j
    SYLV=QCMPLX(0.q0,0.q0)
    DO i = 1,21
       SYLV(i,1) = AA(21-i)
       SYLV(i+1,2) = AA(21-i)
    ENDDO
    DO j = 3,22
       SYLV(j-2,j) = BB(2)
       SYLV(j-1,j) = BB(1)
       SYLV(j,j) = BB(0)
    ENDDO
  END SUBROUTINE compsylv48_QP

! *********************************************
! DOUBLE PRECISION VERSION
! ********************************************* 
  SUBROUTINE compsylv48_DP(AA,BB,SYLV) 
    IMPLICIT NONE 
  INTEGER, PARAMETER :: dkind=kind(1.d0) !for quadruple precision
! Sylvester matrix elements                                         
    COMPLEX(KIND=dkind),INTENT(IN) :: AA(0:20)
    COMPLEX(KIND=dkind),INTENT(IN) :: BB(0:2)
    COMPLEX(KIND=dkind),INTENT(OUT) :: SYLV(1:22,1:22) 
! ======== end interface ==================================
    INTEGER :: i,j
    SYLV=QCMPLX(0.q0,0.q0)
    DO i = 1,21
       SYLV(i,1) = AA(21-i)
       SYLV(i+1,2) = AA(21-i)
    ENDDO
    DO j = 3,22
       SYLV(j-2,j) = BB(2)
       SYLV(j-1,j) = BB(1)
       SYLV(j,j) = BB(0)
    ENDDO
  END SUBROUTINE compsylv48_DP

! ******************************************************************
! ******** computes the determinant of matrix A(22,22) employing *****
! ******** LU factorization (A(22,22) has complex coefficients) ******
! ==================================================================
  SUBROUTINE cdetcomp22_QP(A,detA) 
    IMPLICIT NONE 
    INTEGER, PARAMETER :: qkind=kind(1.q0) !for quadruple precision
! A is a matrix 22x22
    COMPLEX(KIND=qkind),INTENT(IN) :: A(22,22)
    COMPLEX(KIND=qkind),INTENT(OUT) :: detA 
! ============== end interface =================
    INTEGER :: INFO, LDA, M, N 
! IPIV must have the same dimension of N
    INTEGER :: IPIV(22)
    INTEGER :: signdetA 
! loop indexes                                                      
    INTEGER :: i,j 
! ==================================================================
    LDA= 22 
! number of columns of A                                            
    N=22
! number of rows of A                                               
    M=22
! initialization                                                    
    detA=1.q0 
    signdetA=1 
! compute LU factorization                                          
    CALL ZGETF2_QP( M, N, A, LDA, IPIV, INFO ) 
    DO i=1,N 
       IF (IPIV(i).ne.i) THEN 
          signdetA = -signdetA 
       ENDIF
    ENDDO
! compute determinant of A                                          
    DO i=1,N 
       detA=detA*A(i,i) 
    ENDDO   
    detA = signdetA*detA 
  END SUBROUTINE cdetcomp22_QP

! *********************************************
! DOUBLE PRECISION VERSION
! ********************************************* 
  SUBROUTINE cdetcomp22_DP(A,detA) 
    IMPLICIT NONE 
    INTEGER, PARAMETER :: dkind=kind(1.d0) !for quadruple precision
! A is a matrix 22x22
    COMPLEX(KIND=dkind),INTENT(IN) :: A(22,22)
    COMPLEX(KIND=dkind),INTENT(OUT) :: detA 
! ============== end interface =================
    INTEGER :: INFO, LDA, M, N 
! IPIV must have the same dimension of N
    INTEGER :: IPIV(22)
    INTEGER :: signdetA 
! loop indexes                                                      
    INTEGER :: i,j 
! ==================================================================
    LDA= 22 
! number of columns of A                                            
    N=22
! number of rows of A                                               
    M=22
! initialization                                                    
    detA=1.q0 
    signdetA=1 
! compute LU factorization                                          
    CALL ZGETF2( M, N, A, LDA, IPIV, INFO ) 
    DO i=1,N 
       IF (IPIV(i).ne.i) THEN 
          signdetA = -signdetA 
       ENDIF
    ENDDO
! compute determinant of A                                          
    DO i=1,N 
       detA=detA*A(i,i) 
    ENDDO   
    detA = signdetA*detA 
  END SUBROUTINE cdetcomp22_DP


! ==================================================================
! GIVEN rho2 COMPUTE rho1 SUCH THAT (rho1,rho2) 
! SATISFIES THE SYSTEM p = q = 0
! ==================================================================
  SUBROUTINE solvesys_QP(pmod,qmod,ncoe_p,ncoe_q,nroots,roots, &
       & rho1,rho2,nposroots)
    USE output_control, ONLY: verb_2body
    IMPLICIT NONE
    INTEGER, PARAMETER :: qkind=kind(1.q0) !for quadruple precision
    INTEGER,INTENT(IN) :: pmod,qmod
    REAL(KIND=qkind),INTENT(IN) :: ncoe_p(0:24,0:24),ncoe_q(0:2,0:2)
    INTEGER,INTENT(IN) :: nroots 
    REAL(KIND=qkind),INTENT(IN) :: roots(nroots)
    REAL(KIND=qkind),INTENT(OUT) :: rho1(nroots),rho2(nroots)
    INTEGER,INTENT(OUT) :: nposroots 
! ======== end interface ===============================    
    REAL(KIND=qkind) :: rho1sol(2,nroots),rho2sol(nroots),peval(2),qeval(2) ! auxiliary 
    REAL(KIND=qkind) :: bb(0:2)
    REAL(KIND=qkind) :: discrim
    INTEGER :: h,j

    nposroots=0

    DO h=1,nroots
       rho2sol(h)=roots(h)
!       rho2(1)= 3.0340706822328d0
!       DO j=0,2
!          bb(j) = ncoe_q(j,2)*rho2(h)**2 + ncoe_q(j,1)*rho2(h) + ncoe_q(j,0)
!       ENDDO
       bb(0)= ncoe_q(0,2)*rho2sol(h)**2 + ncoe_q(0,1)*rho2sol(h) + ncoe_q(0,0)
       bb(1)= ncoe_q(1,0)
       bb(2)= ncoe_q(2,0)

! 2nd degree polynomial equation
!    bb(2)*rho1**2 + bb(1)*rho1 + bb(0) = 0
       discrim=bb(1)**2-4.d0*bb(0)*bb(2)
       IF(discrim.gt.0.q0)THEN
          rho1sol(1,h) = (-bb(1) + sqrt(bb(1)**2-4.d0*bb(0)*bb(2)))/ &
               & (2.d0*bb(2))
          rho1sol(2,h) = (-bb(1) - sqrt(bb(1)**2-4.d0*bb(0)*bb(2)))/ &
               & (2.d0*bb(2))
!       write(*,*)'rho1, rho2',rho1sol(1,h),rho2(h)
!       write(*,*)'rho1, rho2',rho1sol(2,h),rho2(h)
          CALL poly_eval_QP(pmod,ncoe_p,rho1sol(1,h),rho2sol(h),peval(1))  
          CALL poly_eval_QP(pmod,ncoe_p,rho1sol(2,h),rho2sol(h),peval(2))  
!          write(*,*)'p(rho1_1,rho2)',peval(1)
!          write(*,*)'p(rho1_2,rho2)',peval(2)
          CALL poly_eval_QP(qmod,ncoe_q,rho1sol(1,h),rho2sol(h),qeval(1))  
          CALL poly_eval_QP(qmod,ncoe_q,rho1sol(2,h),rho2sol(h),qeval(2))  

113 FORMAT(A19,2X,3(F14.9,2X))
!          select only positive rho1
          IF(rho1sol(1,h).gt.0.q0.OR.rho1sol(2,h).gt.0.q0)THEN
             nposroots=nposroots+1
             rho2(nposroots) = rho2sol(h)
!             write(*,113)'almeno una positiva',rho2sol(h),rho1sol(1,h),rho1sol(2,h)
             
             IF(rho1sol(1,h).gt.0.q0)THEN
                rho1(nposroots) = rho1sol(1,h)
                IF(abs(peval(2)).lt.abs(peval(1)).AND. &
                     & rho1sol(2,h).gt.0.q0)THEN
                   rho1(nposroots)= rho1sol(2,h)
                ENDIF
             ELSE 
 !            ELSEIF(rho1sol(2,h).gt.0.q0)THEN
                rho1(nposroots) = rho1sol(2,h)
 !               IF(abs(peval(1)).lt.abs(peval(2)).AND. &
 !                    & rho1sol(1,h).gt.0.q0)THEN
 !                  rho1(nposroots)= rho1sol(1,h)
 !               ENDIF
             ENDIF
             
          ELSE
             IF(verb_2body.gt.10)THEN
                WRITE(*,*)'both negative solutions! rho1,rho2', &
                  & rho1sol(1,h),rho1sol(2,h)
             ENDIF
          ENDIF
          
       ELSE
          IF(verb_2body.gt.10)THEN
             write(*,*)'negative discriminant',discrim
          ENDIF
!          rho1sol(1,h) = 0.d0 
!          rho1sol(2,h) = 0.d0 
       ENDIF
    ENDDO

  END SUBROUTINE solvesys_QP


! ==================================================================
! GIVEN rho2 COMPUTE rho1 SUCH THAT (rho1,rho2) 
! SATISFIES THE SYSTEM p = q = 0  (DOUBLE PRECISION)
! ==================================================================
  SUBROUTINE solvesys(pmod,qmod,ncoe_p,ncoe_q,nroots,roots, &
       & rho1,rho2,nposroots)
    USE fund_const, ONLY: dkind
    USE output_control, ONLY: verb_2body
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: pmod,qmod
    REAL(KIND=dkind),INTENT(IN) :: ncoe_p(0:24,0:24),ncoe_q(0:2,0:2)
    INTEGER,INTENT(IN) :: nroots 
    REAL(KIND=dkind),INTENT(IN) :: roots(nroots)
    REAL(KIND=dkind),INTENT(OUT) :: rho1(nroots),rho2(nroots)
    INTEGER,INTENT(OUT) :: nposroots 
! ======== end interface ===============================    
    REAL(KIND=dkind) :: rho1sol(2,nroots),rho2sol(nroots),peval(2),qeval(2) ! auxiliary 
    REAL(KIND=dkind) :: bb(0:2)
    REAL(KIND=dkind) :: discrim
    INTEGER :: h,j

    nposroots=0

    DO h=1,nroots
       rho2sol(h)=roots(h)
!       rho2(1)= 3.0340706822328d0
!       DO j=0,2
!          bb(j) = ncoe_q(j,2)*rho2(h)**2 + ncoe_q(j,1)*rho2(h) + ncoe_q(j,0)
!       ENDDO
       bb(0)= ncoe_q(0,2)*rho2sol(h)**2 + ncoe_q(0,1)*rho2sol(h) + ncoe_q(0,0)
       bb(1)= ncoe_q(1,0)
       bb(2)= ncoe_q(2,0)

! 2nd degree polynomial equation
!    bb(2)*rho1**2 + bb(1)*rho1 + bb(0) = 0
       discrim=bb(1)**2-4.d0*bb(0)*bb(2)
       IF(discrim.gt.0.d0)THEN
          rho1sol(1,h) = (-bb(1) + sqrt(bb(1)**2-4.d0*bb(0)*bb(2)))/ &
               & (2.d0*bb(2))
          rho1sol(2,h) = (-bb(1) - sqrt(bb(1)**2-4.d0*bb(0)*bb(2)))/ &
               & (2.d0*bb(2))
!       write(*,*)'rho1, rho2',rho1sol(1,h),rho2(h)
!       write(*,*)'rho1, rho2',rho1sol(2,h),rho2(h)
          CALL poly_eval(pmod,ncoe_p,rho1sol(1,h),rho2sol(h),peval(1))  
          CALL poly_eval(pmod,ncoe_p,rho1sol(2,h),rho2sol(h),peval(2))  
!          write(*,*)'p(rho1_1,rho2)',peval(1)
!          write(*,*)'p(rho1_2,rho2)',peval(2)
          CALL poly_eval(qmod,ncoe_q,rho1sol(1,h),rho2sol(h),qeval(1))  
          CALL poly_eval(qmod,ncoe_q,rho1sol(2,h),rho2sol(h),qeval(2))  

113 FORMAT(A19,2X,3(F14.9,2X))
!          select only positive rho1
          IF(rho1sol(1,h).gt.0.d0.OR.rho1sol(2,h).gt.0.d0)THEN
             nposroots=nposroots+1
             rho2(nposroots) = rho2sol(h)
!             write(*,113)'almeno una positiva',rho2sol(h),rho1sol(1,h),rho1sol(2,h)
             
             IF(rho1sol(1,h).gt.0.d0)THEN
                rho1(nposroots) = rho1sol(1,h)
                IF(abs(peval(2)).lt.abs(peval(1)).AND. &
                     & rho1sol(2,h).gt.0.d0)THEN
                   rho1(nposroots)= rho1sol(2,h)
                ENDIF
             ELSE 
 !            ELSEIF(rho1sol(2,h).gt.0.d0)THEN
                rho1(nposroots) = rho1sol(2,h)
 !               IF(abs(peval(1)).lt.abs(peval(2)).AND. &
 !                    & rho1sol(1,h).gt.0.d0)THEN
 !                  rho1(nposroots)= rho1sol(1,h)
 !               ENDIF
             ENDIF
             
          ELSE
             IF(verb_2body.gt.10)THEN
                WRITE(*,*)'both negative solutions! rho1,rho2', &
                  & rho1sol(1,h),rho1sol(2,h)
             ENDIF
          ENDIF
          
       ELSE
          IF(verb_2body.gt.10)THEN
             write(*,*)'negative discriminant',discrim
          ENDIF
!          rho1sol(1,h) = 0.d0 
!          rho1sol(2,h) = 0.d0 
       ENDIF
    ENDDO

  END SUBROUTINE solvesys


! ==============================================================
  SUBROUTINE ZGETF2_QP( M, N, A, LDA, IPIV, INFO ) 
    IMPLICIT NONE
    INTEGER, PARAMETER :: qkind=kind(1.q0) !for quadruple precision
!                                                                       
!  -- LAPACK routine (version 2.0) --                                   
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,       
!     Courant Institute, Argonne National Lab, and Rice University      
!     September 30, 1994                                                
!                                                                       
! .. Scalar Arguments ..                                            
    INTEGER :: INFO, LDA, M, N 
! ..                                                                
! .. Array Arguments ..                                             
    INTEGER :: IPIV( * ) 
    COMPLEX(KIND=qkind) ::  A( LDA, * ) 
!  ..                                                                
!                                                                       
!  Purpose                                                              
!  =======                                                              
!                                                                       
!  ZGETF2 computes an LU factorization of a general m-by-n matrix A     
!  using partial pivoting with row interchanges.                        
!                                                                       
!  The factorization has the form                                       
!     A = P * L * U                                                     
!  where P is a permutation matrix, L is lower triangular with unit     
!  diagonal elements (lower trapezoidal if m > n), and U is upper       
!  triangular (upper trapezoidal if m < n).                             
!                                                                       
!  This is the right-looking Level 2 BLAS version of the algorithm.     
!                                                                       
!  Arguments                                                            
!  =========                                                            
!                                                                       
!  M       (input) INTEGER                                              
!          The number of rows of the matrix A.  M >= 0.                 
!                                                                       
!  N       (input) INTEGER                                              
!          The number of columns of the matrix A.  N >= 0.              
!                                                                       
!  A       (input/output) COMPLEX*16 array, dimension (LDA,N)           
!          On entry, the m by n matrix to be factored.                  
!          On exit, the factors L and U from the factorization          
!          A = P*L*U; the unit diagonal elements of L are not stored.   
!                                                                       
!  LDA     (input) INTEGER                                              
!          The leading dimension of the array A.  LDA >= max(1,M).      
!                                                                       
!  IPIV    (output) INTEGER array, dimension (min(M,N))                 
!          The pivot indices; for 1 <= i <= min(M,N), row i of the      
!          matrix was interchanged with row IPIV(i).                    
!                                                                       
!  INFO    (output) INTEGER                                             
!          = 0: successful exit                                         
!          < 0: if INFO = -k, the k-th argument had an illegal value    
!          > 0: if INFO = k, U(k,k) is exactly zero. The factorization  
!               has been completed, but the factor U is exactly         
!               singular, and division by zero will occur if it is used 
!               to solve a system of equations.                         
!                                                                       
!  =====================================================================
!                                                                       
!     .. Parameters ..                                                  
    COMPLEX(KIND=qkind), PARAMETER :: ONE = (1.q0, 0.q0), ZERO = (0.q0, 0.q0)
!     ..                                                                
!     .. Local Scalars ..                                               
    INTEGER            J, JP 
!     ..                                                                
!     .. External Functions ..                                          
    INTEGER            IZAMAX_QP 
    EXTERNAL           IZAMAX_QP 
!     ..                                                                
!     .. External Subroutines ..                                        
    EXTERNAL           XERBLA_QP, ZGERU_QP, ZSCAL_QP, ZSWAP_QP 
!     ..                                                                
!     .. Intrinsic Functions ..                                         
    INTRINSIC          MAX, MIN 
!     ..                                                                
!     .. Executable Statements ..                                       
!                                                                       
!     Test the input parameters.                                        
!                                                                       
    INFO = 0 
    IF( M.LT.0 ) THEN 
       INFO = -1 
    ELSE IF( N.LT.0 ) THEN 
       INFO = -2 
    ELSE IF( LDA.LT.MAX( 1, M ) ) THEN 
       INFO = -4 
    END IF
    IF( INFO.NE.0 ) THEN 
       CALL XERBLA_QP( 'ZGETF2', -INFO ) 
       RETURN 
    END IF
!                                                                       
!     Quick return if possible                                          
!                                                                       
    IF( M.EQ.0 .OR. N.EQ.0 )                                          &
         &   RETURN                                                         
!                                                                       
    DO 10 J = 1, MIN( M, N ) 
!                                                                       
!        Find pivot and test for singularity.                           
!                                                                       
       JP = J - 1 + IZAMAX_QP( M-J+1, A( J, J ), 1 ) 
       IPIV( J ) = JP 
       IF( A( JP, J ).NE.ZERO ) THEN 
!                                                                       
!           Apply the interchange to columns 1:N.                       
!                                                                       
          IF( JP.NE.J )                                               &
               &         CALL ZSWAP_QP( N, A( J, 1 ), LDA, A( JP, 1 ), LDA )         
!                                                                       
!           Compute elements J+1:M of J-th column.                      
!                                                                       
          IF( J.LT.M )                                                &
               &         CALL ZSCAL_QP( M-J, ONE / A( J, J ), A( J+1, J ), 1 )       
!                                                                       
       ELSE IF( INFO.EQ.0 ) THEN 
!                                                                       
          INFO = J 
       END IF
!                                                                       
       IF( J.LT.MIN( M, N ) ) THEN 
!                                                                       
!           Update trailing submatrix.                                  
!                                                                       
          CALL ZGERU_QP( M-J, N-J, -ONE, A( J+1, J ), 1, A( J, J+1 ),    &
               &                  LDA, A( J+1, J+1 ), LDA )                       
       END IF
10  END DO
    RETURN 
!                                                                       
!     End of ZGETF2                                                     
!                                                                       
  END SUBROUTINE ZGETF2_QP

! ==========================================================
  SUBROUTINE ZGERU_QP ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA ) 
    IMPLICIT NONE
    INTEGER, PARAMETER :: qkind=kind(1.q0) !for quadruple precision

!     .. Scalar Arguments ..                                            
      COMPLEX(KIND=qkind) :: ALPHA 
      INTEGER            INCX, INCY, LDA, M, N 
!     .. Array Arguments ..                                             
      COMPLEX(KIND=qkind) :: A( LDA, * ), X( * ), Y( * ) 
!     ..                                                                
!                                                                       
!  Purpose                                                              
!  =======                                                              
!                                                                       
!  ZGERU  performs the rank 1 operation                                 
!                                                                       
!     A := alpha*x*y' + A,                                              
!                                                                       
!  where alpha is a scalar, x is an m element vector, y is an n element 
!  vector and A is an m by n matrix.                                    
!                                                                       
!  Parameters                                                           
!  ==========                                                           
!                                                                       
!  M      - INTEGER.                                                    
!           On entry, M specifies the number of rows of the matrix A.   
!           M must be at least zero.                                    
!           Unchanged on exit.                                          
!                                                                       
!  N      - INTEGER.                                                    
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.                                    
!           Unchanged on exit.                                          
!                                                                       
!  ALPHA  - COMPLEX*16      .                                           
!           On entry, ALPHA specifies the scalar alpha.                 
!           Unchanged on exit.                                          
!                                                                       
!  X      - COMPLEX*16       array of dimension at least                
!           ( 1 + ( m - 1 )*abs( INCX ) ).                              
!           Before entry, the incremented array X must contain the m    
!           element vector x.                                           
!           Unchanged on exit.                                          
!                                                                       
!  INCX   - INTEGER.                                                    
!           On entry, INCX specifies the increment for the elements of  
!           X. INCX must not be zero.                                   
!           Unchanged on exit.                                          
!                                                                       
!  Y      - COMPLEX*16       array of dimension at least                
!           ( 1 + ( n - 1 )*abs( INCY ) ).                              
!           Before entry, the incremented array Y must contain the n    
!           element vector y.                                           
!           Unchanged on exit.                                          
!                                                                       
!  INCY   - INTEGER.                                                    
!           On entry, INCY specifies the increment for the elements of  
!           Y. INCY must not be zero.                                   
!           Unchanged on exit.                                          
!                                                                       
!  A      - COMPLEX*16       array of DIMENSION ( LDA, n ).             
!           Before entry, the leading m by n part of the array A must   
!           contain the matrix of coefficients. On exit, A is           
!           overwritten by the updated matrix.                          
!                                                                       
!  LDA    - INTEGER.                                                    
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least          
!           max( 1, m ).                                                
!           Unchanged on exit.                                          
!                                                                       
!                                                                       
!  Level 2 Blas routine.                                                
!                                                                       
!  -- Written on 22-October-1986.                                       
!     Jack Dongarra, Argonne National Lab.                              
!     Jeremy Du Croz, Nag Central Office.                               
!     Sven Hammarling, Nag Central Office.                              
!     Richard Hanson, Sandia National Labs.                             
!                                                                       
!                                                                       
!     .. Parameters ..                                                  
      COMPLEX(KIND=qkind), PARAMETER :: ZERO = (0.q0, 0.q0)  
!     .. Local Scalars ..                                               
      COMPLEX(KIND=qkind)         TEMP 
      INTEGER            I, INFO, IX, J, JY, KX 
!     .. External Subroutines ..                                        
      EXTERNAL           XERBLA_QP 
!     .. Intrinsic Functions ..                                         
      INTRINSIC          MAX 
!     ..                                                                
!     .. Executable Statements ..                                       
!                                                                       
!     Test the input parameters.                                        
!                                                                       
      INFO = 0 
      IF     ( M.LT.0 )THEN 
         INFO = 1 
      ELSE IF( N.LT.0 )THEN 
         INFO = 2 
      ELSE IF( INCX.EQ.0 )THEN 
         INFO = 5 
      ELSE IF( INCY.EQ.0 )THEN 
         INFO = 7 
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN 
         INFO = 9 
      END IF 
      IF( INFO.NE.0 )THEN 
         CALL XERBLA_QP( 'ZGERU ', INFO ) 
         RETURN 
      END IF 
!                                                                       
!     Quick return if possible.                                         
!                                                                       
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )               &
     &   RETURN                                                         
!                                                                       
!     Start the operations. In this version the elements of A are       
!     accessed sequentially with one pass through A.                    
!                                                                       
      IF( INCY.GT.0 )THEN 
         JY = 1 
      ELSE 
         JY = 1 - ( N - 1 )*INCY 
      END IF 
      IF( INCX.EQ.1 )THEN 
         DO 20, J = 1, N 
            IF( Y( JY ).NE.ZERO )THEN 
               TEMP = ALPHA*Y( JY ) 
               DO 10, I = 1, M 
                  A( I, J ) = A( I, J ) + X( I )*TEMP 
   10          CONTINUE 
            END IF 
            JY = JY + INCY 
   20    CONTINUE 
      ELSE 
         IF( INCX.GT.0 )THEN 
            KX = 1 
         ELSE 
            KX = 1 - ( M - 1 )*INCX 
         END IF 
         DO 40, J = 1, N 
            IF( Y( JY ).NE.ZERO )THEN 
               TEMP = ALPHA*Y( JY ) 
               IX   = KX 
               DO 30, I = 1, M 
                  A( I, J ) = A( I, J ) + X( IX )*TEMP 
                  IX        = IX        + INCX 
   30          CONTINUE 
            END IF 
            JY = JY + INCY 
   40    CONTINUE 
      END IF 
!                                                                       
      RETURN 
!                                                                       
!     End of ZGERU .                                                    
!                                                                       
    END SUBROUTINE ZGERU_QP 

! ===============================================================
      subroutine  zscal_QP(n,za,zx,incx) 
        IMPLICIT NONE
        INTEGER, PARAMETER :: qkind=kind(1.q0) !for quadruple precision

!                                                                       
!     scales a vector by a constant.                                    
!     jack dongarra, 3/11/78.                                           
!     modified 3/93 to return if incx .le. 0.                           
!     modified 12/3/93, array(1) declarations changed to array(*)       
!                                                                       
      complex(KIND=qkind) za,zx(*) 
      integer i,incx,ix,n 
!                                                                       
      if( n.le.0 .or. incx.le.0 )return 
      if(incx.eq.1)go to 20 
!                                                                       
!        code for increment not equal to 1                              
!                                                                       
      ix = 1 
      do 10 i = 1,n 
        zx(ix) = za*zx(ix) 
        ix = ix + incx 
   10 continue 
      return 
!                                                                       
!        code for increment equal to 1                                  
!                                                                       
   20 do 30 i = 1,n 
        zx(i) = za*zx(i) 
   30 continue 
      return 
    END subroutine zscal_QP

! ===============================================================
      subroutine  zswap_QP (n,zx,incx,zy,incy) 
        IMPLICIT NONE
        INTEGER, PARAMETER :: qkind=kind(1.q0) !for quadruple precision
!                                                                       
!     interchanges two vectors.                                         
!     jack dongarra, 3/11/78.                                           
!     modified 12/3/93, array(1) declarations changed to array(*)       
!                                                                       
      complex(KIND=qkind) zx(*),zy(*),ztemp 
      integer i,incx,incy,ix,iy,n 
!                                                                       
      if(n.le.0)return 
      if(incx.eq.1.and.incy.eq.1)go to 20 
!                                                                       
!       code for unequal increments or equal increments not equal       
!         to 1                                                          
!                                                                       
      ix = 1 
      iy = 1 
      if(incx.lt.0)ix = (-n+1)*incx + 1 
      if(incy.lt.0)iy = (-n+1)*incy + 1 
      do 10 i = 1,n 
        ztemp = zx(ix) 
        zx(ix) = zy(iy) 
        zy(iy) = ztemp 
        ix = ix + incx 
        iy = iy + incy 
   10 continue 
      return 
!                                                                       
!       code for both increments equal to 1                             
   20 do 30 i = 1,n 
        ztemp = zx(i) 
        zx(i) = zy(i) 
        zy(i) = ztemp 
   30 continue 
      return 
      END subroutine zswap_QP                                         

! ===============================================================
      SUBROUTINE XERBLA_QP( SRNAME, INFO ) 
        IMPLICIT NONE

!                                                                       
!  -- LAPACK auxiliary routine (version 2.0) --                         
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,       
!     Courant Institute, Argonne National Lab, and Rice University      
!     September 30, 1994                                                
!                                                                       
!     .. Scalar Arguments ..                                            
      CHARACTER*6        SRNAME 
      INTEGER            INFO 
!     ..                                                                
!                                                                       
!  Purpose                                                              
!  =======                                                              
!                                                                       
!  XERBLA  is an error handler for the LAPACK routines.                 
!  It is called by an LAPACK routine if an input parameter has an       
!  invalid value.  A message is printed and execution stops.            
!                                                                       
!  Installers may consider modifying the STOP statement in order to     
!  call system-specific exception-handling facilities.                  
!                                                                       
!  Arguments                                                            
!  =========                                                            
!                                                                       
!  SRNAME  (input) CHARACTER*6                                          
!          The name of the routine which called XERBLA.                 
!                                                                       
!  INFO    (input) INTEGER                                              
!          The position of the invalid parameter in the parameter list  
!          of the calling routine.                                      
!                                                                       
! ===================================================================== 
!                                                                       
!     .. Executable Statements ..                                       
!                                                                       
      WRITE( *, FMT = 9999 )SRNAME, INFO 
!                                                                       
      STOP 
!                                                                       
 9999 FORMAT( ' ** On entry to ', A6, ' parameter number ', I2, ' had ',&
     &      'an illegal value' )                                        
!                                                                       
!     End of XERBLA                                                     
!                                                                       
      END SUBROUTINE XERBLA_QP

! ===============================================================
      integer function izamax_QP(n,zx,incx) 
        IMPLICIT NONE
        INTEGER, PARAMETER :: qkind=kind(1.q0) !for quadruple precision

!                                                                       
!     finds the index of element having max. absolute value.            
!     jack dongarra, 1/15/85.                                           
!     modified 3/93 to return if incx .le. 0.                           
!     modified 12/3/93, array(1) declarations changed to array(*)       
!                                                                       
      complex(KIND=qkind) :: zx(*) 
      REAL(KIND=qkind) :: smax 
      integer i,incx,ix,n 
      REAL(KIND=qkind) :: dcabs1_QP 
!                                                                       
      izamax_QP = 0 
      if( n.lt.1 .or. incx.le.0 )return 
      izamax_QP = 1 
      if(n.eq.1)return 
      if(incx.eq.1)go to 20 
!                                                                       
!        code for increment not equal to 1                              
!                                                                       
      ix = 1 
      smax = dcabs1_QP(zx(1)) 
      ix = ix + incx 
      do 10 i = 2,n 
         if(dcabs1_QP(zx(ix)).le.smax) go to 5 
         izamax_QP = i 
         smax = dcabs1_QP(zx(ix)) 
    5    ix = ix + incx 
   10 continue 
      return 
!                                                                       
!        code for increment equal to 1                                  
!                                                                       
   20 smax = dcabs1_QP(zx(1)) 
      do 30 i = 2,n 
         if(dcabs1_QP(zx(i)).le.smax) go to 30 
         izamax_QP = i 
         smax = dcabs1_QP(zx(i)) 
   30 continue 
      return 
      END function izamax_QP

! ===============================================================
      REAL(KIND=16) function dcabs1_QP(z) 
        IMPLICIT NONE
        INTEGER, PARAMETER :: qkind=kind(1.q0) !for quadruple precision
       complex(KIND=qkind) :: z 
!      complex*16 z,zz                                                  
!      double precision t(2)                                            
!      equivalence (zz,t(1))                                            
!      zz = z                                                           
!      dcabs1 = dabs(t(1)) + dabs(t(2))                                 

!       dcabs1_QP=abs(real(z,qkind))+abs(real(z*(0,1),qkind))
       dcabs1_QP=abs(real(z))+abs(real(z*(0,1)))                          

!       dcabs1=abs(DBLE(z))+abs(DBLE(z*(0,1))) 
                                                                        
       return 
     END function dcabs1_QP

!END MODULE detcomp

! ============================================================
! ******* CODING of the n evaluations of a polynomial ********
! ********* p: p(om(1)), p(om(n)) to apply irvfft ************
! ******** written by GIOVANNI F. GRONCHI (2004) *************
! ============================================================
  SUBROUTINE code_input(N,eval,codeval)
    IMPLICIT NONE
    INTEGER :: N
! evaluation of poly p(x) in the N-th roots of unity
    COMPLEX(KIND=8),DIMENSION(N) :: eval
! -----------------------------------------------------------
! if omega is a primitive complex root of unity then we have
! Re(p(omega^s)) = Re(p(omeg^{-s}))
! Im(p(omega^s)) = -Im(p(omeg^{-s}))   for s=1,..,N/2
! -----------------------------------------------------------
    REAL(KIND=8),DIMENSION(N) :: codeval ! coded evaluations
! ------ end interface ----------------
    INTEGER :: i ! loop indexes
! ============================================================

    DO i = 1,N/2+1
       codeval(i) = DBLE(eval(i)) ! from 1 to N/2+1 take the real part
    ENDDO
! (roots of unity are taken clockwise)
    DO i = 1,N/2-1 
       codeval(N/2+1+i) = DIMAG(eval(N/2-i+1))
    ENDDO
  END SUBROUTINE code_input

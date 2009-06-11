
! Given the multi-index coefficients of two bivariate polinomials 
! calpha,dbeta of degree amod, bmod respectively, computes the
! degree and coefficients of the polynomial sum.

SUBROUTINE poly_sum(amod,bmod,calpha,dbeta,gmod,egamma)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: amod,bmod ! modulus of the multi-index
  DOUBLE PRECISION,INTENT(IN) :: calpha(0:amod,0:amod)
  DOUBLE PRECISION,INTENT(IN) :: dbeta(0:bmod,0:bmod)
  INTEGER,INTENT(OUT) :: gmod ! modulus of the multi-index
!  DOUBLE PRECISION,INTENT(OUT) :: egamma(0:amod+bmod,0:amod+bmod)
  DOUBLE PRECISION :: egamma(0:max(amod,bmod),0:max(amod,bmod))
! ==== end interface ======
  ! auxiliary
  DOUBLE PRECISION :: calpha_aux(0:max(amod,bmod),0:max(amod,bmod))
  DOUBLE PRECISION :: dbeta_aux(0:max(amod,bmod),0:max(amod,bmod))
  INTEGER :: m,n,h,k   ! loop index

  egamma = 0.d0
  calpha_aux=0.d0
  dbeta_aux=0.d0

  calpha_aux(0:amod,0:amod) =  calpha(0:amod,0:amod)
  dbeta_aux(0:bmod,0:bmod) =  dbeta(0:bmod,0:bmod)
  gmod = max(amod,bmod)
  DO h=0,gmod
     DO k = 0,gmod-h
        egamma(h,k)= calpha_aux(h,k)+dbeta_aux(h,k)
     ENDDO
  ENDDO

END SUBROUTINE poly_sum


! Given the multi-index coefficients of two bivariate polinomials 
! calpha,dbeta of degree amod, bmod respectively, computes the
! degree and coefficients of the product polynomial according to the
! Cauchy rule

SUBROUTINE poly_product(amod,bmod,calpha,dbeta,gmod,egamma)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: amod,bmod ! modulus of the multi-index
  DOUBLE PRECISION,INTENT(IN) :: calpha(0:amod,0:amod)
  DOUBLE PRECISION,INTENT(IN) :: dbeta(0:bmod,0:bmod)
  INTEGER,INTENT(OUT) :: gmod ! modulus of the multi-index
  DOUBLE PRECISION,INTENT(OUT) :: egamma(0:amod+bmod,0:amod+bmod)
! ==== end interface ======
  INTEGER :: m,n,h,k   ! loop index

  egamma = 0.d0

  gmod = amod+bmod

  DO h=0,amod
     DO k = 0,amod-h
        DO m = 0,bmod
           DO n = 0,bmod-m
              egamma(h+m,k+n)=egamma(h+m,k+n)+calpha(h,k)*dbeta(m,n)
           ENDDO
        ENDDO
     ENDDO
  ENDDO


END SUBROUTINE poly_product

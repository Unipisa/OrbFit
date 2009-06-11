!*** DOUBLE PRECISION ***
SUBROUTINE poly_eval(pmod,polycoe,rho1,rho2,peval)
  USE fund_const, ONLY: dkind
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: pmod
  REAL(KIND=dkind),INTENT(IN) :: polycoe(0:pmod,0:pmod) !these are qkind
  REAL(KIND=dkind),INTENT(IN) :: rho1,rho2
  REAL(KIND=dkind) ,INTENT(OUT) :: peval 
! ==== end interface ====
  INTEGER :: h,k ! loop indexes

  peval = 0.q0
  DO h=0,pmod
     DO k=0,pmod      
        IF(polycoe(h,k).ne.0.q0)THEN
           peval = peval + polycoe(h,k)*rho1**h*rho2**k
        ENDIF
     ENDDO
  ENDDO

!  WRITE(*,*)'EVALUATION=',peval

END SUBROUTINE poly_eval

!*** QUADRUPLE PRECISION ***
SUBROUTINE poly_eval_QP(pmod,polycoe,rho1,rho2,peval)
  USE fund_const, ONLY: qkind
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: pmod
  REAL(KIND=qkind),INTENT(IN) :: polycoe(0:pmod,0:pmod)
  REAL(KIND=qkind),INTENT(IN) :: rho1,rho2
  REAL(KIND=qkind) ,INTENT(OUT) :: peval 
! ==== end interface ====
  INTEGER :: h,k ! loop indexes

  peval = 0.q0
  DO h=0,pmod
     DO k=0,pmod      
        IF(polycoe(h,k).ne.0.q0)THEN
           peval = peval + polycoe(h,k)*rho1**h*rho2**k
        ENDIF
     ENDDO
  ENDDO

!  WRITE(*,*)'EVALUATION=',peval

END SUBROUTINE poly_eval_QP

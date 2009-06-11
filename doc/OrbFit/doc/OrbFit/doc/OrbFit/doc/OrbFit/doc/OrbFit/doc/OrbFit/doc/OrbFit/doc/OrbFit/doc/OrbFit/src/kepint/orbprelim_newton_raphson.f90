! =================================================================
SUBROUTINE orbprelim_newton_raphson(pmod,qmod,coe_p,coe_q,rho1,rho2, &
     & pk,qk)
  USE fund_const, ONLY: dkind
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: pmod,qmod ! modulus of the multi-index
  REAL(KIND=dkind),INTENT(IN) :: coe_p(0:pmod,0:pmod),coe_q(0:qmod,0:qmod)
!
  REAL(KIND=dkind),INTENT(INOUT) :: rho1,rho2
  REAL(KIND=dkind),INTENT(OUT) :: pk,qk
! ----- end interface ----
  REAL(KIND=dkind) :: r1k,r2k   ! auxiliary
  REAL(KIND=dkind) :: coe_p_r1(0:pmod-1,0:pmod-1),coe_p_r2(0:pmod-1,0:pmod-1)
  REAL(KIND=dkind) :: coe_q_r1(0:qmod-1,0:qmod-1),coe_q_r2(0:qmod-1,0:qmod-1)
  REAL(KIND=dkind) :: pkr1,pkr2,qkr1,qkr2  ! auxiliary
  REAL(KIND=dkind) :: detJ
  REAL(KIND=dkind),DIMENSION(2,2) :: Jinv
! ====================================================================

  r1k=rho1
  r2k=rho2

! compute values of p,q in r1k,r2k
  CALL poly_eval(pmod,coe_p,r1k,r2k,pk)  
  CALL poly_eval(qmod,coe_q,r1k,r2k,qk)  

! compute Jacobian matrix J = der(p,q)/der(rho1,rho2)
  CALL pq_derpar(pmod,qmod,coe_p,coe_q,coe_p_r1,coe_p_r2,coe_q_r1,coe_q_r2)

! coefficients of J
  CALL poly_eval(pmod-1,coe_p_r1,r1k,r2k,pkr1)  
  CALL poly_eval(pmod-1,coe_p_r2,r1k,r2k,pkr2)  
  CALL poly_eval(qmod-1,coe_q_r1,r1k,r2k,qkr1)  
  CALL poly_eval(qmod-1,coe_q_r2,r1k,r2k,qkr2)  

! determinant  of J
  detJ = pkr1*qkr2-qkr1*pkr2
!  write(*,*)'detJ',detJ

! inverse of J
  Jinv(1,1) =  qkr2/detJ
  Jinv(1,2) = -pkr2/detJ
  Jinv(2,1) = -qkr1/detJ
  Jinv(2,2) =  pkr1/detJ

! iteration of Newton-Raphson
  rho1 = r1k - Jinv(1,1)*pk - Jinv(1,2)*qk
  rho2 = r2k - Jinv(2,1)*pk - Jinv(2,2)*qk

END SUBROUTINE orbprelim_newton_raphson

! =================================================================
SUBROUTINE orbprelim_newton_raphson_QP(pmod,qmod,coe_p,coe_q,rho1,rho2, &
     & pk,qk)
  USE fund_const, ONLY: qkind
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: pmod,qmod ! modulus of the multi-index
  REAL(KIND=qkind),INTENT(IN) :: coe_p(0:pmod,0:pmod),coe_q(0:qmod,0:qmod)
!
  REAL(KIND=qkind),INTENT(INOUT) :: rho1,rho2
  REAL(KIND=qkind),INTENT(OUT) :: pk,qk
! ----- end interface ----
  REAL(KIND=qkind) :: r1k,r2k   ! auxiliary
  REAL(KIND=qkind) :: coe_p_r1(0:pmod-1,0:pmod-1),coe_p_r2(0:pmod-1,0:pmod-1)
  REAL(KIND=qkind) :: coe_q_r1(0:qmod-1,0:qmod-1),coe_q_r2(0:qmod-1,0:qmod-1)
  REAL(KIND=qkind) :: pkr1,pkr2,qkr1,qkr2  ! auxiliary
  REAL(KIND=qkind) :: detJ
  REAL(KIND=qkind),DIMENSION(2,2) :: Jinv
! ====================================================================

  r1k=rho1
  r2k=rho2

! compute values of p,q in r1k,r2k
  CALL poly_eval_QP(pmod,coe_p,r1k,r2k,pk)  
  CALL poly_eval_QP(qmod,coe_q,r1k,r2k,qk)  

! compute Jacobian matrix J = der(p,q)/der(rho1,rho2)
  CALL pq_derpar_QP(pmod,qmod,coe_p,coe_q,coe_p_r1,coe_p_r2,coe_q_r1,coe_q_r2)

! coefficients of J
  CALL poly_eval_QP(pmod-1,coe_p_r1,r1k,r2k,pkr1)  
  CALL poly_eval_QP(pmod-1,coe_p_r2,r1k,r2k,pkr2)  
  CALL poly_eval_QP(qmod-1,coe_q_r1,r1k,r2k,qkr1)  
  CALL poly_eval_QP(qmod-1,coe_q_r2,r1k,r2k,qkr2)  

! determinant  of J
  detJ = pkr1*qkr2-qkr1*pkr2
!  write(*,*)'detJ',detJ

! inverse of J
  Jinv(1,1) =  qkr2/detJ
  Jinv(1,2) = -pkr2/detJ
  Jinv(2,1) = -qkr1/detJ
  Jinv(2,2) =  pkr1/detJ

! iteration of Newton-Raphson
  rho1 = r1k - Jinv(1,1)*pk - Jinv(1,2)*qk
  rho2 = r2k - Jinv(2,1)*pk - Jinv(2,2)*qk

END SUBROUTINE orbprelim_newton_raphson_QP

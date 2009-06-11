
PROGRAM test_pq_derpar
  IMPLICIT NONE
  INTEGER, PARAMETER :: qkind=kind(1.q0) !for quadruple precision
  INTEGER,PARAMETER :: pmod=24,qmod=2
  REAL(KIND=qkind) ::  coe_p(0:pmod,0:pmod),coe_q(0:qmod,0:qmod)

  REAL(KIND=qkind) :: coe_p_r1(0:pmod-1,0:pmod-1)
  REAL(KIND=qkind) :: coe_p_r2(0:pmod-1,0:pmod-1)
  REAL(KIND=qkind) :: coe_q_r1(0:qmod-1,0:qmod-1)
  REAL(KIND=qkind) :: coe_q_r2(0:qmod-1,0:qmod-1)
  INTEGER :: i,j ! loop indexes

  coe_p=0.d0;coe_q=0.d0

  coe_p(20,0)= 1.d0
  coe_p(4,20)= 2.d0
  coe_p(10,10)= -3.d0
  coe_p(5,0)= 5.d0
  coe_p(0,3)= -2.d0

  coe_q(0,0)= 0.d0
  coe_q(0,1)= -1.d0
  coe_q(0,2)= 5.d0
  coe_q(1,0)= 0.d0
  coe_q(1,1)= 2.d0
  coe_q(2,0)= 3.d0

  CALL pq_derpar(pmod,qmod,coe_p,coe_q,coe_p_r1,coe_p_r2,coe_q_r1,coe_q_r2)

  DO i=0,pmod-1
     DO j=0,pmod-i-1
        IF(coe_p_r1(i,j).ne.0.d0)THEN
           WRITE(*,*)'coe_p_r1(',i,j,')=',coe_p_r1(i,j)
        ENDIF
        IF(coe_p_r2(i,j).ne.0.d0)THEN
           WRITE(*,*)'coe_p_r2(',i,j,')=',coe_p_r2(i,j)
        ENDIF
     ENDDO
  ENDDO

  DO i=0,qmod-1
     DO j=0,qmod-i-1
        IF(coe_q_r1(i,j).ne.0.d0)THEN
           WRITE(*,*)'coe_q_r1(',i,j,')=',coe_q_r1(i,j)
        ENDIF
        IF(coe_q_r2(i,j).ne.0.d0)THEN        
           WRITE(*,*)'coe_q_r2(',i,j,')=',coe_q_r2(i,j)
        ENDIF
     ENDDO
  ENDDO

END PROGRAM test_pq_derpar

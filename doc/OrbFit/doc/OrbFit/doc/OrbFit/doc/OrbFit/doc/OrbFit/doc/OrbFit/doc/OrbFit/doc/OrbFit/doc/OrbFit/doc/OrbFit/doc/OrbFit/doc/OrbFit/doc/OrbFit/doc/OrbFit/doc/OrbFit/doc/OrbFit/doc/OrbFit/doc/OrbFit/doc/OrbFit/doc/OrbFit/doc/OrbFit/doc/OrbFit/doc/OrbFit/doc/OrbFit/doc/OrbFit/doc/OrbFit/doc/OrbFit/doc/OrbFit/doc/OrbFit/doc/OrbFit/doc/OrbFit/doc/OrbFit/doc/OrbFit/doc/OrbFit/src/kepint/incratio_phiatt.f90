SUBROUTINE incratio_phiatt(att1,att2,rrd1,rrd2,qqd1,qqd2,incr_phiatt)
  USE attributable
  IMPLICIT NONE
  !INTEGER, PARAMETER :: qkind=kind(1.q0) !for quadruple precision
  TYPE(attrib),INTENT(IN) :: att1,att2
  REAL(KIND=qkind), INTENT(IN), DIMENSION(2) :: rrd1, rrd2 !radial dist/vel
  REAL(KIND=qkind), INTENT(IN), DIMENSION(3,2) :: qqd1,qqd2 !obs pos/vel
  
  REAL(KIND=qkind), INTENT(OUT), DIMENSION(4,8) :: incr_phiatt
  !end interface

  REAL(KIND=qkind), DIMENSION(8) :: att, att_prime ! att1, att2
  REAL(KIND=qkind), DIMENSION(4) :: rrd ! rrd1, rrd2
  REAL(KIND=qkind),DIMENSION(3,2) :: rhohat1,rhohat1_a,rhohat1_d
  REAL(KIND=qkind),DIMENSION(3,2) :: rhohat1,rhohat1_a,rhohat1_d
  REAL(KIND=qkind), DIMENSION(4) :: phi1,phi2
  REAL(KIND=qkind), DIMENSION(3,2) :: Dvec1,Evec1,Fvec1,Gvec1
  REAL(KIND=qkind), DIMENSION(3,2) :: Dvec2,Evec2,Fvec2,Gvec2 
  REAL(KIND=qkind),DIMENSION(2) :: c0_pr,c1_pr,c2_pr,c3_pr,c31_pr,c32_pr,c4_pr,c5_pr
  REAL(KIND=qkind),DIMENSION(2) :: c0_sec,c1_sec,c2_sec,c3_sec,c31_sec,c32_sec
  REAL(KIND=qkind), DIMENSION(2) :: c4_sec,c5_sec 
  REAL(KIND=qkind),DIMENSION(3,2) :: pe,pde
  REAL(KIND=qkind),DIMENSION(2) :: alpha1,delta1,ad1,dd1
  REAL(KIND=qkind),DIMENSION(2) :: alpha2,delta2,ad2,dd2
  REAL(KIND=qkind),DIMENSION(2) :: cosa1,sina1,cosd1,sind1
  REAL(KIND=qkind),DIMENSION(2) :: cosa2,sina2,cosd2,sind2
  REAL(KIND=qkind),DIMENSION(3,2) :: rhohat1,rhohat1_a,rhohat1_d
  REAL(KIND=qkind),DIMENSION(3,2) :: rhohat2,rhohat2_a,rhohat2_d
  REAL(KIND=qkind),DIMENSION(3,2) :: pe_pde,pe_rhat1,pe_rhat2
  REAL(KIND=qkind),DIMENSION(3,2) :: pe_rhat1a,pe_rhat1d,rhat1_rhat1a,rhat1_rhat1d
  REAL(KIND=qkind),DIMENSION(3,2) :: pe_rhat2a,pe_rhat2d,rhat2_rhat2a,rhat2_rhat2d
  REAL(KIND=qkind),DIMENSION(3,2) :: rhat1_pde, rhat2_pde

  REAL(KIND=qkind) :: incr
  INTEGER :: i,j,k
  
  att=
  att_prime=

  rrd(1:2)=rrd1
  rrd(3:4)=rrd2

  pe(1:3,1)= qqd1(1:3,1)
  pde(1:3,1)=qqd1(1:3,2)
  pe(1:3,2)=qqd2(1:3,1)
  pde(1:3,2)=qqd2(1:3,2)

  alpha1(1) = att(1)
  delta1(1) = att(2)
  ad1(1) = att(3)
  dd1(1) = att(4)
  alpha1(2) = att(5)
  delta1(2) = att(6)
  ad1(2) = att(7)
  dd1(2) = att(8)

  alpha2(1) = att_prime(1)
  delta2(1) = att_prime(2)
  ad2(1) = att_prime(3)
  dd2(1) = att_prime(4)
  alpha2(2) = att_prime(5)
  delta2(2) = att_prime(6)
  ad2(2) = att_prime(7)
  dd2(2) = att_prime(8)

  DO i=1,2
     cosa1(i)=qcos(alpha1(i))
     sina1(i)=qsin(alpha1(i))
     cosd1(i)=qcos(delta1(i))
     sind1(i)=qsin(delta1(i))
  ENDDO

  DO i=1,2
     cosa2(i)=qcos(alpha2(i))
     sina2(i)=qsin(alpha2(i))
     cosd2(i)=qcos(delta2(i))
     sind2(i)=qsin(delta2(i))
  ENDDO
  
  DO i=1,2
     rhohat1(1,i) = cosa1(i)*cosd1(i)
     rhohat1(2,i) = sina1(i)*cosd1(i)
     rhohat1(3,i) = sind1(i)
     rhohat1_a(1,i) = -sina1(i)*cosd1(i)
     rhohat1_a(2,i) = cosa1(i)*cosd1(i)
     rhohat1_a(3,i) = 0.q0
     rhohat1_d(1,i) = -cosa1(i)*sind1(i)
     rhohat1_d(2,i) = -sina1(i)*sind1(i)
     rhohat1_d(3,i) = cosd1(i)
  ENDDO

  DO i=1,2
     rhohat2(1,i) = cosa2(i)*cosd2(i)
     rhohat2(2,i) = sina2(i)*cosd2(i)
     rhohat2(3,i) = sind2(i)
     rhohat2_a(1,i) = -sina2(i)*cosd2(i)
     rhohat2_a(2,i) = cosa2(i)*cosd2(i)
     rhohat2_a(3,i) = 0.q0
     rhohat2_d(1,i) = -cosa2(i)*sind2(i)
     rhohat2_d(2,i) = -sina2(i)*sind2(i)
     rhohat2_d(3,i) = cosd2(i)
  ENDDO
  
  DO j=1,2
     CALL cross_prod_QP(pe(1:3,j),pde(1:3,j),pe_pde(1:3,j))
     CALL cross_prod_QP(pe(1:3,j),rhohat1_a(1:3,j),pe_rhat1a(1:3,j))
     CALL cross_prod_QP(pe(1:3,j),rhohat1_d(1:3,j),pe_rhat1d(1:3,j))
     CALL cross_prod_QP(rhohat1(1:3,j),pde(1:3,j),rhat1_pde(1:3,j))
     CALL cross_prod_QP(rhohat1(1:3,j),rhohat1_a(1:3,j),rhat1_rhat1a(1:3,j))
     CALL cross_prod_QP(rhohat1(1:3,j),rhohat1_d(1:3,j),rhat1_rhat1d(1:3,j))
     CALL cross_prod_QP(pe(1:3,j),rhohat1(1:3,j),pe_rhat1(1:3,j))
  ENDDO

  DO j=1,2
     CALL cross_prod_QP(pe(1:3,j),pde(1:3,j),pe_pde(1:3,j))
     CALL cross_prod_QP(pe(1:3,j),rhohat2_a(1:3,j),pe_rhat2a(1:3,j))
     CALL cross_prod_QP(pe(1:3,j),rhohat2_d(1:3,j),pe_rhat2d(1:3,j))
     CALL cross_prod_QP(rhohat2(1:3,j),pde(1:3,j),rhat2_pde(1:3,j))
     CALL cross_prod_QP(rhohat2(1:3,j),rhohat2_a(1:3,j),rhat2_rhat2a(1:3,j))
     CALL cross_prod_QP(rhohat2(1:3,j),rhohat2_d(1:3,j),rhat2_rhat2d(1:3,j))
     CALL cross_prod_QP(pe(1:3,j),rhohat2(1:3,j),pe_rhat2(1:3,j))
  ENDDO

! Admissible Region defining coefficients
  DO j=1,2
     c0_pr(j) = DOT_PRODUCT(pe(1:3,j),pe(1:3,j))
     c1_pr(j) = 2.q0*DOT_PRODUCT(pde(1:3,j),rhohat1(1:3,j))
     c2_pr(j) = (ad1(j)*cos(delta1(j)))**2 + dd1(j)**2
     c31_pr(j) = 2.q0*DOT_PRODUCT(pde(1:3,j),rhohat1_a(1:3,j))
     c32_pr(j) = 2.q0*DOT_PRODUCT(pde(1:3,j),rhohat1_d(1:3,j))
     c3_pr(j) = ad1(j)*c31_pr(j) + dd1(j)*c32_pr(j);
     c4_pr(j) = DOT_PRODUCT(pde(1:3,j),pde(1:3,j))
     c5_pr(j) = 2.q0*DOT_PRODUCT(pe(1:3,j),rhohat1(1:3,j))
  ENDDO

  DO j=1,2
     c0_sec(j) = DOT_PRODUCT(pe(1:3,j),pe(1:3,j))
     c1_sec(j) = 2.q0*DOT_PRODUCT(pde(1:3,j),rhohat2(1:3,j))
     c2_sec(j) = (ad2(j)*cos(delta2(j)))**2 + dd2(j)**2
     c31_sec(j) = 2.q0*DOT_PRODUCT(pde(1:3,j),rhohat2_a(1:3,j))
     c32_sec(j) = 2.q0*DOT_PRODUCT(pde(1:3,j),rhohat2_d(1:3,j))
     c3_sec(j) = ad2(j)*c31_sec(j) + dd2(j)*c32_sec(j);
     c4_sec(j) = DOT_PRODUCT(pde(1:3,j),pde(1:3,j))
     c5_sec(j) = 2.q0*DOT_PRODUCT(pe(1:3,j),rhohat2(1:3,j))
  ENDDO

! Angular Momentum curve
  DO j = 1,2
     Gvec1(1:3,j)= pe_pde(1:3,j)
     Fvec1(1:3,j)= ad1(j)*pe_rhat1a(1:3,j)+dd1(j)*pe_rhat1d(1:3,j)+rhat1_pde(1:3,j)
     Evec1(1:3,j)= ad1(j)*rhat1_rhat1a(1:3,j)+dd1(j)*rhat1_rhat1d(1:3,j)
     Dvec1(1:3,j)= pe_rhat1(1:3,j)
  ENDDO
  DO j = 1,2
     Gvec2(1:3,j)= pe_pde(1:3,j)
     Fvec2(1:3,j)= ad2(j)*pe_rhat2a(1:3,j)+dd2(j)*pe_rhat2d(1:3,j)+rhat2_pde(1:3,j)
     Evec2(1:3,j)= ad2(j)*rhat2_rhat2a(1:3,j)+dd2(j)*rhat2_rhat2d(1:3,j)
     Dvec2(1:3,j)= pe_rhat2(1:3,j)
  ENDDO
 

  phi1(1:3)=Dvec1(1:3,1)*rrd(2)-Dvec1(1:3,1)*rrd(4)-Evec1(1:3,2)*(rrd(3)**2)+&
       & Evec1(1:3,1)*(rrd(1)**2)-Fvec1(1:3,2)*rrd(3)+Fvec1(1:3,1)*rrd(1)- &
       & Gvec1(1:3,2)+Gvec1(1:3,1)
  phi2(1:3)=Dvec2(1:3,1)*rrd(2)-Dvec2(1:3,1)*rrd(4)-Evec2(1:3,2)*(rrd(3)**2)+&
       & Evec2(1:3,1)*(rrd(1)**2)-Fvec2(1:3,2)*rrd(3)+Fvec2(1:3,1)*rrd(1)- &
       & Gvec2(1:3,2)+Gvec2(1:3,1)
  phi1(4)=rrd(2)**2+c1_pr(1)*rrd(2)+c2_pr(1)*(rrd(1)**2)+c3_pr(1)*rrd(1)+c4_pr(1)-&
       & 2*gms*(rrd(1)**2+c5_pr(1)*rrd(1)+c0_pr(1))**(-1/2)-&
       & (rrd(4)**2+c1_pr(2)*rrd(4)+c2_pr(2)*(rrd(3)**2)+c3_pr(2)*rrd(3)+c4_pr(2)-&
       & 2*gms*(rrd(2)**2+c5_pr(2)*rrd(3)+c0_pr(2))**(-1/2))
  phi2(4)=rrd(2)**2+c1_sec(1)*rrd(2)+c2_sec(1)*(rrd(1)**2)+c3_sec(1)*rrd(1)+c4_sec(1)-&
       & 2*gms*(rrd(1)**2+c5_sec(1)*rrd(1)+c0_sec(1))**(-1/2)-&
       & (rrd(4)**2+c1_sec(2)*rrd(4)+c2_sec(2)*(rrd(3)**2)+c3_sec(2)*rrd(3)+c4_sec(2)-&
       & 2*gms*(rrd(2)**2+c5_sec(2)*rrd(3)+c0_sec(2))**(-1/2))
  DO i=1,4
     DO j=1,8
        incr_phiatt(i,j)=(phi1(i)-phi2(i))/(att(j)-att_prime(j))
     ENDDO
  ENDDO
END SUBROUTINE incratio_phiatt

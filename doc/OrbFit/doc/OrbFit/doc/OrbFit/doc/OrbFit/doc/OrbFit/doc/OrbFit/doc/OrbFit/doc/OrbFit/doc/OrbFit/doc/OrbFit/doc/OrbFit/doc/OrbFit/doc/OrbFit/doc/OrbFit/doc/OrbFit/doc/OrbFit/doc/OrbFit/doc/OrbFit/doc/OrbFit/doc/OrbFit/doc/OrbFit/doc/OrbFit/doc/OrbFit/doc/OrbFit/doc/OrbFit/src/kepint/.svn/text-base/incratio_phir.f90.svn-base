SUBROUTINE incratio_phir(att1,att2,rrd1,rrd2,qqd1,qqd2,incr_phir)
  USE attributable
  IMPLICIT NONE
  !  INTEGER, PARAMETER :: qkind=kind(1.q0) !for quadruple precision
  
  TYPE(attrib),INTENT(IN) :: att1,att2
  REAL(KIND=qkind), INTENT(IN), DIMENSION(2) :: rrd1, rrd2 !radial dist/vel
  REAL(KIND=qkind), INTENT(IN), DIMENSION(3,2) :: qqd1,qqd2 !obs pos/vel
  
  REAL(KIND=qkind), INTENT(OUT), DIMENSION(4,4) :: incr_phir
  !end interface
  REAL(KIND=qkind), DIMENSION(4) :: r1,r2 !rho1,rhod1,rho2,rhod2 
  REAL(KIND=qkind), DIMENSION(4) :: phi1,phi2
  REAL(KIND=qkind) :: incr
  INTEGER :: i,j
  !REAL(KIND=qkind), DIMENSION(3,2) :: Dvec,Evec,Fvec,Gvec !ang mom coeffs
  !REAL(KIND=qkind),DIMENSION(2) :: c0,c1,c2,c3,c31,c32,c4,c5 !energy coeffs

  r1= 
  r2=
  phi1(1:3)=Dvec(1:3,1)*r1(2)-Dvec(1:3,1)*r1(4)-Evec(1:3,2)*(r1(3)**2)+&
       & Evec(1:3,1)*(r1(1)**2)-Fvec(1:3,2)*r1(3)+Fvec(1:3,1)*r1(1)- &
       & Gvec(1:3,2)+Gvec(1:3,1)
  phi2(1:3)=Dvec(1:3,1)*r2(2)-Dvec(1:3,1)*r2(4)-Evec(1:3,2)*(r2(3)**2)+&
       & Evec(1:3,1)*(r2(1)**2)-Fvec(1:3,2)*r2(3)+Fvec(1:3,1)*r1(1)- &
       & Gvec(1:3,2)+Gvec(1:3,1)
  phi1(4)=r1(2)**2+c1(1)*r1(2)+c2(1)*(r1(1)**2)+c3(1)*r1(1)+c4(1)-&
       & 2*gms*(r1(1)**2+c5(1)*r1(1)+c0(1))**(-1/2)-&
       & (r1(4)**2+c1(2)*r1(4)+c2(2)*(r1(3)**2)+c3(2)*r1(3)+c4(2)-&
       & 2*gms*(r1(2)**2+c5(2)*r1(3)+c0(2))**(-1/2))
  phi2(4)=r2(2)**2+c1(1)*r2(2)+c2(1)*(r2(1)**2)+c3(1)*r2(1)+c4(1)-&
       & 2*gms*(r2(1)**2+c5(1)*r2(1)+c0(1))**(-1/2)-&
       & (r2(4)**2+c1(2)*r2(4)+c2(2)*(r2(3)**2)+c3(2)*r2(3)+c4(2)-&
       & 2*gms*(r2(2)**2+c5(2)*r2(3)+c0(2))**(-1/2))
  DO i=1,4
     DO j=1,4
        incr_phir(i,j)=(phi1(i)-phi2(i))/(r1(j)-r2(j))
     ENDDO
  ENDDO
END SUBROUTINE incratio_phir

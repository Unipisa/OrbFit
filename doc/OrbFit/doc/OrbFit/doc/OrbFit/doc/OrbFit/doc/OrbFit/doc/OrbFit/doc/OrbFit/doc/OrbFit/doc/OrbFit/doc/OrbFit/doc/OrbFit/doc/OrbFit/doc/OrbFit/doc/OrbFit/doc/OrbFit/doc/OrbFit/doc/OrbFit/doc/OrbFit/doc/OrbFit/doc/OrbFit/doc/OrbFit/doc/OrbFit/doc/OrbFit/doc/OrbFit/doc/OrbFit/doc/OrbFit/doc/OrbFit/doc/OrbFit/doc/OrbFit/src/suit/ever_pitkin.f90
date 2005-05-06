MODULE ever_pitkin
USE fund_const  
IMPLICIT NONE

PRIVATE

PUBLIC fser_propag, fser_propag_der, solve_peri

PUBLIC s_funct, r_of_psi !maybe to be made private later

CONTAINS
  SUBROUTINE fser_propag(x0,y0,t0,t,mu,x,y)
    DOUBLE PRECISION, DIMENSION(3), INTENT(IN):: x0,y0
    DOUBLE PRECISION, INTENT(IN):: t,t0,mu
    DOUBLE PRECISION, DIMENSION(3), INTENT(OUT):: x,y
    DOUBLE PRECISION :: s0,s1,s2,s3
    DOUBLE PRECISION r0, v0, vsize, prscal, alpha, sig0
    DOUBLE PRECISION psi, r,f,g,fdot,gdot
    DOUBLE PRECISION eps,period, dt
    r0=vsize(x0)
    sig0=prscal(x0,y0)
    v0=vsize(y0)
    alpha=v0**2-2*mu/r0 ! 2*E
    eps=100*epsilon(1.d0)
!   IF(alpha.lt.-eps)THEN
!       period=dpig*mu/(-alpha)**(1.5d0)
!      write(*,*)' period ', period
!      dt=t-t0-anint((t-t0)/period)*period
!     ELSE
      dt=t-t0
!     ENDIF
    CALL solve_kepuniv(dt,r0,sig0,mu,alpha,psi,s0,s1,s2,s3)
    f=1.d0-mu*s2/r0
    g=dt-mu*s3
    r=r0*s0+sig0*s1+mu*s2
    fdot=-mu*s1/(r0*r)
    gdot=1.d0-mu*s2/r
    x=x0*f+y0*g
    y=x0*fdot+y0*gdot
    RETURN
  END SUBROUTINE fser_propag

  SUBROUTINE fser_propag_der(x0,y0,t0,t,mu,x,y,dxydxy0)
    DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: x0,y0
    DOUBLE PRECISION, INTENT(IN) :: t,t0,mu
    DOUBLE PRECISION, DIMENSION(3), INTENT(OUT) :: x,y
    DOUBLE PRECISION, DIMENSION(6,6), INTENT(OUT) :: dxydxy0
    DOUBLE PRECISION :: s0,s1,s2,s3
    DOUBLE PRECISION r0, v0, vsize, prscal, alpha, sig0
    DOUBLE PRECISION psi, r,f,g,fdot,gdot
    DOUBLE PRECISION eps,period, dt
    DOUBLE PRECISION :: s4,s5,u_good
    DOUBLE PRECISION, DIMENSION(3,3) :: unit3,dxdx0,dxdy0,dydx0,dydy0
    DOUBLE PRECISION, DIMENSION(3,1) :: xdd, xd
    DOUBLE PRECISION, DIMENSION(3,2) :: x_xd, tmp 
    DOUBLE PRECISION, DIMENSION(1,3) :: xdd0, xd0
    DOUBLE PRECISION, DIMENSION(2,3) :: x0_xd0
    DOUBLE PRECISION, DIMENSION(2,2) :: ff
    r0=vsize(x0)
    sig0=prscal(x0,y0)
    v0=vsize(y0)
    alpha=v0**2-2*mu/r0 ! 2*E
    eps=100*epsilon(1.d0)
!    IF(alpha.lt.-eps)THEN
!      period=dpig*mu/(-alpha)**(1.5d0)
!      write(*,*)' period ', period
!      dt=t-t0-anint((t-t0)/period)*period
!    ELSE
      dt=t-t0
!    ENDIF
    CALL solve_kepuniv(dt,r0,sig0,mu,alpha,psi,s0,s1,s2,s3)
    f=1.d0-mu*s2/r0
    g=dt-mu*s3
    r=r0*s0+sig0*s1+mu*s2
    fdot=-mu*s1/(r0*r)
    gdot=1.d0-mu*s2/r
    x=x0*f+y0*g
    y=x0*fdot+y0*gdot
! additional functions for derivatives
    s4=(s2-psi**2/2)/alpha
    s5=(s3-psi**3/6)/alpha
    u_good=s2*dt+mu*(psi*s4-3*s5)
    CALL eye(3,unit3)
    xd(1:3,1)=y
    xd0(1,1:3)=y0
    xdd(1:3,1)=-x*mu/r**3
    xdd0(1,1:3)=-x0*mu/r0**3
    x_xd(1:3,1)=x
    x_xd(1:3,2)=y
    x0_xd0(1,1:3)=x0
    x0_xd0(2,1:3)=y0
    ff(1,1)= -(fdot*s1+(f-1.d0)/r0)/r0 
    ff(1,2)= -fdot*s2
    ff(2,1)= ((f-1.d0)*s1/r0)
    ff(2,2)= (f-1.d0)*s2 
    tmp= MATMUL(x_xd, ff)
    dxdx0=unit3*f + MATMUL(xd,xdd0)*u_good + MATMUL(tmp,x0_xd0)
    ff(1,1)= -fdot*s2 
    ff(1,2)= -(gdot-1.d0)*s2 
    ff(2,1)= (f-1.d0)*s2
    ff(2,2)= g*s2 
    tmp= MATMUL(x_xd, ff)
    dxdy0=unit3*g - MATMUL(xd,xd0)*u_good +  MATMUL(tmp,x0_xd0)
    ff(1,1)= -fdot*(s0/(r*r0)+ 1/r**2 +1/r0**2)
    ff(1,2)= -(fdot*s1+(gdot-1.d0)/r)/r
    ff(2,1)= (fdot*s1+(f-1.d0)/r0)/r0  
    ff(2,2)= fdot*s2
    tmp= MATMUL(x_xd, ff)
    dydx0=unit3*fdot + MATMUL(xdd,xdd0)*u_good + MATMUL(tmp,x0_xd0)
    ff(1,1)=-(fdot*s1+(gdot-1.d0)/r)/r
    ff(1,2)= -(gdot-1.d0)*s1/r
    ff(2,1)=fdot*s2 
    ff(2,2)=(gdot-1.d0)*s2
    tmp= MATMUL(x_xd, ff)
    dydy0=unit3*gdot - MATMUL(xdd,xd0)*u_good + MATMUL(tmp,x0_xd0)
    dxydxy0(1:3,1:3)=dxdx0
    dxydxy0(1:3,4:6)=dxdy0
    dxydxy0(4:6,1:3)=dydx0
    dxydxy0(4:6,4:6)=dydy0
    RETURN
  END SUBROUTINE fser_propag_der


  SUBROUTINE solve_kepuniv(dt,r0,sig0,mu,alpha,psi,s0,s1,s2,s3,conv_contr) 
    DOUBLE PRECISION, INTENT(IN):: dt,r0,sig0,mu,alpha
    DOUBLE PRECISION, INTENT(IN), OPTIONAL :: conv_contr ! covergence control
    DOUBLE PRECISION, INTENT(OUT):: psi, s0, s1, s2, s3
    DOUBLE PRECISION contr,dpsi,psi1,fun, funp
    DOUBLE PRECISION :: a0, u0, e0, cosu, enne, du
    INTEGER j
    INTEGER, PARAMETER :: jmax=100, itx=10
! control
    IF(PRESENT(conv_contr))THEN
       contr=conv_contr
    ELSE
       contr=100*epsilon(1.d0)
    ENDIF
    IF(alpha.lt.0.d0)THEN
       psi=dt/r0 !initial guess: only the sign matters
! +dt**2*sig0/(2*r0**3) !initial guess
       IF(abs(sig0).lt.contr)THEN
          a0=-mu/alpha
          e0=1.d0-r0/a0
          enne=sqrt(-alpha**3)/mu
          u0=pig
          if(enne*dt.lt.0.d0)u0=-pig
          DO j=1,itx
            du=-(u0-e0*sin(u0)-enne*dt)/(1.d0-e0*cos(u0))
            u0=u0+du
            IF(abs(du).lt.contr)EXIT
          ENDDO
          psi=u0/sqrt(-alpha)
       ENDIF
    ELSE
       psi=dt/r0 !initial guess: only the sign matters
!       psi=psi+dt**2*sig0/(2*r0**3) !improved initial guess
    ENDIF
    DO j=1,jmax
      CALL s_funct(psi,alpha,s0,s1,s2,s3)
      fun=r0*s1+sig0*s2+mu*s3-dt
      funp=r0*s0+sig0*s1+mu*s2
      dpsi=-fun/funp
      psi1=psi+dpsi
      IF(psi1*psi.lt.0.d0)THEN
         psi=psi/2.d0
      ELSE
         psi=psi1
      ENDIF
      IF(abs(dpsi).lt.contr.or.abs(dpsi).lt.contr*10*abs(psi)) RETURN
    ENDDO
    IF(abs(dpsi).gt.10*contr.or.abs(dpsi).gt.contr*100*abs(psi))THEN
       WRITE(*,*)' solve_kepuniv: poor convergence ', jmax,dpsi,psi,contr
    ENDIF
  END SUBROUTINE solve_kepuniv

  SUBROUTINE solve_peri(r0,sig0,peri,mu,alpha,psi,dt,conv_contr)
    DOUBLE PRECISION, INTENT(IN):: r0,sig0,mu,alpha,peri
    DOUBLE PRECISION, INTENT(IN), OPTIONAL :: conv_contr ! covergence control
    DOUBLE PRECISION, INTENT(OUT):: psi, dt
    DOUBLE PRECISION :: s0, s1, s2, s3
    DOUBLE PRECISION contr,dpsi, rp, rpp, psi0, psiper, r, a0, u0, e0, cosu
    INTEGER j, icount
    INTEGER, PARAMETER ::jmax=60
! control
    IF(PRESENT(conv_contr))THEN
       contr=conv_contr
    ELSE
       contr=100*epsilon(1.d0)
    ENDIF
    IF(alpha.le.0.d0)THEN
       a0=-mu/alpha
       e0=1.d0-peri/a0
       cosu=(1.d0-r0/a0)/e0
       IF(cosu.le.-1.d0)THEN
          u0=pig
       ELSEIF(cosu.ge.1.d0)THEN
          u0=0.d0
       ELSE
          u0=acos(cosu)
       ENDIF
       psi=-sign(1.d0,sig0)*u0/sqrt(-alpha)
    ELSE
       IF(abs(sig0).gt.contr)THEN
          psi=(peri-r0)/sig0
       ELSE
          psi=0.d0
       ENDIF
    ENDIF
!    WRITE(*,*)' psi0, psiper ', psi0, psiper
    icount=1  
!    psi=psi0
    DO j=1,jmax
      CALL s_funct(psi,alpha,s0,s1,s2,s3)
      r=r0*s0+sig0*s1+mu*s2
      rp=(r0*alpha+mu)*s1+sig0*s0
      rpp=(r0*alpha+mu)*s0+sig0*alpha*s1
      IF(ABS(rp).lt.contr)THEN
!         WRITE(*,*)j, psi,r,rp,rpp,' no correction'
         GOTO 1
      ENDIF
      dpsi=-rp/rpp
!      WRITE(*,*)j, psi,r,rp,rpp,dpsi
      psi=psi+dpsi
      IF(ABS(dpsi).lt.contr.or.abs(rp).lt.contr) GOTO 1
!      IF(ABS(psi).gt.psiper.or.psi*psi0.lt.0.d0)THEN
!         WRITE(*,*)' solve_peri: unstable newton ', icount, j, psi
!         icount=icount+1
!         psi=psi0/icount
!      ENDIF
    ENDDO
    WRITE(*,*)' solve_peri: poor convergence ', jmax,dpsi,contr
1   CONTINUE
    dt=r0*S1+sig0*s2+mu*s3
!    WRITE(*,*)' dt ', dt
    RETURN
  END SUBROUTINE solve_peri

  DOUBLE PRECISION FUNCTION r_of_psi(psi,alpha,r0,sig0,mu)
  DOUBLE PRECISION, INTENT(IN) :: psi, alpha, r0,sig0,mu
  DOUBLE PRECISION :: s0,s1,s2,s3
  CALL s_funct(psi,alpha,s0,s1,s2,s3)
  r_of_psi=r0*s0+sig0*s1+mu*s2
  RETURN
  END FUNCTION r_of_psi

  SUBROUTINE s_funct(psi,alpha,s0,s1,s2,s3,conv_contr)
    DOUBLE PRECISION, INTENT(IN) :: psi, alpha ! univ. ecc. anom., 2*energ 
    DOUBLE PRECISION, INTENT(OUT) :: s0,s1,s2,s3 ! Stumpff functions
    DOUBLE PRECISION, INTENT(IN), OPTIONAL :: conv_contr ! covergence control
    DOUBLE PRECISION beta, term2, term3, term0, term1, psi2, s02, s12
    DOUBLE PRECISION contr, dif1, dif2, sval
    INTEGER j,jj,nhalf
    INTEGER, PARAMETER:: jmax=70, halfmax=8
    DOUBLE PRECISION, PARAMETER :: betacontr=1.d2 !to make it not operational
! control
    IF(PRESENT(conv_contr))THEN
       contr=conv_contr
    ELSE
       contr=100*epsilon(1.d0)
    ENDIF
! compute s2,s3
    beta=alpha*psi**2 
!    WRITE(*,*)psi,beta,alpha
! if beta>1000 then...
    IF(abs(beta).lt.betacontr)THEN
       term2=psi**2/2
       term3=term2*psi/3
       s2=term2
       s3=term3
       DO j=1,jmax
          term2=term2*beta/((2*j+1)*(2*j+2))
          s2=s2+term2
          IF(abs(term2).lt.contr)EXIT
       ENDDO
       IF(j.gt.jmax)WRITE(*,*)' j, term2, beta ', j-1, term2,beta
       DO j=1,jmax
          term3=term3*beta/((2*j+2)*(2*j+3))
          s3=s3+term3
          IF(abs(term3).lt.contr)EXIT
       ENDDO
       IF(j.gt.jmax)WRITE(*,*)' j, term3,beta ', j-1, term3, beta
! recursive formulae
       s1=psi+alpha*s3
       s0=1.d0+alpha*s2
!       WRITE(*,*)s0,s1,s2,s3
    ELSE
       psi2=psi
       nhalf=0
       DO jj=1,halfmax
         psi2=psi2*0.5d0
         nhalf=nhalf+1
         beta=alpha*psi2**2
         IF(abs(beta).lt.betacontr)EXIT
       ENDDO
       WRITE(*,*)' halfing', psi,nhalf,beta
       term0=1.d0
       term1=psi2
       s0=1.d0
       s1=psi2
       DO j=1,jmax
          term0=term0*beta/((2*j-1)*(2*j))
          s0=s0+term0
          IF(abs(term0).lt.contr)EXIT
       ENDDO
       IF(j.gt.jmax)WRITE(*,*)' j, term0, beta ', j-1, term0,beta
       DO j=1,jmax
          term1=term1*beta/((2*j)*(2*j+1))
          s1=s1+term1
          IF(abs(term1).lt.contr)EXIT
       ENDDO
       IF(j.gt.jmax)WRITE(*,*)' j, term1, beta ', j-1, term1,beta
! duplication formula
       DO jj=1,nhalf
         s02=2*s0**2-1.d0
         s12=2*s0*s1
         s0=s02
         s1=s12
       ENDDO
! recursive formulae
       s3=(s1-psi)/alpha
       s2=(s0-1.d0)/alpha
!       WRITE(*,*)s0,s1,s2,s3
    ENDIF
! controls
    dif1=2*s2+alpha*s2**2-s1**2
    dif2=s2+s0*s2-s1**2
    sval=s1**2+s2**2
    IF(abs(dif1).gt.100*contr*sval.or.abs(dif2).gt.100*contr*sval)THEN
    IF(abs(dif1).gt.1000*contr.or.abs(dif2).gt.1000*contr)THEN
       WRITE(*,100)psi,alpha,s0,s1,s2,s3
  100  FORMAT('psi, alpha, s functions ', 1p, 6d13.5)
       WRITE(*,101)dif1,dif2,contr*sval
  101  FORMAT('check id1, check id2, contr ', 1p, 3d13.5)
    ENDIF
    ENDIF
  END SUBROUTINE s_funct
END MODULE ever_pitkin

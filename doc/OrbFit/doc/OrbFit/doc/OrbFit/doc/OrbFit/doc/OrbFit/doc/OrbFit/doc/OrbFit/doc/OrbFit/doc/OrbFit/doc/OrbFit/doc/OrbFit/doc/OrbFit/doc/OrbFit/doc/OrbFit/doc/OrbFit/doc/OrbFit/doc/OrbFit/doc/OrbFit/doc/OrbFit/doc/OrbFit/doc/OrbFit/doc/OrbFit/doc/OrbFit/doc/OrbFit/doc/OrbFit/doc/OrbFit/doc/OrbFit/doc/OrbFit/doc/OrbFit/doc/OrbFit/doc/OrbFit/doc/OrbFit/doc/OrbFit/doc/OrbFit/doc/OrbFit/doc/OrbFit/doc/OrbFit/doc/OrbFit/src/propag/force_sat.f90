MODULE force_sat
  USE fund_const
  USE planet_masses, ONLY: gmearth,gmoon,dmea
  IMPLICIT NONE 
  PRIVATE

! public routines
  PUBLIC :: eamoon_mass, forcesat
! shared data
  

CONTAINS

! ===========================================================           
! FORCESAT : accelerations acting on a satellite of Earth                   
! ===========================================================           
! version 3.6.1; 29 October 2008     
  SUBROUTINE forcesat(x,v,t0,f,nd,idc,xxpla,ips,imem, derf0)
    USE perturbations
    USE spher_harm
    USE iers_ser
! ======INPUT===================                                        
! dimension of position vector                                          
    INTEGER, INTENT(IN) :: nd 
! Position, velocity,  time                                             
    REAL(KIND=dkind), INTENT(IN) :: x(nd), v(nd) ,t0 
! flag for recomputation, memory location                               
    INTEGER, INTENT(IN) :: ips,imem 
! WARNING: for now storage/partial recomputation not implemented
! ======OUTPUT===================                                       
! acceleration                                                          
    REAL(KIND=dkind), INTENT(OUT) :: f(nd) 
! Positions and vel of the planet involved in  close-app                
! stored only if idc.ne.0                                               
    INTEGER, INTENT(OUT) :: idc 
    REAL(KIND=dkind), INTENT(OUT) :: xxpla(6) 
    REAL(KIND=dkind), INTENT(OUT), OPTIONAL :: derf0(3,3)
! ======END INTERFACE============== 
! to store data; not used yet
!    INTEGER, PARAMETER:: memx=30 ! stored planets/rotations array 
! for spherical harmonics 
    REAL(KIND=dkind), DIMENSION(3) :: y ! body-fixed position
    INTEGER :: lmi, lmax ! index in array of harmonics, max value

    REAL(KIND=dkind) rot(3,3),rot1(3,3),rot2(3,3) ! rotation matrices to BF
    REAL(KIND=dkind) rott(3,3) ! ! rotation matrix from BF
    REAL(KIND=dkind) derf(3,3),dfb(3,3),dfc(3,3) ! jacobian, d2B/dt2, d2C/dt2
    LOGICAL partials ! need for partial derivatives
! temporaries to store harmonics                                       
    REAL(KIND=dkind) :: army(3,4,narmx),arm(3,4,narmx), tmp3(3,3), plmi
    REAL(KIND=dkind) :: armacc(3,4), de, e, g, armaccy(3,4)
! conversion of times 
    REAL(KIND=dkind) :: sec1
    INTEGER mjd1
! two body part
    REAL(KIND=dkind) :: r2,r3 !, de
! for sun-moon perturbations
    REAL(KIND=dkind), DIMENSION(6) :: xsun,xmoon
    REAL(KIND=dkind) :: ef(3,2),pk2cur
    LOGICAL :: equ
! perturbative acceleration (and partials): planetary and tidal effects 
    REAL(KIND=dkind), DIMENSION(3,4) :: dplan 
! partial derivatives of the acceleration w.r.t. to the Love number $k_2$  
    REAL(KIND=dkind), DIMENSION(3) :: dlove
! non gravitational perturbations
    REAL(KIND=dkind), DIMENSION(3) :: accrad, accradsec
    REAL(KIND=dkind) :: vsize
! loop indexes
    INTEGER i,j,k
! initialization
    idc=0  ! close approach flag
    xxpla=0.d0
! are partial required?
    IF(nd.eq.3)THEN
       partials=.false.
    ELSEIF(nd.gt.3)THEN
       partials=.true.
    ELSE
       STOP '***** nd strange in force_sat ******'
    ENDIF
! two-body term
    r2=x(1)**2+x(2)**2+x(3)**2 
    r3=r2*dsqrt(r2) 
    de=-gmearth/r3
    f(1:3)=de*x(1:3) 
    IF(partials)THEN
!  partials of the two body acceleration                                
       e=-3.d0*de/r2 
       DO i=1,3
          DO j=1,i 
             g=e*x(i)*x(j) 
             IF(i.eq.j) g=g+de 
             derf(i,j)=g
             derf(j,i)=g
          ENDDO
       ENDDO
!       f(13)=1.d0
!       f(17)=1.d0
!       f(21)=1.d0            
    ENDIF
! luni-solar perturbations
    IF(ipla.eq.2)THEN
       equ=.TRUE.
       CALL sunmoon_car(t0,equ,xsun,xmoon) 
       ef(1:3,1)=xsun(1:3)
       ef(1:3,2)=xmoon(1:3)
! love number for Earth
       pk2cur=0.3
       CALL sunmoon_pert(x,ef,dplan,dlove,pk2cur,partials)
       f(1:3)=f(1:3)+dplan(1:3,1)
       IF(partials)THEN
          derf(1:3,1:3)=derf(1:3,1:3)+ dplan(1:3,2:4)  
       ENDIF
    ELSEIF(ipla.gt.0)THEN
       WRITE(*,*)' force_sat: unknown option ipla=', ipla
       STOP
    ENDIF
! non gravitational perturbations
    IF(irad.gt.0)THEN
       IF(ipla.eq.0)THEN
          equ=.TRUE.
          CALL sunmoon_car(t0,equ,xsun,xmoon) 
       ENDIF 
      IF(irad.eq.1.or.irad.eq.3)THEN
         CALL radp(x(1:3),xsun(1:3),accrad) 
         f(1:3)=f(1:3)+accrad*amrat
      ENDIF
      IF(irad.eq.2.or.irad.eq.3)THEN
         CALL secacc(x(1:3),v(1:3),xsun(1:3),accradsec) 
         f(1:3)=f(1:3)+accradsec*amratsec
      ENDIF
! WARNING: no partial derivatives from these (possible problem in penumbra)
! no contributions to variational equation
    ENDIF
! spherical harmonics 
! rotation matrix (and time derivatives) to transform                 
! from the inertial system of the JPL ephemerides                     
! to the body fixed system used in the integration 
    IF(ites.ge.2)THEN
       mjd1=t0
       sec1=(t0-mjd1)*86400.d0                   
       CALL rotsys('MEAN',mj2000,s2000,'BF  ',mjd1,sec1,rot,rot1,rot2,1) 
       rott=TRANSPOSE(rot)
       y=matmul(rot,x(1:3))
! spherical harmonics
       CALL parm10(y,ites,army,partials)
       lmax=numcoe(ites)
! adding up spherical harmonics
       armaccy=0.d0
       DO lmi=5,lmax
          plmi=harmco(lmi) 
          armaccy(1:3,1)=armaccy(1:3,1)+plmi*army(1:3,1,lmi)
          IF(partials)THEN
             armaccy(1:3,2:4)=armaccy(1:3,2:4)+plmi*army(1:3,2:4,lmi)
          ENDIF
       ENDDO
 ! rotation of acceleration
       CALL prodmv(armacc(1:3,1),rott,armaccy(1:3,1))
       f(1:3)=f(1:3)+armacc(1:3,1)
       IF(partials)THEN
           tmp3=matmul(armaccy(1:3,2:4),rot) 
           derf(1:3,1:3)=derf(1:3,1:3)+ matmul(rott,tmp3)
       ENDIF
!       WRITE(25,100)t0,y(1:3),armaccy(1:3,1)
!       WRITE(26,100)t0,x(1:3),armacc(1:3,1)
!100    FORMAT(F12.6,1P,3(1X,D15.8),3(1X,D15.8))
! rotation of acceleration
!       armacc=0.d0
!       DO lmi=5,lmax 
! harmonic coefficient
!          plmi=harmco(lmi) 
!          arm(1:3,1,lmi)=matmul(rott,army(1:3,1,lmi)) 
!          armacc(1:3,1)=armacc(1:3,1)+plmi*arm(1:3,1,lmi)
!          f(1:3)=f(1:3)+armacc(1:3,1)
!          IF(partials)THEN 
! rotation of second derivatives
!             tmp3=matmul(army(1:3,2:4,lmi),rot) 
!             arm(1:3,2:4,lmi)=matmul(rott,tmp3) 
! update of jacobian
!             derf(1:3,1:3)=derf(1:3,1:3)+ plmi* arm(1:3,2:4,lmi)
!          ENDIF
!       ENDDO

    ENDIF

! right hand side of variational equation
    IF(partials)THEN
       IF(PRESENT(derf0))THEN
! output matrix of partials
          derf0=derf
       ENDIF
       DO k=1,3 ! k=column index
          CALL prodmv(dfb(1:3,k), derf,x(3*k+1:3*k+3))
! right hand side of equation for B=dx/dx0
          f(k*3+1:k*3+3)=dfb(1:3,k)
          CALL prodmv(dfc(1:3,k), derf,x(3*k+10:3*k+12))
! right hand side of equation for C=dv/dx0
          f(k*3+10:k*3+12)=dfc(1:3,k)
       ENDDO
    ENDIF
  END SUBROUTINE forcesat

    SUBROUTINE eamoon_mass
! ======== JPL EPHEM ============================                       
! standard JPL header                                                   
      INCLUDE 'jplhdr.h90' 
! ===============================================
! mass of Earth needs to be known anyway
      gmearth= cval(11)/(1.d0+1.d0/cval(8))
! moon mass in solar masses 
!      gmoon=cval(11)/((1.d0+emrat)*gmsun)
      gmoon=gmearth/cval(8)  
! If you want GM in solar masses you need to divide by gmsun
!     gmearth=gmearth/gms
!     gmoon=gmoon/gms
      dmea=reau*1.02d0
    END SUBROUTINE eamoon_mass



END MODULE force_sat

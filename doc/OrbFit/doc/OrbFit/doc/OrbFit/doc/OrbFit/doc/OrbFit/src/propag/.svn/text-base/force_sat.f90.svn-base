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
  SUBROUTINE forcesat(x,v,t0,f,nd,idc,xxpla,ips,imem)
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
! ======END INTERFACE============== 
! to store data; not used yet
!    INTEGER, PARAMETER:: memx=30 ! stored planets/rotations array 
! for spherical harmonics 
    REAL(KIND=dkind), DIMENSION(3) :: y ! body-fixed position
    INTEGER :: lmi, lmax ! index in array of harmonics, max value

    REAL(KIND=dkind) rot(3,3),rot1(3,3),rot2(3,3) ! rotation matrices to BF
    REAL(KIND=dkind) rott(3,3) ! ! rotation matrix from BF
    LOGICAL partials ! need for partial derivatives
! temporaries to store harmonics                                       
    REAL(KIND=dkind) :: army(3,4,narmx),arm(3,4,narmx), tmp3(3,3), plmi
    REAL(KIND=dkind) :: armacc(3,4), de, e, g
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
    REAL(KIND=dkind), DIMENSION(3) :: accrad
! loop indexes
    INTEGER i,j 
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
    de=gmearth/r3
    f(1:3)=-de*x(1:3) 
    IF(partials)THEN
!  partials of the two body acceleration                                
       e=-3.d0*de/r2 
       DO i=1,3
          DO j=1,i 
             g=e*x(i)*x(j) 
             IF(i.eq.j) g=g+de 
             f(3*i+j)=g
             f(3*j+i)=g
          ENDDO
       ENDDO
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
          DO j=1,3
             f(3*j+1:3*j+3)= f(3*j+1:3*j+3)+dplan(1:3,j+1)
          ENDDO
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
      CALL radp(x,xsun(1:3),accrad) 
      f(1:3)=f(1:3)+accrad
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
 ! rotation of acceleration
       armacc=0.d0
       DO lmi=5,lmax 
          arm(1:3,1,lmi)=matmul(rott,army(1:3,1,lmi)) 
          armacc(1:3,1)=armacc(1:3,1)+plmi*arm(1:3,1,lmi)
          IF(partials)THEN 
! rotation of second derivatives
             tmp3=matmul(army(1:3,2:4,lmi),rot) 
             arm(1:3,2:4,lmi)=matmul(rott,tmp3) 
! harmonic coefficient
             plmi=harmco(lmi) 
! update of jacobian
             DO j=1,3
                f(3*j+1:3*j+3)=f(3*j+1:3*j+3)+plmi*arm(1:3,j+1,lmi)
             ENDDO
          ENDIF
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

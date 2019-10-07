! ******************************************************************
! ********** COMPUTATION of USEFUL ORBIT COEFFICIENTS  *************
! ******************************************************************
! ********** written by GIOVANNI F. GRONCHI (Oct. 2004) ************
! ==================================================================

MODULE critical_points_shift
  IMPLICIT NONE
  PRIVATE
! elements in the first 2 rows of the cracovians (see Sitarski 1968)
  REAL(KIND=8) :: PPx,PPy,PPz,QQx,QQy,QQz
  REAL(KIND=8) :: px,py,pz,qx,qy,qz
  REAL(KIND=8) :: A1x,A1y,A1z,B1x,B1y,B1z,a2x,a2y,a2z,b2x,b2y,b2z 
!                
  REAL(KIND=8) :: KK,LL,MM,NN,Ppl,pcom,Epl,ecom
  REAL(KIND=8) :: a2Q1,b2Q1,B1a2,B1b2,A1a2,A1b2
  REAL(KIND=8) :: A1q2,B1q2,b2A1,b2B1,a2A1,a2B1
! ---------------- Sylvester matrix elements ----------
  INTEGER, PARAMETER :: N = 16 ! number of evaluations of the 15th deg poly
  INTEGER, PARAMETER :: expo = 4 ! 2^expo = N
 INTEGER, PARAMETER :: poldeg = 16 ! polynomial degree
  REAL(KIND=8),DIMENSION(N) :: p0,p1,p2
  REAL(KIND=8),DIMENSION(N) :: q0,q1,q2,q3,q4
  REAL(KIND=8),DIMENSION(N) :: r31,r32,r33,r34,r35
  REAL(KIND=8),DIMENSION(N) :: r36,r41,r45,r46
! to write all the complex roots of 1-D poly r(w)
  LOGICAL :: follow_roots ! follow_roots = .false. do not write
                  ! follow_roots = .true. write the roots of r(w)
!  SAVE
! -----------------------------------------------------
  PUBLIC :: PPx,PPy,PPz,QQx,QQy,QQz,px,py,pz,qx,qy,qz
  PUBLIC :: A1x,A1y,A1z,B1x,B1y,B1z,a2x,a2y,a2z,b2x,b2y,b2z
  PUBLIC :: KK,LL,MM,NN,Ppl,pcom,Epl,ecom
  PUBLIC :: a2Q1,b2Q1,B1a2,B1b2,A1a2,A1b2,A1q2,B1q2,b2A1,b2B1,a2A1,a2B1
  PUBLIC :: p0,p1,p2,q0,q1,q2,q3,q4
  PUBLIC :: N,expo,poldeg
  PUBLIC :: r31,r32,r33,r34,r35,r36,r41,r45,r46
  PUBLIC :: follow_roots
! subroutines
  PUBLIC :: orbitcoe_shift,matrixdat_ta_shift,compmodsylv16_shift
  PUBLIC :: solvesystem_ta_shift,hess_ta,int_eval_ta
  PUBLIC :: d2eval_ta
CONTAINS

! =====================================================
  SUBROUTINE orbitcoe_shift(alpha,beta,elpl,elcom) 
! HINT: angles must be passed in radians 
    IMPLICIT NONE 
! shift angles
    REAL(KIND=8),INTENT(IN) :: alpha,beta
! cometary elements of the planet and of the comet
    REAL(KIND=8),INTENT(IN),DIMENSION(5) :: elpl
    REAL(KIND=8),INTENT(IN),DIMENSION(5) :: elcom 
! ---------------- end interface ---------------------
    REAL(KIND=8) :: Qpl,i1,Omnod1,omega1 ! planet cometary elements
    REAL(KIND=8) :: q,i2,Omnod2,omega2 ! comet cometary elements
    REAL(KIND=8) :: deltaOm ! diff of node longitudes
! elements in the first 2 rows of the cracovians (see Sitarski 1968)
    REAL(KIND=8) :: sa,ca,sb,cb
! ==================================================================
 
    Qpl    = elpl(1)
    Epl    = elpl(2)
    i1     = elpl(3)
    Omnod1 = elpl(4)
    omega1 = elpl(5)
    q      = elcom(1)
    ecom   = elcom(2)
    i2     = elcom(3)
    Omnod2 = elcom(4)
    omega2 = elcom(5)
    deltaOm = Omnod2-Omnod1

    PPx = cos(omega1) 
    PPy = sin(omega1)*cos(i1)
    PPz = sin(omega1)*sin(i1)
!    write(*,*)'PPx',PPx,'PPy',PPy,'PPz',PPz
    QQx = -sin(omega1) 
    QQy = cos(omega1)*cos(i1)
    QQz = cos(omega1)*sin(i1)
!    write(*,*)'QQx',QQx,'QQy',QQy,'QQz',QQz
    px = cos(omega2)*cos(deltaOm) - &
         & sin(omega2)*cos(i2)*sin(deltaOm)
    py = cos(omega2)*sin(deltaOm) + &
         & sin(omega2)*cos(i2)*cos(deltaOm)
    pz = sin(omega2)*sin(i2)
!    write(*,*)'px',px,'py',py,'pz',pz
    qx = -sin(omega2)*cos(deltaOm) - &
         & cos(omega2)*cos(i2)*sin(deltaOm)
    qy = -sin(omega2)*sin(deltaOm) + &
         & cos(omega2)*cos(i2)*cos(deltaOm)
    qz = cos(omega2)*sin(i2)
!    write(*,*)'qx',qx,'qy',qy,'qz',qz

    KK = PPx*px + PPy*py + PPz*pz
    LL = QQx*px + QQy*py + QQz*pz
    MM = PPx*qx + PPy*qy + PPz*qz
    NN = QQx*qx + QQy*qy + QQz*qz
    Ppl = Qpl*(1.d0+Epl)
    pcom = q*(1.d0+ecom)
!    write(*,*)'KK=',KK,'LL=',LL,'MM=',MM,'NN=',NN
!    write(*,*)'Ppl=',Ppl,'pcom=',pcom,'Epl=',Epl,'ecom=',ecom

! -------- SHIFT COEFFICIENTS ---------
    ca = cos(alpha)
    sa = sin(alpha)
    cb = cos(beta)
    sb = sin(beta)
    A1x = ca*PPx + sa*QQx
    A1y = ca*PPy + sa*QQy
    A1z = ca*PPz + sa*QQz
!    write(*,*)'Ax',A1x,'Ay',A1y,'Az',A1z
    B1x = -sa*PPx + ca*QQx
    B1y = -sa*PPy + ca*QQy
    B1z = -sa*PPz + ca*QQz
!    write(*,*)'Bx',B1x,'By',B1y,'Bz',B1z
    a2x = cb*px + sb*qx
    a2y = cb*py + sb*qy
    a2z = cb*pz + sb*qz
!    write(*,*)'ax',a2x,'ay',a2y,'az',a2z
    b2x = -sb*px + cb*qx
    b2y = -sb*py + cb*qy
    b2z = -sb*pz + cb*qz
!    write(*,*)'bx',b2x,'by',b2y,'bz',b2z
!    A1Q1 = A1x*QQx + A1y*QQy + A1z*QQz
!    B1Q1 = B1x*QQx + B1y*QQy + B1z*QQz
    a2Q1 = a2x*QQx + a2y*QQy + a2z*QQz
    b2Q1 = b2x*QQx + b2y*QQy + b2z*QQz
    B1a2 = B1x*a2x + B1y*a2y + B1z*a2z
    B1b2 = B1x*b2x + B1y*b2y + B1z*b2z
    A1a2 = A1x*a2x + A1y*a2y + A1z*a2z 
    A1b2 = A1x*b2x + A1y*b2y + A1z*b2z 
!    a2q2 = a2x*qx + a2y*qy + a2z*qz
!    b2q2 = b2x*qx + b2y*qy + b2z*qz
    A1q2 = A1x*qx + A1y*qy + A1z*qz
    B1q2 = B1x*qx + B1y*qy + B1z*qz
    b2A1 = A1b2
    b2B1 = B1b2
    a2A1 = A1a2
    a2B1 = B1a2
  END SUBROUTINE orbitcoe_shift
  
! ******************************************************************
! ***** Given the orbcoe's, get the coeff's of the polynomials *****
! ***** of the modified Sylvester matrix ***************************
! *********** written by GIOVANNI F. GRONCHI (Oct.2004) ************
! ==================================================================
  SUBROUTINE matrixdat_ta_shift(alpha,beta)
    IMPLICIT NONE 
    REAL(KIND=8),INTENT(IN) :: alpha,beta ! shift angles                                             
! ========== end interface =========================================
    REAL(KIND=8) :: AE,BE,CE,DE,EE,FE ! coeff of the linear combination 
                                      ! of the rows of Sylvester's matrix
    INTEGER :: i,j,l ! loop indexes   
    REAL(KIND=8) :: sa,ca,sb,cb ! auxiliary
! ==================================================================
 
    sa=sin(alpha)
    ca=cos(alpha)
    sb=sin(beta)
    cb=cos(beta)

    AE = Epl*sa/(1.d0+Epl*ca)
    BE = Epl*sa/(1.d0-Epl*ca)
    CE = (Epl**2-1.d0+Epl**2*sa**2)/((1.d0+Epl*ca)**2)
    DE = (Epl**2-1.d0+Epl**2*sa**2)/((1.d0-Epl*ca)**2)
    EE = (Epl*sa/((1.d0+Epl*ca)**3))*(3.d0*(Epl**2-1.d0) + Epl**2*sa**2)
    FE = (Epl*sa/((1.d0-Epl*ca)**3))*(3.d0*(Epl**2-1.d0) + Epl**2*sa**2)

! p0(t) = p0(5)*t^4 + p0(4)*t^3 + p0(3)*t^2 + p0(2)*t + p0(1)
    p0(5) = -Ppl*A1b2+ecom*pcom*sb-Ppl*ecom**2*cb*A1q2+ &
         &ecom*pcom*sb*Epl*ca+Ppl*ecom*cb*A1b2+Ppl*ecom*A1q2
    p0(4) = -2*Ppl*A1a2+2*Ppl*ecom*cb*A1a2-2*Ppl*ecom**2*sb*A1q2- &
         & 2*ecom*pcom*cb*Epl*ca-2*ecom*pcom*cb+2*Ppl*ecom*sb*A1b2
    p0(3) = 2*Ppl*ecom*(-cb*A1b2+2*sb*A1a2+A1q2)
    p0(2) = -2*Ppl*A1a2-2*Ppl*ecom*cb*A1a2-2*Ppl*ecom**2*sb*A1q2- &
         & 2*ecom*pcom*cb-2*ecom*pcom*cb*Epl*ca-2*Ppl*ecom*sb*A1b2
    p0(1) = -ecom*pcom*sb+Ppl*A1b2+Ppl*ecom*A1q2+Ppl*ecom**2*cb*A1q2+ &
         & Ppl*ecom*cb*A1b2-ecom*pcom*sb*Epl*ca

!    write(*,*)'pcom, Ppl',pcom, Ppl
!    write(*,*)'p0',p0(1:5)

! p1(t) = p1(5)*t^4 + p1(4)*t^3 + p1(3)*t^2 + p1(2)*t + p1(1)
    p1(5) = -2*ecom*pcom*sb*Epl*sa-2*Ppl*ecom**2*cb*B1q2+ &
         & 2*Ppl*ecom*cb*B1b2+2*Ppl*ecom*B1q2-2*Ppl*B1b2
    p1(4) = 4*ecom*pcom*cb*Epl*sa- &
         & 4*Ppl*ecom**2*sb*B1q2+4*Ppl*ecom*sb*B1b2+4*Ppl*ecom*cb*a2B1-4*Ppl*a2B1
    p1(3) = 4*Ppl*ecom*(2*sb*a2B1-cb*B1b2+B1q2)
    p1(2) = -4*Ppl*a2B1+4*ecom*pcom*cb*Epl*sa- &
         & 4*Ppl*ecom**2*sb*B1q2-4*Ppl*ecom*sb*B1b2-4*Ppl*ecom*cb*a2B1
    p1(1) = 2*Ppl*B1b2+2*Ppl*ecom*B1q2+2*Ppl*ecom**2*cb*B1q2+ &
         & 2*Ppl*ecom*cb*B1b2+2*ecom*pcom*sb*Epl*sa

!    write(*,*)'p1',p1(1:5)

! p2(t) = p2(5)*t^4 + p2(4)*t^3 + p2(3)*t^2 + p2(2)*t + p2(1)
    p2(5) = Ppl*A1b2+ecom*pcom*sb+Ppl*ecom**2*cb*A1q2- &
         & ecom*pcom*sb*Epl*ca-Ppl*ecom*cb*A1b2-Ppl*ecom*A1q2
    p2(4) = 2*Ppl*A1a2-2*Ppl*ecom*cb*A1a2+2*Ppl*ecom**2*sb*A1q2 &
         & +2*ecom*pcom*cb*Epl*ca-2*ecom*pcom*cb-2*Ppl*ecom*sb*A1b2
    p2(3) = -2*Ppl*ecom*(-cb*A1b2+2*sb*A1a2+A1q2)
    p2(2) = 2*Ppl*A1a2+2*Ppl*ecom*cb*A1a2+2*Ppl*ecom**2*sb*A1q2 &
         & -2*ecom*pcom*cb+2*ecom*pcom*cb*Epl*ca+2*Ppl*ecom*sb*A1b2
    p2(1) = -ecom*pcom*sb-Ppl*A1b2-Ppl*ecom*A1q2-Ppl*ecom**2*cb*A1q2 &
         & -Ppl*ecom*cb*A1b2+ecom*pcom*sb*Epl*ca

!    write(*,*)'p2',p2(1:5)

! q0(t) = q0(3)*t^2 + q0(2)*t + q0(1)
    q0(3) = pcom*Epl*a2Q1+Epl*Ppl*sa+pcom*Epl*ca*B1a2+ &
         & pcom*Epl**2*ca*a2Q1-Epl*Ppl*sa*ecom*cb+pcom*B1a2
    q0(2) = -2*pcom*Epl*b2Q1-2*pcom*Epl*ca*B1b2-2*pcom*Epl**2*ca*b2Q1 &
         & -2*Epl*Ppl*sa*ecom*sb-2*pcom*B1b2
    q0(1) = Epl*Ppl*sa-pcom*B1a2-pcom*Epl*ca*B1a2-pcom*Epl**2*ca*a2Q1 &
           & +Epl*Ppl*sa*ecom*cb-pcom*Epl*a2Q1
            
!   write(*,*)'matrixdat: q0=',q0(1:3)

! q1(t) = q1(3)*t^2 + q1(2)*t + q1(1)
   q1(3) = 2*Epl*Ppl*ca-2*pcom*Epl**2*sa*a2Q1-2*Epl*Ppl*ca*ecom*cb- &
        &2*pcom*Epl*ca*A1a2-2*pcom*Epl*sa*B1a2-2*pcom*A1a2 
   q1(2) = 4*pcom*A1b2-4*Epl*Ppl*ca*ecom*sb+4*pcom*Epl**2*sa*b2Q1+ &
        & 4*pcom*Epl*sa*B1b2+4*pcom*Epl*ca*A1b2
   q1(1) = 2*Epl*Ppl*ca+2*pcom*Epl**2*sa*a2Q1+2*pcom*Epl*ca*A1a2+ &
        & 2*pcom*Epl*sa*B1a2+2*Epl*Ppl*ca*ecom*cb+2*pcom*A1a2
                     
!   write(*,*)'matrixdat: q1=',q1(1:3)

! q2(t) = q2(3)*t^2 + q2(2)*t + q2(1)
   q2(3) = 2*pcom*Epl*(a2Q1-ca*B1a2+2*sa*A1a2)      
   q2(2) = -4*pcom*Epl*(b2Q1+2*sa*A1b2-ca*B1b2)
   q2(1) = -2*pcom*Epl*(a2Q1-ca*B1a2+2*sa*A1a2)
                      
!   write(*,*)'q2',q2(1:3)

! q3(t) = q3(3)*t^2 + q3(2)*t + q3(1)
   q3(3) = -2*pcom*A1a2-2*Epl*Ppl*ca*ecom*cb+2*Epl*Ppl*ca+2*pcom*Epl*sa*B1a2+ &
        & 2*pcom*Epl*ca*A1a2-2*pcom*Epl**2*sa*a2Q1 
   q3(2) = 4*pcom*A1b2+4*pcom*Epl**2*sa*b2Q1-4*Epl*Ppl*ca*ecom*sb- &
        & 4*pcom*Epl*ca*A1b2-4*pcom*Epl*sa*B1b2
   q3(1) = 2*Epl*Ppl*ca+2*pcom*A1a2+2*Epl*Ppl*ca*ecom*cb-2*pcom*Epl*sa*B1a2 &
        & -2*pcom*Epl*ca*A1a2+2*pcom*Epl**2*sa*a2Q1

!   write(*,*)'q3',q3(1:3)
            
! q4(t) = q4(3)*t^2 + q4(2)*t + q4(1)
   q4(3) = pcom*Epl*a2Q1-pcom*B1a2-Epl*Ppl*sa+Epl*Ppl*sa*ecom*cb- &
        & pcom*Epl**2*ca*a2Q1+pcom*Epl*ca*B1a2
   q4(2) = 2*pcom*B1b2-2*pcom*Epl*b2Q1+2*pcom*Epl**2*ca*b2Q1+ &
        & 2*Epl*Ppl*sa*ecom*sb-2*pcom*Epl*ca*B1b2
   q4(1) = pcom*B1a2-pcom*Epl*a2Q1-Epl*Ppl*sa+pcom*Epl**2*ca*a2Q1 &
        & -Epl*Ppl*sa*ecom*cb-pcom*Epl*ca*B1a2

!   write(*,*)'q4',q4(1:3)
 
! ===============================================================
   
   r31(3) = -2*Ppl*(ecom*A1q2-Epl**2*ecom*A1q2+ecom*A1q2*Epl**2*ca**2- &
        & ecom*A1q2*Epl*ca-Epl**2*sa*ecom*B1q2*ca+Epl*sa*ecom*B1q2+A1b2*Epl*ca- &
        & Epl*sa*B1b2-A1b2*Epl**2*ca**2+Epl**2*sa*B1b2*ca-A1b2+Epl**2*A1b2)/(-1+Epl*ca)**2

   r31(2) = -4*Ppl*(-A1a2+Epl**2*sa*a2B1*ca-A1a2*Epl**2*ca**2- &
        & Epl*sa*a2B1+A1a2*Epl*ca+Epl**2*A1a2)/(-1+Epl*ca)**2

   r31(1) = -2*Ppl*(A1b2*Epl**2*ca**2+ecom*A1q2-Epl**2*A1b2-A1b2*Epl*ca- &
        & Epl**2*sa*ecom*B1q2*ca+ecom*A1q2*Epl**2*ca**2+Epl*sa*ecom*B1q2+Epl*sa*B1b2- &
        & Epl**2*sa*B1b2*ca-ecom*A1q2*Epl*ca+A1b2-Epl**2*ecom*A1q2)/(-1+Epl*ca)**2

   r32(3) = -2*Ppl*(-ecom*B1q2+ecom*B1q2*Epl**2*ca**2+Epl**2*sa*ecom*A1q2*ca- &
        & Epl**2*sa*A1b2*ca-B1b2*Epl**2*ca**2+B1b2)/(-1+Epl*ca)/(Epl*ca+1)

   r32(2) = 4*Ppl*(Epl**2*sa*A1a2*ca-a2B1+a2B1*Epl**2*ca**2)/(-1+Epl*ca)/(Epl*ca+1)
   
   r32(1) = -2*Ppl*(-ecom*B1q2+ecom*B1q2*Epl**2*ca**2+Epl**2*sa*ecom*A1q2*ca+ &
        & Epl**2*sa*A1b2*ca-B1b2+B1b2*Epl**2*ca**2)/(-1+Epl*ca)/(Epl*ca+1)
   
   r33(3) = 2*Ppl*(ecom*A1q2-Epl**2*ecom*A1q2+ecom*A1q2*Epl**2*ca**2+ &
        & ecom*A1q2*Epl*ca-Epl**2*sa*ecom*B1q2*ca-Epl*sa*ecom*B1q2-A1b2*Epl*ca+ &
        & Epl*sa*B1b2-A1b2*Epl**2*ca**2+Epl**2*sa*B1b2*ca-A1b2+Epl**2*A1b2)/(Epl*ca+1)**2

   r33(2) = 4*Ppl*(-A1a2+Epl**2*sa*a2B1*ca-A1a2*Epl**2*ca**2+Epl*sa*a2B1-A1a2*Epl*ca+ &
        & Epl**2*A1a2)/(Epl*ca+1)**2
   
   r33(1) = 2*Ppl*(A1b2*Epl**2*ca**2+ecom*A1q2-Epl**2*A1b2+A1b2*Epl*ca- &
        & Epl**2*sa*ecom*B1q2*ca+ecom*A1q2*Epl**2*ca**2-Epl*sa*ecom*B1q2-Epl*sa*B1b2- &
        & Epl**2*sa*B1b2*ca+ecom*A1q2*Epl*ca+A1b2-Epl**2*ecom*A1q2)/(Epl*ca+1)**2

   r34(3) = 2*Ppl*(-B1b2*Epl**2*ca**2-Epl**2*sa*A1b2*ca+Epl**3*sa*ecom*A1q2*ca**2- &
        & 2*Epl**2*ecom*B1q2+Epl**3*ecom*B1q2*ca**3-2*Epl**3*sa*ecom*A1q2+2*Epl*sa*ecom*A1q2- &
        & 2*Epl**3*ecom*B1q2*ca+ecom*B1q2*Epl*ca+2*Epl**2*B1b2+ecom*B1q2-B1b2- &
        & Epl**3*sa*A1b2*ca**2-2*Epl*sa*A1b2+2*Epl**3*B1b2*ca-B1b2*Epl*ca-Epl**3*B1b2*ca**3+ &
        & 2*Epl**3*sa*A1b2+Epl**2*sa*ecom*A1q2*ca+ecom*B1q2*Epl**2*ca**2)/(Epl*ca+1)**3

   r34(2) = -4*Ppl*(-2*Epl**2*a2B1-2*Epl**3*a2B1*ca+a2B1*Epl*ca+2*Epl*sa*A1a2- &
        & 2*Epl**3*sa*A1a2+Epl**3*a2B1*ca**3+Epl**3*sa*A1a2*ca**2+a2B1+ &
        & a2B1*Epl**2*ca**2+Epl**2*sa*A1a2*ca)/(Epl*ca+1)**3

   r34(1) = 2*Ppl*(-2*Epl**2*B1b2-2*Epl**2*ecom*B1q2+Epl**3*sa*A1b2*ca**2- &
        & 2*Epl**3*ecom*B1q2*ca+ecom*B1q2*Epl*ca+2*Epl*sa*ecom*A1q2- &
        & 2*Epl**3*sa*ecom*A1q2+B1b2*Epl*ca+Epl**3*B1b2*ca**3+Epl**3*ecom*B1q2*ca**3- &
        & 2*Epl**3*B1b2*ca+2*Epl*sa*A1b2-2*Epl**3*sa*A1b2+Epl**3*sa*ecom*A1q2*ca**2+ &
        & ecom*B1q2+B1b2*Epl**2*ca**2+Epl**2*sa*A1b2*ca+ecom*B1q2*Epl**2*ca**2+ &
        & Epl**2*sa*ecom*A1q2*ca+B1b2)/(Epl*ca+1)**3
   
   r35(1) = 8*Epl**2*sa*Ppl*ca*(Epl-1)*(Epl+1)/(-1+Epl*ca)**2/(Epl*ca+1)**2
   
   r36(1) = -4*Epl**2*Ppl*(Epl-1)*(Epl+1)*(Epl*ca-1+2*ca**2)/(-1+Epl*ca)/(Epl*ca+1)**3
   
   r41(3) = 2*Ppl*(B1b2*Epl**2*ca**2+Epl**2*sa*A1b2*ca+Epl**3*sa*ecom*A1q2*ca**2+ &
        & 2*Epl**2*ecom*B1q2+Epl**3*ecom*B1q2*ca**3-2*Epl**3*sa*ecom*A1q2+ &
        & 2*Epl*sa*ecom*A1q2-2*Epl**3*ecom*B1q2*ca+ecom*B1q2*Epl*ca-2*Epl**2*B1b2- &
        & ecom*B1q2+B1b2-Epl**3*sa*A1b2*ca**2-2*Epl*sa*A1b2+2*Epl**3*B1b2*ca-B1b2*Epl*ca- &
        & Epl**3*B1b2*ca**3+2*Epl**3*sa*A1b2-Epl**2*sa*ecom*A1q2*ca-ecom*B1q2*Epl**2*ca**2)/(-1+Epl*ca)**3

   r41(2) = -4*Ppl*(2*Epl**2*a2B1-2*Epl**3*a2B1*ca+a2B1*Epl*ca+2*Epl*sa*A1a2- &
        & 2*Epl**3*sa*A1a2+Epl**3*a2B1*ca**3+Epl**3*sa*A1a2*ca**2-a2B1- &
        & a2B1*Epl**2*ca**2-Epl**2*sa*A1a2*ca)/(-1+Epl*ca)**3

   r41(1) = 2*Ppl*(2*Epl**2*B1b2+2*Epl**2*ecom*B1q2+Epl**3*sa*A1b2*ca**2- &
        & 2*Epl**3*ecom*B1q2*ca+ecom*B1q2*Epl*ca+2*Epl*sa*ecom*A1q2-2*Epl**3*sa*ecom*A1q2 &
        & +B1b2*Epl*ca+Epl**3*B1b2*ca**3+Epl**3*ecom*B1q2*ca**3-2*Epl**3*B1b2*ca+2*Epl*sa*A1b2 &
        & -2*Epl**3*sa*A1b2+Epl**3*sa*ecom*A1q2*ca**2-ecom*B1q2-B1b2*Epl**2*ca**2- &
        & Epl**2*sa*A1b2*ca-ecom*B1q2*Epl**2*ca**2-Epl**2*sa*ecom*A1q2*ca-B1b2)/(-1+Epl*ca)**3

   r45(1) = 4*Epl**2*Ppl*(Epl-1)*(Epl+1)*(-Epl*ca+2*ca**2-1)/(-1+Epl*ca)**3/(Epl*ca+1)
     
   r46(1) = 8*Epl**2*sa*Ppl*ca*(Epl-1)*(Epl+1)/(-1+Epl*ca)**2/(Epl*ca+1)**2

!   write(*,*)'r31(1:3)',r31(1),r31(2),r31(3)

! set to zero the coefficients of higher degree terms
    p0(6:N) = 0.d0 
    p1(6:N) = 0.d0 
    p2(6:N) = 0.d0
 
    q0(4:N) = 0.d0 
    q1(4:N) = 0.d0 
    q2(4:N) = 0.d0 
    q3(4:N) = 0.d0 
    q4(4:N) = 0.d0 
    r31(4:N) = 0.d0
    r32(4:N) = 0.d0
    r33(4:N) = 0.d0
    r34(4:N) = 0.d0
    r41(4:N) = 0.d0

    r35(2:N) = 0.d0
    r36(2:N) = 0.d0
    r45(2:N) = 0.d0
    r46(2:N) = 0.d0

  END SUBROUTINE matrixdat_ta_shift

! ******************************************************************
! ***** COMPUTE COEFFICIENTS OF MODIFIED SYLVESTER MATRIX **********
! ******************************************************************
! *********** written by GIOVANNI F. GRONCHI (Oct.2004) ************
! ==================================================================
  SUBROUTINE compmodsylv16_shift(pp0,pp1,pp2,qq0,qq1,qq2,qq3,qq4,rr31,rr32,rr33, &
         & rr34,rr35,rr36,rr41,rr45,rr46,SSYLV) 
    IMPLICIT NONE 
! Sylvester matrix elements                                         
    COMPLEX*16,INTENT(IN) :: pp0,pp1,pp2
    COMPLEX*16,INTENT(IN) :: qq0,qq1,qq2,qq3,qq4
    COMPLEX*16,INTENT(IN) :: rr31,rr32,rr33,rr34,rr35
    COMPLEX*16,INTENT(IN) :: rr36,rr41,rr45,rr46
    COMPLEX*16,INTENT(OUT) :: SSYLV(6,6) 
! ==================================================================

    SSYLV(1,1) = pp2 
    SSYLV(1,2) = (0.d0,0.d0) 
    SSYLV(1,3) = (0.d0,0.d0) 
    SSYLV(1,4) = (0.d0,0.d0) 
    SSYLV(1,5) = qq4 
    SSYLV(1,6) = (0.d0,0.d0) 
!                                                                       
    SSYLV(2,1) = pp1 
    SSYLV(2,2) = pp2 
    SSYLV(2,3) = (0.d0,0.d0) 
    SSYLV(2,4) = (0.d0,0.d0) 
    SSYLV(2,5) = qq3
    SSYLV(2,6) = qq4
!                                                                       
    SSYLV(3,1) = rr31
    SSYLV(3,2) = rr32
    SSYLV(3,3) = rr33 
    SSYLV(3,4) = rr34
    SSYLV(3,5) = rr35
    SSYLV(3,6) = rr36
!                                                                       
    SSYLV(4,1) = rr41
    SSYLV(4,2) = rr31
    SSYLV(4,3) = rr32
    SSYLV(4,4) = rr33
    SSYLV(4,5) = rr45
    SSYLV(4,6) = rr46
!                                                                       
    SSYLV(5,1) = (0.d0,0.d0) 
    SSYLV(5,2) = (0.d0,0.d0) 
    SSYLV(5,3) = pp0 
    SSYLV(5,4) = pp1 
    SSYLV(5,5) = qq0 
    SSYLV(5,6) = qq1 
!                                                                       
    SSYLV(6,1) = (0.d0,0.d0) 
    SSYLV(6,2) = (0.d0,0.d0) 
    SSYLV(6,3) = (0.d0,0.d0) 
    SSYLV(6,4) = pp0 
    SSYLV(6,5) = (0.d0,0.d0) 
    SSYLV(6,6) = qq0 
    
  END SUBROUTINE compmodsylv16_shift

! ******************************************************************
! GIVEN t COMPUTE s SUCH THAT (s,t) SATISFIES THE SYSTEM p = q = 0
! ******************************************************************
! ************ written by GIOVANNI F. GRONCHI (Sept.2004) **********
! ==================================================================
  SUBROUTINE solvesystem_ta_shift(Vpltil,vcomtil,nroots,wzero,zzero,sflag,hwflag) 
    USE output_control
    IMPLICIT NONE 
    REAL(KIND=8),INTENT(IN) :: Vpltil,vcomtil
    INTEGER,INTENT(INOUT) :: nroots 
    REAL(KIND=8),INTENT(INOUT), DIMENSION(nroots) :: wzero 
! is inout to eliminate components if the corresponding zzero is complex
    REAL(KIND=8),INTENT(OUT), DIMENSION(nroots) :: zzero 
    LOGICAL,INTENT(INOUT) :: sflag(6) ! solving system messages:
!         sflag(1) = true    OK
!         sflag(1) = false   there are two good solutions             
!         sflag(2) = true    OK                                               
!         sflag(2) = false   neither of the two evaluation is close to 0
!         sflag(3) = true    OK(the discriminant of the 2nd-deg poly is non-negative)
!         sflag(3) = false   the z-component of the solution is complex
!         sflag(4) = true    OK    
!         sflag(4) = false   leading coeff. of the 2nd degree pol. is small
!         sflag(5) = true    OK    
!         sflag(5) = false   1st and 2nd coeffs except of 2nd degree poly
!                            are small                       
!         sflag(6) = true    OK    
!         sflag(6) = false   2nd degree poly have all coefficients small
    LOGICAL,INTENT(INOUT) :: hwflag ! hwflag = .true.  OK!
                                    ! hwflag = .false. abs(root)>10^5
!   2nd degree poly smallness parameter
    REAL(KIND=8),PARAMETER :: eps_2deg=1.d-25
!   2nd degree poly bigness parameter
    REAL(KIND=8),PARAMETER :: big_root=1.d25
!   evaluation smallness parameter
!    REAL(KIND=8),PARAMETER :: eps_eval=1.d-3  
    REAL(KIND=8),PARAMETER :: eps_eval=1.d-5  
!   --------- end interface ------------------------------------------
    REAL(KIND=8) :: svcomtil,cvcomtil,sVpltil,cVpltil
    REAL(KIND=8) :: sx,cx,sy,cy
    REAL(KIND=8) :: stmp(2)
    REAL(KIND=8) :: sf,cf,vcom,sfpl(2),cfpl(2),Vpl(2)
    REAL(KIND=8) :: sfplchk(N),cfplchk(N)
    REAL(KIND=8) :: sfplchktmp(N,N),cfplchktmp(N,N)
    REAL(KIND=8) :: sfchk(N),cfchk(N)
    REAL(KIND=8) :: vvcom(N),VVpl
    REAL(KIND=8) :: deg2z,deg1z,deg0z
!    REAL(KIND=8) :: deg2zold,deg1zold,deg0wold
    REAL(KIND=8) :: deg2ztmp,deg0wtmp
!                        
    REAL(KIND=8) :: z,w,evalst(2),evalstdeg2w(2)
!   choice flag ( choice = 0 ------> init value
!               ( choice = 1 ------> stmp(1)
!               ( choice = 2 ------> stmp(2)
    INTEGER :: choice
!   auxiliary
    REAL(KIND=8), DIMENSION(20) :: wzerotmp(N),zzerotmp(N) 
    LOGICAL :: discheck(N) ! check the discriminant of the 2nd-deg poly for each root
                           ! discheck(i) = .true. OK
                           ! discheck(i) = .false. negative discriminant 
    INTEGER :: nrealroots
!   loop indexes                                                      
    INTEGER :: i,j,k 
!   ==================================================================

!    verb_moid=20

! initialization
    discheck(1:N)=.true.
    sfplchktmp(1:N,1:N)=0.d0
    cfplchktmp(1:N,1:N)=0.d0
! initialization
    sfplchk(1:N)=0.d0
    cfplchk(1:N)=0.d0
    sfchk(1:N)=0.d0
    cfchk(1:N)=0.d0                        
    vvcom(1:N)=0.d0
    
    svcomtil = sin(vcomtil)
    cvcomtil = cos(vcomtil)
    sVpltil = sin(Vpltil)
    cVpltil = cos(Vpltil)
    
    DO i = 1,nroots     
       w=wzero(i) 
!       write(*,*)'root number',i,' w=',w
! ==============================================                    
! deg2z(w)*z**2 + deg1z(w)*z + deg0z(w) = 0 
! (p2(w)*z**2 + p1(w)*z + p0(w) = 0)                     
       deg2z = p2(5)*w**4 + p2(4)*w**3 + p2(3)*w**2 + p2(2)*w + p2(1)
       deg1z = p1(5)*w**4 + p1(4)*w**3 + p1(3)*w**2 + p1(2)*w + p1(1)
       deg0z = p0(5)*w**4 + p0(4)*w**3 + p0(3)*w**2 + p0(2)*w + p0(1)
!       write(*,*)'deg2z,deg1z,deg0z',deg2z,deg1z,deg0z

! CHECK SMALLNESS of the COEFFICIENTS
       IF(abs(deg2z).lt.eps_2deg)THEN 
          ! linear equation
          sflag(4) = .false.
          if(verb_moid.ge.20) then
             WRITE(*,*)'small leading coeff of 2nd degree poly',deg2z 
          endif
          IF(abs(deg1z).ge.eps_2deg) THEN
             zzero(i) = -deg0z/deg1z
             discheck(i) = .true. ! dummy value
             GOTO 133
          ELSEIF(abs(deg1z).lt.eps_2deg)THEN 
             sflag(5) = .false.
             if(verb_moid.ge.20) then
                WRITE(*,*)'small also the coeff of the linear term',deg1z 
             endif
             IF(abs(deg0z).lt.eps_2deg)THEN 
                sflag(6) = .false.
                if(verb_moid.ge.20) then
                   WRITE(*,*)'almost completely degenerate &
                        & 2nd degree poly',deg0z,deg1z,deg2z
                endif
             ENDIF
             ! to skip the computation of this root
             GOTO 124 !actually should compute the roots of the 4th-deg poly...
          ENDIF
       ENDIF

! NORMALIZING COEFFICIENTS (the order is IMPORTANT! deg2z must be the last)
       deg0z = deg0z/deg2z
       deg1z = deg1z/deg2z
       deg2z = deg2z/deg2z
!       write(*,*)'deg2z,deg1z,deg0z after normalization:',deg2z,deg1z,deg0z

! check the discriminant
       IF(deg1z**2 - 4.d0*deg0z*deg2z.lt.0.d0) THEN
          discheck(i) = .false.
          sflag(3) = .false.
          if (verb_moid.ge.20) then
             write(*,*)'negative discriminant',deg1z**2 - 4.d0*deg0z*deg2z
          endif
! *****************************************
          GOTO 133 ! skip this root
! *****************************************
       ELSE
          discheck(i) = .true.
          sflag(3) = .true.
       ENDIF

! solutions                                                         
       stmp(1) = (-deg1z + dsqrt(deg1z**2 - 4.d0*deg0z*deg2z))/       &
            &        (2.d0*deg2z)                                              
       stmp(2) = (-deg1z - dsqrt(deg1z**2 - 4.d0*deg0z*deg2z))/       &
            &        (2.d0*deg2z)                                              
                                                                        
       IF((abs(stmp(1)).gt.big_root).or.(abs(stmp(2)).gt.big_root)) THEN
          hwflag = .false.
          if(verb_moid.ge.20) then
             write(*,*)'solvesystem: large roots of 2nd deg poly!',abs(stmp(1:2))
          endif
          GOTO 124
       ENDIF

!     ======================================================            
!     deg4z(w)*z**4 + deg3z(w)*z**3 + deg2z(w)*z**2 + deg1z(w)*z + deg0z(w) = 0
!     (q4(w)*z**4 + q3(w)*z**3 + q2(w)*z**2 + q1(w)*z + q0(w) = 0)

!     conversion into true anomalies                               
       sx = 2.d0*wzero(i)/(1.d0+(wzero(i))**2) 
       cx = (1.d0-(wzero(i))**2)/(1.d0+(wzero(i))**2) 
       sf = sx*cvcomtil + cx*svcomtil
       cf = cx*cvcomtil - sx*svcomtil
       vcom = datan2(sf,cf) 
       vvcom(i) = vcom ! unshifted true anomaly of the comet   
       sfchk(i) = sx 
       cfchk(i) = cx 
!                                                                       
       DO j = 1,2 
          sy =  2.d0*stmp(j)/(1.d0+(stmp(j))**2) 
          cy = (1.d0-(stmp(j))**2)/(1.d0+(stmp(j))**2) 
          sfpl(j) = sy*cVpltil + cy*sVpltil
          cfpl(j) = cy*cVpltil - sy*sVpltil
          Vpl(j) = datan2(sfpl(j),cfpl(j)) 
          VVpl = Vpl(j) ! unshifted true anomaly of the planet
          sfplchktmp(i,j) =  sy 
          cfplchktmp(i,j) =  cy 

          evalstdeg2w(j) = -Ppl*(1.d0+ecom*cos(vvcom(i)))* &
               & ( -ecom*MM*cos(VVpl) - ecom*NN*sin(VVpl) + &
               & sin(vvcom(i))*KK*cos(VVpl) + sin(vvcom(i))*LL*sin(VVpl) - &
               & cos(vvcom(i))*MM*cos(VVpl) - cos(vvcom(i))*NN*sin(VVpl)) -&
               & ecom*pcom*sin(vvcom(i))*(1.d0+Epl*cos(VVpl))
!
          evalst(j) = -pcom*(1.d0+Epl*cos(VVpl))* &
               & ( Epl*LL*cos(vvcom(i)) + Epl*NN*sin(vvcom(i)) + &
               & cos(VVpl)*LL*cos(vvcom(i)) + cos(VVpl)*NN*sin(vvcom(i)) -&
               & sin(VVpl)*KK*cos(vvcom(i)) - sin(VVpl)*MM*sin(vvcom(i))) +&
               & Epl*Ppl*sin(VVpl)*(1.d0+ecom*cos(vvcom(i)))
       ENDDO

!       write(*,*)'************* root(',i,')=',w
!       WRITE(*,*)'EVALSTdeg2w(1)=',evalstdeg2w(1), &
!            &'  EVALSTdeg2w(2)=',evalstdeg2w(2)
!       WRITE(*,*)'EVALST(1)=',evalst(1),'  EVALST(2)=',evalst(2)
!       write(*,*)'--------------------------------------------'

! -------------------------------------------
!     SELECTING the CORRESPONDING SOLUTION                              
       choice = 0 ! initialization 
!     ************ checking smallness of the evaluations *************  
       IF((abs(evalst(1)).lt.eps_eval).or.(abs(evalst(2)).lt.eps_eval))THEN 

!     choosing the one that gives the minimum value                     
          IF (abs(evalst(1)).le.abs(evalst(2))) THEN 
             choice = 1 
             zzero(i) = stmp(1) 
             sfplchk(i) = sfplchktmp(i,1) 
             cfplchk(i) = cfplchktmp(i,1)  

          ELSEIF (abs(evalst(1)).gt.abs(evalst(2))) THEN 
             choice = 2 
             zzero(i) = stmp(2) 
             sfplchk(i) = sfplchktmp(i,2) 
             cfplchk(i) = cfplchktmp(i,2)  

          ENDIF

!     if they are BOTH CLOSE TO ZERO                                    
          IF((abs(evalst(1)).lt.eps_eval).and.(abs(evalst(2)).lt.eps_eval))THEN
!     **************************                                         
             sflag(1) = .false.
!     **************************                                           
             if(verb_moid.ge.20) then
!     writing on screen                                                 
               WRITE(*,*)'SOLVING SYSTEM WARNING: TWO GOOD SOLUTIONS'        
               WRITE(*,*)'EVALST(1)=',evalst(1),'  EVALST(2)=',evalst(2)
            endif
         ENDIF
                                                                        
!     ================== DOUBLE CORRESPONDENCE CHECK ===================
!     if sflag(1)-.false., then do this check for all the following values of i
!     -------------------------------------------------------------     
!     CASE 1 (i=1)                
          IF((.not.sflag(1)).and.(i.eq.1)) THEN 
!     definitively choosing previously chosen value of z
!     -------------------------------------------------------------     
!     CASE 2 (i>1)                    
          ELSEIF((.not.sflag(1)).and.(i.gt.1))THEN 
!     at first choosing previously chosen value of z
             DO k = 1,i-1 
                IF( (abs(sfchk(i)-sfchk(k)).lt.1.d-1).and.     &
                     &  (abs(cfchk(i)-cfchk(k)).lt.1.d-1).and.         &
                     &  (abs(sfplchk(i)-sfplchk(k)).lt.1.d-1).and.     &
                     &  (abs(cfplchk(i)-cfplchk(k)).lt.1.d-1) )THEN      
!                   write(*,*)'OK, MODIFYING zzero to',stmp(2)
                   IF(choice.eq.1) THEN 
                      zzero(i) = stmp(2) 
                      sfplchk(i) = sfpl(2) 
                      cfplchk(i) = cfpl(2) 
                   ELSEIF(choice.eq.2) THEN 
                      zzero(i) = stmp(1) 
                      sfplchk(i) = sfpl(1) 
                      cfplchk(i) = cfpl(1) 
                   ELSE
                      WRITE(*,*)'ERROR: choice=',choice
                      STOP
                   ENDIF
                ELSE 
!     select previous value of z                                              
                ENDIF
             ENDDO
          ENDIF
          
!     neither of the evaluation is close to zero                        
       ELSE 
!     ***********************                                           
          sflag(2) = .false. 
!     ***********************                                           
          if(verb_moid.ge.20) then
             WRITE(*,*)' NEITHER OF THE EVALUATIONS IS CLOSE TO ZERO'   
             WRITE(*,*)' SELECTING THE ONE WITH THE SMALLEST VALUE '    
             WRITE(*,*)'EVALST(1)=',evalst(1),'EVALST(2)=',evalst(2)    
!            WRITE(*,*)'THE SMALLEST IS',min(abs(evalst(1)),            
!     *           abs(evalst(2)))                                       
          endif
!     choosing the one that gives the minimum value                     
          IF (abs(evalst(1)).le.abs(evalst(2))) THEN 
             zzero(i) = stmp(1) 
          ELSEIF (abs(evalst(1)).gt.abs(evalst(2))) THEN 
             zzero(i) = stmp(2) 
          ENDIF
                                                                        
       ENDIF

133    CONTINUE
       
    ENDDO

!     selecting roots with both real components
    nrealroots = 0
    DO i=1,nroots
       IF(discheck(i)) THEN
          nrealroots=nrealroots+1
          wzerotmp(nrealroots)=wzero(i)
          zzerotmp(nrealroots)=zzero(i)
       ENDIF
    ENDDO

! renumbering wzero,zzero,nroots
    DO i = 1,nrealroots
       wzero(i) = wzerotmp(i)
       zzero(i) = zzerotmp(i)
    ENDDO
    nroots = nrealroots
                                                                        
124 CONTINUE ! to skip computations in case of degeneration of 2nd deg poly
             ! or in case of large roots

  END SUBROUTINE solvesystem_ta_shift

! ******************************************************************
! ********** SELECT AMONG MINIMA, MAXIMA and SADLE POINTS **********
! ******************************************************************
! *********** written by GIOVANNI F. GRONCHI (Sept. 2004) **********
! ==================================================================
  SUBROUTINE hess_ta(Vpl,vcom,ans) 
    USE fund_const
    IMPLICIT NONE 
!   Vpl,vcom are passed in radians                                       
    REAL(KIND=8),INTENT(IN) :: Vpl,vcom
    INTEGER,INTENT(OUT) :: ans !  ans = 1: maximum 
!                                 ans = 0: saddle 
!                                 ans= -1: minimum 
!                                 ans= -2: cannot decide
!   -------- end interface -----------------------------------------
    REAL(KIND=8) :: H11,H12,H21,H22,trH,detH 
    REAL(KIND=8) :: Rpl,rcom
    REAL(KIND=8) :: Xpl,Ypl,xcom,ycom
!   eigenvalues                                                       
    REAL(KIND=8) :: lam1,lam2 
    REAL(KIND=8), PARAMETER :: eps=1.d-20 !tolerance parameter
!   ==================================================================
                                                                        
!   HESSIAN MATRIX 
!(we have eliminated some terms that are zero at critical points)

      Rpl = Ppl/(1.d0 + Epl*cos(Vpl))
      rcom = pcom/(1.d0 + ecom*cos(vcom))
      Xpl = Rpl*cos(Vpl)
      Ypl = Rpl*sin(Vpl)
      xcom = rcom*cos(vcom)
      ycom = rcom*sin(vcom)
      
      H11 = (Rpl/Ppl)*( Epl*Rpl*(Rpl/Ppl)*(Epl*Rpl+Xpl) + &
           & Xpl*(KK*xcom+MM*ycom) + Ypl*(LL*xcom+NN*ycom) )

      H12 = (Rpl*rcom/(Ppl*pcom))*( (ecom*rcom+xcom)*(NN*(Epl*Rpl+Xpl) - &
           & MM*Ypl) - ycom*(LL*(Epl*Rpl+Xpl) - KK*Ypl))

      H21 = H12

      H22 = (rcom/pcom)*( ecom*rcom*(rcom/pcom)*(ecom*rcom+xcom) + &
           & xcom*(KK*Xpl+LL*Ypl) + ycom*(MM*Xpl+NN*Ypl) )

!      write(*,*)'H11,H12,H22',H11,H12,H22

!   TRACE OF THE HESSIAN                                              
    trH = H11 + H22 
!   DETERMINANT OF THE HESSIAN                                        
    detH = H11*H22 - H12*H21 
!    write(*,*)'detH,trH',detH,trH

!   check                                                             
    IF (abs(detH).le.0.d0) THEN 
       WRITE(*,*)'!!!!!!!! det(H) = 0 !!!!!!!' 
       ans = -2
       GOTO 10 
    ENDIF
!   compute eigenvalues                                               
    IF ((trH**2-4.d0*detH).ge.0.d0) THEN 
       lam1 = (trH + dsqrt(trH**2 - 4.d0*detH))/2.d0 
       lam2 = (trH - dsqrt(trH**2 - 4.d0*detH))/2.d0 
       IF ((lam1.gt.eps).and.(lam2.gt.eps)) THEN 
          ans = -1 
       ELSEIF ((lam1.lt.-eps).and.(lam2.lt.-eps)) THEN 
          ans = 1 
       ELSEIF ((lam1.gt.eps).and.(lam2.lt.-eps)) THEN 
          ans = 0 
       ELSEIF ((lam1.lt.-eps).and.(lam2.gt.eps)) THEN 
          ans = 0 
       ELSE 
          ans = -2
       ENDIF
    ELSE 
       WRITE(*,*)'negative discriminant: complex eigenvalues!' 
       ans = -2
    ENDIF
    
10  CONTINUE 
    
  END SUBROUTINE hess_ta
  
! ******************************************************************
! ********** SELECT AMONG MINIMA, MAXIMA and SADLE POINTS **********
! **** by EVALUATING D2 in a NEIGHBORHOOD of the CRITICAL POiNT ****
! ******************************************************************
! ************* written by GIOVANNI F. GRONCHI (2003) **************
! ********** Department of Mathematics, UNIVERSITY of PISA *********
! ==================================================================
  SUBROUTINE int_eval_ta(Vpl,vcom,ans)
    USE fund_const
    IMPLICIT NONE 
!   u,upl are passed in radians                                       
    REAL(KIND=8),INTENT(IN) :: Vpl,vcom 
    INTEGER,INTENT(OUT) :: ans !  ans = 1: maximum 
!                                 ans = 0: saddle 
!                                 ans= -1: minimum 
!                                 ans= -2: cannot decide
! ----------- end interface ---------------------------------------
    INTEGER,PARAMETER :: n = 30 
    REAL(KIND=8) :: fpl,fcom
    REAL(KIND=8) :: vD2
    REAL(KIND=8) :: fD2
    INTEGER :: count ! counter
    INTEGER :: j
! =================================================================

    count = 0

    CALL d2eval_ta(Vpl,vcom,vD2) 
!    write(*,*)'uD2',uD2

    DO j = 0,n-1
       fpl = Vpl + radeg*sin(j*dpig/n)
       fcom = vcom + radeg*cos(j*dpig/n)
       CALL d2eval_ta(fpl,fcom,fD2) 
!    write(*,*)'vD2',vD2

       IF(vD2.lt.fD2) THEN
          count = count - 1
       ELSEIF(vD2.gt.fD2) THEN
          count = count + 1
       ELSE
          WRITE(*,*)'equal values of D2'
       ENDIF
    ENDDO

    IF(count.eq.-n) THEN
       ans = -1
    ELSEIF(count.eq.n) THEN
       ans = 1
    ELSE
       ans = 0
    ENDIF

  END SUBROUTINE int_eval_ta

! ******************************************************************
! ****** EVALUATE THE SQUARED DISTANCE IN THE POINTS fpl,fcom ******
! ******************************************************************
! *********** written by GIOVANNI F. GRONCHI (Sept.2004) ***********
! ==================================================================
  SUBROUTINE d2eval_ta(Vpl,vcom,D2) 
    USE fund_const
    IMPLICIT NONE                            
    DOUBLE PRECISION,INTENT(IN) :: Vpl,vcom ! true anomalies (in rad.)
    DOUBLE PRECISION D2 ! SQUARED DISTANCE function 
!   ------------- end interface --------------------------------------
    DOUBLE PRECISION :: Rpl,rcom
    DOUBLE PRECISION :: Xpl,Ypl,xcom,ycom
    DOUBLE PRECISION :: Xpl1,Xpl2,Xpl3
    DOUBLE PRECISION :: xcom1,xcom2,xcom3
    DOUBLE PRECISION :: XX,YY,ZZ
!   ==================================================================

    Rpl = Ppl/(1.d0 + Epl*cos(Vpl))
    rcom = pcom/(1.d0 + ecom*cos(vcom))
    Xpl = Rpl*cos(Vpl)
    Ypl = Rpl*sin(Vpl)
    xcom = rcom*cos(vcom)
    ycom = rcom*sin(vcom)
! space coords of the planet orbit
    Xpl1 = Xpl*PPx + Ypl*QQx
    Xpl2 = Xpl*PPy + Ypl*QQy
    Xpl3 = Xpl*PPz + Ypl*QQz
! space coords of the comet orbit
    xcom1 = xcom*px + ycom*qx
    xcom2 = xcom*py + ycom*qy
    xcom3 = xcom*pz + ycom*qz
! difference vector
    XX = Xpl1 - xcom1
    YY = Xpl2 - xcom2
    ZZ = Xpl3 - xcom3
    D2 = XX**2 + YY**2 + ZZ**2
                                                                        
!   check if D2 i positive; otherwise set it to zero                  
    IF (D2.ge.0) THEN 
!   do nothing                                                        
    ELSEIF (D2.lt.0) THEN 
       D2 = 0.d0 
    ENDIF
    
  END SUBROUTINE d2eval_ta
  
END MODULE critical_points_shift

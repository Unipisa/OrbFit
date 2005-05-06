SUBROUTINE rotmhalfgam1(sgmez,bJ,tau,dtaudsgm)
  DOUBLE PRECISION,INTENT(IN) :: sgmez !sinus(gamma/2)
  DOUBLE PRECISION,DIMENSION(3),INTENT(IN) :: bJ !ang. momentum unit vector
! rotation of -gamma/2 around bJ and its derivative w.r.t. sinus(gamma/2)
  DOUBLE PRECISION,DIMENSION(3,3),INTENT(OUT) :: tau,dtaudsgm 
! ------------- end interface -------------------------------------------
  DOUBLE PRECISION :: normJ,HRoy
  DOUBLE PRECISION :: cI,sI,cO,sO, Incl, Omeg, gmez
  DOUBLE PRECISION :: cgmez,tgmez ! cosinus and tangent of gamma/2
! rotation matrices
  DOUBLE PRECISION, DIMENSION(3,3) :: ROmeg
  DOUBLE PRECISION, DIMENSION(3,3) :: RmOmeg
  DOUBLE PRECISION, DIMENSION(3,3) :: RIncl
  DOUBLE PRECISION, DIMENSION(3,3) :: RmIncl
  DOUBLE PRECISION, DIMENSION(3,3) :: Rmhgam

  DOUBLE PRECISION, DIMENSION(3,3) :: dRmhgam
  DOUBLE PRECISION, DIMENSION(3,3) :: RORI,RmIRmO !auxiliary
  DOUBLE PRECISION :: vsize !functions
! ---------------------------------------------

!  normJ=sqrt(DOT_PRODUCT(bJ,bJ))
  normJ=vsize(bJ)
  HRoy=sqrt(bJ(1)**2+bJ(2)**2)
  IF(HRoy/normJ.lt.epsilon(1.d0)*100)THEN
     cI=1.d0
     sI=0.d0
     cO=1.d0
     sO=0.d0
  ELSE
     cI = bJ(3)/normJ
     sI = HRoy/normJ
     cO = -bJ(2)/HRoy
     sO = bJ(1)/HRoy
  ENDIF

  cgmez = sqrt(1.d0-sgmez**2)
  tgmez = sgmez/cgmez

  Incl=atan2(sI,cI)
  Omeg=atan2(sO,cO)
  gmez=atan2(sgmez,cgmez)
  CALL rotmt(-Omeg,ROmeg,3)
  CALL rotmt(Omeg,RmOmeg,3)
  CALL rotmt(-Incl,RIncl,1)
  CALL rotmt(Incl,RmIncl,1)
  CALL rotmt(gmez,Rmhgam,3)
  RORI=matmul(ROmeg,RIncl)
  RmIRmO=matmul(RmIncl,RmOmeg)
  tau = matmul(Rmhgam,RmIRmO) ! temporary step
  tau = matmul(RORI,tau)
! -----------------------------------------------------------------------
! tau = matmul(ROmeg,matmul(RIncl,matmul(Rmhgam,matmul(RmIncl,RmOmeg))))
! -----------------------------------------------------------------------

! derivative of Rmhgam w.r.t. sgmez 
  dRmhgam(1:3,1:3)=0.d0! initialization
  dRmhgam(1,1)=-tgmez
  dRmhgam(1,2)=1.d0
  dRmhgam(2,1)=-1.d0
  dRmhgam(2,2)=-tgmez
  dtaudsgm = matmul(dRmhgam, RmIRmO) ! temporary step
  dtaudsgm = matmul(RORI, dtaudsgm)
! -----------------------------------------------------------------------------
!  dtaudsgm = matmul(ROmeg,matmul(RIncl,matmul(dRmhgam,matmul(RmIncl,RmOmeg))))
! -----------------------------------------------------------------------------

END SUBROUTINE rotmhalfgam1

SUBROUTINE rotmhalfgam2(sgmez,bJ,tau,dtaudsgm)
  DOUBLE PRECISION,INTENT(IN) :: sgmez !sinus(gamma/2)
  DOUBLE PRECISION,DIMENSION(3),INTENT(IN) :: bJ !ang. momentum unit vector
! rotation of -gamma/2 around bJ and its derivative w.r.t. sinus(gamma/2)
  DOUBLE PRECISION,DIMENSION(3,3),INTENT(OUT) :: tau,dtaudsgm 
! ------------- end interface -------------------------------------------
  DOUBLE PRECISION :: normJ,HRoy
  DOUBLE PRECISION :: cI,sI,cO,sO
  DOUBLE PRECISION :: cgmez,tgmez ! cosinus and tangent of gamma/2
! rotation matrices
  DOUBLE PRECISION, DIMENSION(3,3) :: ROmeg
  DOUBLE PRECISION, DIMENSION(3,3) :: RmOmeg
  DOUBLE PRECISION, DIMENSION(3,3) :: RIncl
  DOUBLE PRECISION, DIMENSION(3,3) :: RmIncl
  DOUBLE PRECISION, DIMENSION(3,3) :: Rmhgam
  DOUBLE PRECISION, DIMENSION(3,3) :: dRmhgam
! ---------------------------------------------

  normJ=sqrt(DOT_PRODUCT(bJ,bJ))
  !normJ=size(bJ)  ????? 
  HRoy=sqrt(bJ(1)**2+bJ(2)**2)
  IF(HRoy/normJ.lt.epsilon(1.d0)*100)THEN
     cI=1.d0
     sI=0.d0
     cO=1.d0
     sO=0.d0
  ELSE
     cI = bJ(3)/normJ
     sI = HRoy/normJ
     cO = -bJ(2)/HRoy
     sO = bJ(1)/HRoy
  ENDIF
  cgmez = sqrt(1.d0-sgmez**2)
  tgmez = sgmez/cgmez

  ROmeg(1:3,1:3)=0.d0
  RmOmeg(1:3,1:3)=0.d0
  RIncl(1:3,1:3)=0.d0
  RmIncl(1:3,1:3)=0.d0
  Rmhgam(1:3,1:3)=0.d0
  dRmhgam(1:3,1:3)=0.d0

  ROmeg(1,1)=cO
  ROmeg(1,2)=-sO
  ROmeg(2,1)=sO
  ROmeg(2,2)=cO
  ROmeg(3,3)=1.d0

  RmOmeg(1,1)=cO
  RmOmeg(1,2)=sO
  RmOmeg(2,1)=-sO
  RmOmeg(2,2)=cO
  RmOmeg(3,3)=1.d0

  RIncl(2,2)=cI
  RIncl(2,3)=-sI
  RIncl(3,2)=sI
  RIncl(3,3)=cI
  RIncl(1,1)=1.d0

  RmIncl(2,2)=cI
  RmIncl(2,3)=sI
  RmIncl(3,2)=-sI
  RmIncl(3,3)=cI
  RmIncl(1,1)=1.d0

  Rmhgam(1,1)=cgmez
  Rmhgam(1,2)=sgmez
  Rmhgam(2,1)=-sgmez
  Rmhgam(2,2)=cgmez
  Rmhgam(3,3)=1.d0

  dRmhgam(1,1)=-tgmez
  dRmhgam(1,2)=1.d0
  dRmhgam(2,1)=-1.d0
  dRmhgam(2,2)=-tgmez

 tau = matmul(ROmeg,matmul(RIncl,matmul(Rmhgam,matmul(RmIncl,RmOmeg))))
 dtaudsgm = matmul(ROmeg,matmul(RIncl,matmul(dRmhgam,matmul(RmIncl,RmOmeg))))

END SUBROUTINE rotmhalfgam2

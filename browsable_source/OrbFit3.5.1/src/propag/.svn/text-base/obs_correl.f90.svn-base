MODULE obs_correl
IMPLICIT NONE

PUBLIC :: obscor

! Models of autocorrelation functions for all the observatories
!
! Max number of parameters/functions (all observatories together)
INTEGER, PARAMETER :: nparx=500
! Min and max value for observatory code
INTEGER,PARAMETER :: obsc1=0,obsc2=2000
! Models of autocorrelation functions for all the observatories
! pto2f       -  pointer from observatory code to function
! pto2fm      -  pointer from observatory code to function for mixed class
! nfo         -  number of functions for the model of each station
! nfom        -  number of functions for the model of mixed class
! nfunt       -  total number of functions
! npart       -  total number of parameters
! nparf       -  number of parameters for each function
! nparo       -  number of parameters for each observatory
! nparom      -  number of parameters for mixed class
! kfun        -  function integer code
! ptf2p       -  pointer from function to parameters
! par         -  parameter values
! aprmx       -  max a-priori RMS (arcsec)
! minapw      -  min a-priori weight (rad**(-2))
! maxdst      -  max "decision step" (1=specific obs; 2=mixed class)
! iiccor      -  initialization check
!
INTEGER nfunt,npart,pto2f(obsc1:obsc2,2,2),nfo(obsc1:obsc2,2)
INTEGER nparf(nparx),kfun(nparx),ptf2p(2,nparx)
INTEGER nparo(obsc1:obsc2,2),maxdst,iiccor
INTEGER pto2fm(2,2),nparom(2),nfom(2)
DOUBLE PRECISION par(nparx),aprmx,minapw

CONTAINS

! Copyright (C) 2000 by Mario Carpino (carpino@brera.mi.astro.it)
! Version: November 28, 2000
! ---------------------------------------------------------------------
!
!  *****************************************************************
!  *                                                               *
!  *                         O B S C O R                           *
!  *                                                               *
!  * Correlation between two observations of the same observatory  *
!  *                                                               *
!  *****************************************************************
!
! INPUT:    ERRMOD    -  Error model
!           OBS1      -  Observation 1
!           W1        -  A-priori weight for RA/DEC (rad**(-2)) (obs 1)
!           OBS2      -  Observation 2
!           W2        -  A-priori weight for RA/DEC (rad**(-2)) (obs 2)
!
! OUTPUT:   COVRA     -  Off-diagonal covariance in RA (rad**2)
!           COVDEC    -  Off-diagonal covariance in DEC (rad**2)
!
SUBROUTINE obscor(errmod,obs1,w1,obs2,w2,covra,covdec)
USE astrometric_observations
IMPLICIT NONE

CHARACTER(LEN=*),               INTENT(IN)  :: errmod
TYPE(ast_obs),                  INTENT(IN)  :: obs1,obs2
DOUBLE PRECISION, DIMENSION(2), INTENT(IN)  :: w1,w2
DOUBLE PRECISION,               INTENT(OUT) :: covra,covdec

DOUBLE PRECISION :: rmsc1,rmsc2,bias1,bias2,dt
INTEGER, DIMENSION(2) :: idst1,idst2
LOGICAL,SAVE :: first=.true.,doit
CHARACTER*80 file
INTEGER :: lf,ic,ip1,if1
INTEGER, EXTERNAL :: lench
CHARACTER*20,SAVE :: errmodkeep
LOGICAL :: mixcl

! Initialization
IF(first) THEN
    first=.false.
    doit=(errmod.ne.' ')
    errmodkeep=errmod
! Input of correlation model
    IF(doit) THEN
        file=errmod//'.cor'
        CALL rmsp(file,lf)
        WRITE(*,210) file(1:lf)
210 FORMAT('INFO(obscor): reading error correlation model from file "',A,'"')
        CALL rdcorm(file)
    ELSE
        WRITE(*,211)
    END IF
END IF

IF(errmod.ne.errmodkeep)THEN
   WRITE(*,*)' obscor: change in errmod from ',errmodkeep,' to ', errmod,' not allowed'
   STOP
ENDIF
211 FORMAT('INFO(obscor): NO error correlation model is used')

covra=0
covdec=0
IF(.NOT.doit) RETURN
IF(iiccor.NE.36) STOP '**** obscor: internal error (01) ****'
IF(obs1%obscod_i /= obs2%obscod_i) RETURN
IF(obs1%type /= obs2%type) RETURN
IF(obs1%type /= 'O') RETURN
IF(obs1%obscod_i.LT.obsc1.OR.obs1%obscod_i.GT.obsc2) STOP '**** obscor: internal error (02) ****'
dt=ABS(obs2%time_tdt-obs1%time_tdt)

! These calls are needed only for computing idst1 and idst2
CALL astrow(errmod,obs1%tech,obs1%time_tdt,obs1%obscod_i,obs1%acc_coord(1),obs1%acc_coord(2),  &
            rmsc1,rmsc2,bias1,bias2,idst1)
CALL astrow(errmod,obs2%tech,obs2%time_tdt,obs2%obscod_i,obs2%acc_coord(1),obs2%acc_coord(2),  &
            rmsc1,rmsc2,bias1,bias2,idst2)

! Right ascension
ic=1
IF(idst1(ic).EQ.0 .OR. idst2(ic).EQ.0) GOTO 1
IF(idst1(ic).GT.maxdst) GOTO 1
IF(idst2(ic).GT.maxdst) GOTO 1
IF(w1(1).LT.minapw) GOTO 1
IF(w2(1).LT.minapw) GOTO 1
if1=pto2f(obs1%obscod_i,ic,1)
mixcl=(if1.LE.0)
IF(mixcl) THEN
   if1=pto2fm(ic,1)
   IF(if1.LE.0) GOTO 1
   ip1=ptf2p(1,if1)
   CALL fcorob(kfun(if1),nfom(ic),par(ip1),nparf(if1),dt,covra)
ELSE
   ip1=ptf2p(1,if1)
   CALL fcorob(kfun(if1),nfo(obs1%obscod_i,ic),par(ip1),nparf(if1),dt,covra)
END IF
covra=covra/SQRT(w1(1)*w2(1))
1 CONTINUE

! Declination
ic=2
IF(idst1(ic).EQ.0 .OR. idst2(ic).EQ.0) GOTO 2
IF(idst1(ic).GT.maxdst) GOTO 2
IF(idst2(ic).GT.maxdst) GOTO 2
IF(w1(2).LT.minapw) GOTO 2
IF(w2(2).LT.minapw) GOTO 2
if1=pto2f(obs1%obscod_i,ic,1)
mixcl=(if1.LE.0)
IF(mixcl) THEN
   if1=pto2fm(ic,1)
   IF(if1.LE.0) GOTO 2
   ip1=ptf2p(1,if1)
   CALL fcorob(kfun(if1),nfom(ic),par(ip1),nparf(if1),dt,covdec)
ELSE
   ip1=ptf2p(1,if1)
   CALL fcorob(kfun(if1),nfo(obs1%obscod_i,ic),par(ip1),nparf(if1),dt,covdec)
END IF
covdec=covdec/SQRT(w1(2)*w2(2))
2 CONTINUE

END SUBROUTINE obscor

! Copyright (C) 2000 by Mario Carpino (carpino@brera.mi.astro.it)
! Version: November 16, 2000
! ---------------------------------------------------------------------
!
!  *****************************************************************
!  *                                                               *
!  *                         F C O R O B                           *
!  *                                                               *
!  *              Computation of correlation function              *
!  *                                                               *
!  *****************************************************************
!
! INPUT:    KFUN      -  Function integer codes
!           NFUN      -  Number of elementary functions
!           PAR       -  Parameter values
!           NPARF     -  Number of parameters per function
!           DT        -  Time lag
!
! OUTPUT:   FUN       -  Function value
!
      SUBROUTINE fcorob(kfun,nfun,par,nparf,dt,fun)
      IMPLICIT NONE

      INTEGER nfun
      INTEGER nparf(nfun),kfun(nfun)
      DOUBLE PRECISION par(*),dt,fun

      INTEGER if1,ip1,kf1
      DOUBLE PRECISION par1,ff1,fd1,term

      fun=0
      term=0

      ip1=1
      DO 1 if1=1,nfun
      kf1=kfun(if1)
      par1=par(ip1)

      INCLUDE 'fcfund.h90'

      ELSE
          STOP '**** fcorob: internal error (01) ****'
      END IF

      IF(kf1.EQ.1) THEN
          fun=fun+term
          term=ff1
      ELSE
          term=term*ff1
      END IF
      ip1=ip1+nparf(if1)

    1 CONTINUE
      fun=fun+term

      END SUBROUTINE fcorob
! Copyright (C) 2000 by Mario Carpino (carpino@brera.mi.astro.it)
! Version: November 28, 2000
! ---------------------------------------------------------------------
!
!  *****************************************************************
!  *                                                               *
!  *                         R D C O R M                           *
!  *                                                               *
!  *   Read the autocorrelation models for all the observatories   *
!  *                                                               *
!  *****************************************************************
!
! INPUT:    FILE      -  Input file
!
! OUTPUT:
!
      SUBROUTINE rdcorm(file)
      USE fund_const
      IMPLICIT NONE

      CHARACTER*(*) file

      INTEGER unit,io,line,ic,lf,ip,ivers,lr
      CHARACTER(LEN=80) :: rec,rec1
      LOGICAL :: mixcl

      INTEGER lench
      EXTERNAL lench

! Initialization
      DO 1 io=obsc1,obsc2
      pto2f(io,1,1)=0
      pto2f(io,2,1)=0
    1 END DO
      pto2fm(1,1)=0
      pto2fm(2,1)=0
      nfunt=0
      npart=0
      line=0

! Reading model
      lf=lench(file)
      CALL filopl(unit,file)
      READ(unit,*) ivers
      IF(ivers.NE.1) THEN
          WRITE(*,201) file(1:lf)
          STOP '**** Abnormal end ****'
      END IF
  201 FORMAT('ERROR(rdcorm): unsupported version of file "',A,'"')
      READ(unit,*) aprmx
      READ(unit,*) maxdst
      minapw=1.D0/((aprmx*radsec)**2)
      line=3

    2 CONTINUE
! Read a new observatory/coordinate
      READ(unit,100,END=10) rec
  100 FORMAT(A)
      line=line+1
! IC = coordinate (1=RA, 2=DEC)
      IF(rec(1:8).EQ.'TIME RA ') THEN
          ic=1
      ELSEIF(rec(1:8).EQ.'TIME DEC') THEN
          ic=2
      ELSE
          GOTO 20
      END IF
! IO = observatory code
      rec1=rec(9:)
      CALL norstr(rec1,lr)
      mixcl=(rec1.EQ.'MIX')
      IF(mixcl) THEN
          pto2fm(ic,1)=nfunt+1
          nfom(ic)=0
          nparom(ic)=0
      ELSE
          READ(rec1,*,ERR=20) io
          IF(io.LT.obsc1) STOP '**** rdcorm: obsc < obsc1 ****'
          IF(io.GT.obsc2) STOP '**** rdcorm: obsc > obsc2 ****'
          pto2f(io,ic,1)=nfunt+1
          nfo(io,ic)=0
          nparo(io,ic)=0
      END IF

    3 CONTINUE
! Read the model for one coordinate (IC) of one observatory (IO)
      READ(unit,100,ERR=20) rec
      line=line+1
      IF(rec.EQ.'END') THEN
          IF(mixcl) THEN
              IF(nfom(ic).LE.0) GOTO 20
          ELSE
              IF(nfo(io,ic).LE.0) GOTO 20
          END IF
          GOTO 2
      END IF
      nfunt=nfunt+1
      IF(nfunt.GT.nparx) STOP '**** rdcorm: nfunt > nparx ****'
      IF(mixcl) THEN
          pto2fm(ic,2)=nfunt
          nfom(ic)=nfom(ic)+1
      ELSE
          pto2f(io,ic,2)=nfunt
          nfo(io,ic)=nfo(io,ic)+1
      END IF
! Determine function integer code (KP1) and number of parameters (NP1)
      CALL fcsfun(rec,kfun(nfunt),nparf(nfunt))
      IF(kfun(nfunt).LE.0) GOTO 20
      IF(mixcl) THEN
          nparom(ic)=nparom(ic)+nparf(nfunt)
      ELSE
          nparo(io,ic)=nparo(io,ic)+nparf(nfunt)
      END IF
! Reading and storing parameter values and properties
      IF(npart+nparf(nfunt).GT.nparx) STOP '**** rdcorm: npart > nparx ****'
      ptf2p(1,nfunt)=npart+1
      DO 4 ip=1,nparf(nfunt)
      npart=npart+1
      READ(unit,*,ERR=20) par(npart)
    4 END DO
      ptf2p(2,nfunt)=npart
      GOTO 3

! Regular end
   10 CONTINUE
      CALL filclo(unit,' ')
      iiccor=36
      RETURN

! Error termination
   20 CONTINUE
      WRITE(*,200) file(1:lf),line
  200 FORMAT('rdcorm: INPUT ERROR from file "',A,'" at record',I5)
      STOP '**** rdcorm: abnormal end ****'

      END SUBROUTINE rdcorm
! Copyright (C) 2000 by Mario Carpino (carpino@brera.mi.astro.it)
! Version: November 3, 2000
! ---------------------------------------------------------------------
!
!  *****************************************************************
!  *                                                               *
!  *                         F C S F U N                           *
!  *                                                               *
!  *              LS fit of covariance functions:                  *
!  *       function integer codes and number of parameters         *
!  *                                                               *
!  *****************************************************************
!
! INPUT:    NAME      -  Function name
!
! OUTPUT:   KFUN      -  Function integer identificator (0=don't know)
!           NPAR      -  Number of parameters
!
      SUBROUTINE fcsfun(name,kfun,npar)
      IMPLICIT NONE

      INTEGER kfun,npar
      CHARACTER*(*) name

      kfun=0
      npar=0
! Multiplicative coefficient (constant)
      IF(name.EQ.'+COEF') THEN
          kfun=1
          npar=1
! Exponential function
      ELSEIF(name.EQ.'*EXP') THEN
          kfun=2
          npar=1
! Normal function
      ELSEIF(name.EQ.'*NORM') THEN
          kfun=3
          npar=1
! Parabola (to be multiplied by exponential or normal functions)
      ELSEIF(name.EQ.'*PARAB') THEN
          kfun=4
          npar=1
! ADD HERE NEW FUNCTIONS, using increasing integer identificator (KFUN)
! Remember to update accordingly also subroutines FCFUND and FCWPAR
      END IF

      END SUBROUTINE fcsfun

END MODULE obs_correl

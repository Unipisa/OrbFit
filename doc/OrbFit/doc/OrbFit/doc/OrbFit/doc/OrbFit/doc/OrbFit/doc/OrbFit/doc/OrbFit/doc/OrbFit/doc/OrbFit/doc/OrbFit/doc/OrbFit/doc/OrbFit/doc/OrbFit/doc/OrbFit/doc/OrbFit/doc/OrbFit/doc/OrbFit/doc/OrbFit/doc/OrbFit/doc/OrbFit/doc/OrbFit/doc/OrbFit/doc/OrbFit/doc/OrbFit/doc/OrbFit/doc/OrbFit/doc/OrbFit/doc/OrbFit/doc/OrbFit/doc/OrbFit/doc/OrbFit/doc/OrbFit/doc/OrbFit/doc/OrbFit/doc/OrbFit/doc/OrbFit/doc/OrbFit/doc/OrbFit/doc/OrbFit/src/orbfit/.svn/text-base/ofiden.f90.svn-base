!
!  *****************************************************************
!  *                                                               *
!  *                         O F I D E N                           *
!  *                                                               *
!  *                Auxiliary routine for ORBFIT:                  *
!  *                     orbit identification                      *
!  *                                                               *
!  *****************************************************************
!
! IN/OUT:   UNIREP    -  FORTRAN unit for report
!           UNIELE    -  FORTRAN unit for output of orbital elements
!           UNIDIF    -  FORTRAN unit for logfilr of DIFCOR
!           OPELE     -  Flag (need to open UNIELE)
!           ELEOUT    -  Name of file for output of orbital elements
!           MAGOPT    -  Magnitude determination option
!           NAME      -  Object names
!           NAMOF     -  Names of observation files (w/o extension)
!           DEFORB    -  Orbit definition flag
!           DEFCN     -  Tells whether covariance/normal matrices
!                            are defined
!           ELEM      -  Orbital elements
!           ELEM_UNC  -  Orbital element uncertainty
!           MASS      -  Object masses
!           COMELE    -  Comment on orbital elements
!           DIR       -  Directory containing observation/residual file
!           NOBJ      -  Number of objects
!           OBS, OBSW -  Observations
!           N         -  Number of observations for each object
!           NT        -  Total number of observations
!           IP1       -  Pointer to first observation for each object
!           IP2       -  Pointer to first weight for each object
!           GMSUN     -  G*Mass(Sun)
!           OETYPE    -  Type of output elements (CAR/EQU/KEP/EQP)
!
      SUBROUTINE ofiden(unirep,uniele,unidif,opele,eleout,magopt,name,namof,deforb,defcn,elem,elem_unc,         &
     &                  mass,comele,dir,nobj,obs,obsw,n,nt,ip1,ip2,gmsun,oetype,error_model)
      USE astrometric_observations
      USE least_squares
      USE orbit_elements
      IMPLICIT NONE

      INCLUDE 'parnob.h90'
      INCLUDE 'parobx.h90'
      INCLUDE 'parcmc.h90'

      INTEGER unirep,uniele,unidif,nobj,magopt
      INTEGER ip1(nobjx),ip2(nobjx)
! OBSERVATIONS
      INTEGER nt,n(nobj)
! new data types
      TYPE(ast_obs),DIMENSION(nt) :: obs
      TYPE(ast_wbsr),DIMENSION(nt) :: obsw
      TYPE(orbit_elem),DIMENSION(nobj1x) :: elem
      TYPE(orb_uncert),DIMENSION(nobj1x) :: elem_unc

      CHARACTER(LEN=*) ::  error_model ! weighing model

      DOUBLE PRECISION gmsun
      DOUBLE PRECISION mass(nobj1x)
      LOGICAL opele,deforb(nobj1x),defcn(nobj1x)
      CHARACTER*(*) eleout,name(nobj1x),namof(nobj1x)
      CHARACTER*(*) dir(nobj1x),comele(nobj1x),oetype

      INTEGER i,k,ng,lnt,ln1,icor(6),nsm,fail_flag
      TYPE(orbit_elem),DIMENSION(nobjx) :: elemp
      DOUBLE PRECISION hnew,rmsh
      DOUBLE PRECISION csinor,delnor,elemo(6),telemo
      DOUBLE PRECISION telemt,elemt(6),masst,ht,gt,ennet,gma1,enne
      DOUBLE PRECISION gtwg(6,6),gmagc
      TYPE(orbit_elem) :: elem1
      CHARACTER namet*100,titnam*80,file*150
      LOGICAL ok,autrep

! NEEDED common blocks:
     INCLUDE 'comidn.h90'
!      INCLUDE 'comrej.h'

! TUNING PARAMETERS
! Small time (d)
      DOUBLE PRECISION epst
      PARAMETER (epst=1.0D-9)

      INTEGER lench
      EXTERNAL lench

      IF(nobj.LT.2) RETURN
      DO 1 i=1,nobj
      IF(.NOT.(deforb(i).AND.defcn(i))) RETURN
    1 END DO

      nobj=2
      IF(nobj.GT.nobjx) STOP '**** ofiden: internal error (03) ****'
      IF(nobj1x.LT.3) STOP '**** ofiden: internal error (04) ****'

      WRITE(unirep,207)
  207 FORMAT('Orbit identification:')

      DO 3 i=1,nobj
!      gma1=gmsun*(1+mass(i))
      CALL coo_cha(elem(i),'EQU',elemp(i),fail_flag)
      ln1=lench(name(i))
      WRITE(unirep,200) name(i)(1:ln1)
      CALL outele(unirep,elemp(i)%coord,'EQU',elemp(i)%t,' ',.true.,.false.)
      IF(magopt.NE.0 .AND. elem(i)%mag_set .AND. elem(i)%h_mag.GT.-100.D0) WRITE(unirep,140) elem(i)%h_mag
    3 END DO
  200 FORMAT(5X,'Starting orbital elements for object ',A,':')
  140 FORMAT(10X,'Absolute mag. H    =',F9.2)

! Computation of starting elements (average of single objects)
      telemt=0
      masst=0
      DO 4 k=2,5
      elemt(k)=0
    4 END DO
      ht=0
      gt=0
      nsm=0
      DO 5 i=1,nobj
      telemt=telemt+elemp(i)%t
      masst=masst+mass(i)
      DO 6 k=2,5
      elemt(k)=elemt(k)+elemp(i)%coord(k)
    6 END DO
      IF(elem(i)%mag_set) THEN
         IF(elem(i)%h_mag.GT.-100.d0) THEN
            ht=ht+elem(i)%h_mag
            gt=gt+elem(i)%g_mag
            nsm=nsm+1
         END IF
      END IF
    5 END DO
      telemt=telemt/nobj
      masst=masst/nobj
      DO 7 k=2,5
      elemt(k)=elemt(k)/nobj
    7 END DO
      IF(nsm.GT.0) THEN
          ht=ht/nsm
          gt=gt/nsm
      ELSE
          ht=-1.D9
          gt=0.d0
      END IF
      gma1=gmsun*(1+masst)
      CALL start(elemp(1),elemp(2),1,ng,ennet,elemt(1),elemt(6),ok)
      IF(.not.ok) THEN
         WRITE(*,*)' problem in start (initial guess) '
         STOP
      ENDIF
      CALL titast(3,namof(1),namof(2),titnam,namet,lnt)
      WRITE(unirep,200) namet(1:lnt)
      CALL outele(unirep,elemt,'EQU',telemt,' ',.true.,.false.)
      IF(magopt.NE.0 .AND. ht.GT.-100.D0) WRITE(unirep,140) ht
      WRITE(unirep,206) namet(1:lnt)
  206 FORMAT('Differential correction for object ',A,':')
      dir(3)=' '

! Preliminary 2-parameter fit (only semimajor axis and mean anomaly)
      IF(amfit) THEN
          DO 11 i=2,5
          icor(i)=0
   11     CONTINUE
          icor(1)=1
          icor(6)=1
! Disables temporarily outlier rejection
          autrep=autrej
          autrej=.false.
          IF(magopt.NE.0 .AND. ht.GT.-100.D0) gmagc=gt
          elem1=undefined_orbit_elem
          elem1%t=telemt
          elem1%coo='EQU'
          elem1%coord=elemt
          elem1%h_mag=ht
          elem1%g_mag=gmagc
          elem1%mag_set=(elem1%h_mag > -100.D0)
          CALL diff_cor(nt,obs,obsw,elem1,icor,unidif,elem(3),elem_unc(3),csinor,delnor,ok)
          autrej=autrep
          IF(ok) THEN
              elemt(1)=elem(3)%coord(1)
              elemt(6)=elem(3)%coord(6)
              WRITE(unirep,125)
              CALL outele(unirep,elemt,'EQU',telemt,' ',.true.,.false.)
              WRITE(unirep,230) csinor,delnor
          ELSE
              WRITE(unirep,202)
              RETURN
          END IF
      END IF

      DO 10 i=1,6
      icor(i)=1
   10 END DO
      IF(magopt.NE.0 .AND. ht.GT.-100.D0) gmagc=gt
      elem1=undefined_orbit_elem
      elem1%t=telemt
      elem1%coo='EQU'
      elem1%coord=elemt
      elem1%h_mag=ht
      elem1%g_mag=gmagc
      elem1%mag_set=(elem1%h_mag > -100.D0)
      CALL diff_cor(nt,obs,obsw,elem1,icor,unidif,elem(3),elem_unc(3),csinor,delnor,ok)
      IF(ok) THEN
          deforb(3)=.true.
          defcn(3)=.true.
          comele(3)='computed with least squares fit'
! Magnitude estimation
          IF(magopt.NE.0) THEN
              CALL mag_est(nt,obs,obsw,hnew,rmsh)
              IF(rmsh.GT.0.d0) ht=hnew
          END IF
          elem(3)%h_mag=ht
          elem(3)%g_mag=gt
          elem(3)%mag_set=(elem(3)%h_mag > -100.D0)
          name(3)=namet
          namof(3)=namet
          mass(3)=masst
          nobj=3
          IF(opele) THEN
              OPEN(uniele,FILE=eleout,STATUS='UNKNOWN')
              opele=.false.
              WRITE(uniele,301) comcha
              CALL wromlh(uniele,'ECLM','J2000')
          END IF
          CALL coocha(elem(3)%coord,elem(3)%coo,gma1,elemo,oetype,enne)
          telemo=elem(3)%t
          WRITE(unirep,124)
          CALL outele(unirep,elemo,oetype,telemo,' ',.true.,.false.)
          IF(magopt.NE.0 .AND. elem(3)%h_mag.GT.-100.D0) WRITE(unirep,140) elem(3)%h_mag
          WRITE(unirep,230) csinor,delnor
          CALL wromlr(uniele,name(3),elem(3)%coord,elem(3)%coo,elem(3)%t,     &
     &                elem_unc(3)%g,defcn(3),elem_unc(3)%c,defcn(3),        &
     &                elem(3)%h_mag,elem(3)%g_mag,mass(3))
          file=namet(1:lnt)//'.rwo'

          CALL write_rwo(file,obs,obsw,nt,error_model,csinor,rmsh)
      ELSE
          WRITE(unirep,202)
      END IF
  202 FORMAT(5X,'FAILED')
  301 FORMAT(A,' Orbits computed with differential corrections')
  124 FORMAT(5X,'Corrected orbital elements:')
  125 FORMAT(5X,'Intermediate orbital elements (a-M correction):')
  230 FORMAT(5X,'Residual norm   =',1P,E11.3/                           &
     &       5X,'Correction norm =',1P,E11.3)

      END SUBROUTINE ofiden

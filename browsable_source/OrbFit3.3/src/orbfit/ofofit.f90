!
!  *****************************************************************
!  *                                                               *
!  *                         O F O F I T                           *
!  *                                                               *
!  *                Auxiliary routine for ORBFIT:                  *
!  *                  least squares orbital fit                    *
!  *                                                               *
!  *****************************************************************
!
! IN/OUT:   UNIREP    -  FORTRAN unit for report
!           UNIELE    -  FORTRAN unit for output of orbital elements
!           UNIDIF    -  FORTRAN unit for logfilr of DIFCOR
!           OPELE     -  Flag (need to open UNIELE)
!           ELEOUT    -  Name of file for output of orbital elements
!           ORBOPT    -  Orbit determination option
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
!           ERROR_MODEL - Error model
!           OEPTIM    -  Epoch of output elements (MJD, TDT)
!           OEPSET    -  Flag stating that an output epoch is requested
!           OETYPE    -  Type of output elements (CAR/EQU/KEP/EQP)
!
      SUBROUTINE ofofit(unirep,uniele,unidif,opele,eleout,orbopt,magopt,name,namof,deforb,defcn,elem,elem_unc, &
     &                  mass,comele,dir,nobj,obs,obsw,n,nt,ip1,ip2,gmsun,error_model,oeptim,oepset,oetype)
      USE fund_const
      USE astrometric_observations
      USE orbit_elements
      USE least_squares
      USE propag_state
      IMPLICIT NONE

! objects
      INTEGER unirep,uniele,unidif,orbopt,magopt,nobj,nt
      INTEGER ip1(nobj),ip2(nobj),sel(nt),iobs(nt),obscod(nt)
      DOUBLE PRECISION mass(nobj),gmsun
      DOUBLE PRECISION oeptim
      LOGICAL deforb(nobj),defcn(nobj),opele,oepset
      INTEGER fail
      CHARACTER(LEN=*) name(nobj),namof(nobj),eleout,dir(nobj)
      CHARACTER(LEN=*) comele(nobj),oetype
! OBSERVATIONS
      INTEGER n(nobj),nfer
! new data types
      TYPE(ast_obs),DIMENSION(nt) :: obs
      TYPE(ast_wbsr),DIMENSION(nt) :: obsw
      TYPE(orbit_elem),DIMENSION(nobj) :: elem
      TYPE(orb_uncert),DIMENSION(nobj) :: elem_unc

      CHARACTER(LEN=*) ::  error_model ! weighing model

      INCLUDE 'parobx.h90'
      INCLUDE 'parcmc.h90'

! NEEDED common blocks:
      INCLUDE 'comlsf.h90'

! TUNING PARAMETERS
! Max RMS of a "good" observation (arcsec)
      DOUBLE PRECISION maxrms
      PARAMETER (maxrms=5.d0)
! Small time (d)
      DOUBLE PRECISION epst
      PARAMETER (epst=1.0D-9)

      INTEGER i,ln,icor(6),j1,j2,ld
      DOUBLE PRECISION gma1,elemt(6),enne,csinor,delnor,rms1,t1,t2
      DOUBLE PRECISION dxde(6,6),telemp
!      DOUBLE PRECISION lsfres(nob2x),w(nob2x)
      DOUBLE PRECISION covep(6,6),norep(6,6),newep,elemp(6),rmsh
      DOUBLE PRECISION hnew
      CHARACTER file*150
      LOGICAL doit,ok,chep,error,defcov,defnor

      TYPE(orbit_elem) :: elem0,elem1,elem_print
      TYPE(orb_uncert) :: elem_unc_print

!      INCLUDE 'mag.h'

      INTEGER lench
      EXTERNAL lench

      IF(iiclsf.NE.36) STOP '**** ofofit: internal error (01) ****'
      rms1=maxrms*radsec

      DO 10 i=1,6
      icor(i)=1
   10 END DO

      DO 1 i=1,nobj
      IF(.NOT.deforb(i)) GOTO 1
! This check is not enough (many things could have changed since the las
! orbit determination)
!**   IF(orbopt.EQ.1) THEN
!**       doit=(.NOT.defcn(i))
!**   ELSE
          doit=.true.
!**   END IF
      IF(.NOT.doit) GOTO 1
      CALL errmod_set(error_model)
      gma1=gmsun*(1+mass(i))

! Selection of a suitable time of initial conditions for the fit
! Useful timespan of observations (discarding observation of low
! accuracy)
      CALL ustsp(obs(ip1(i):ip1(i)+n(i)-1)%time_tdt,  &
     &           obsw(ip1(i):ip1(i)+n(i)-1)%rms_coord(1), &
     &           obsw(ip1(i):ip1(i)+n(i)-1)%rms_coord(2), &
     &           obsw(ip1(i):ip1(i)+n(i)-1)%sel_coord,     &
     &           n(i),rms1,t1,t2)
      chep=.false.
      newep=elem(i)%t
! Use the time requested for output of orbital elements, but only
! if it is within the useful timespan of observations
      IF(oepset) THEN
          IF(oeptim.GE.t1 .AND. oeptim.LE.t2) THEN
              chep=.true.
              newep=oeptim
          END IF
      END IF
! Otherwise, use the end-point of the observation timespan nearest
! to the requested time
      IF(.NOT.chep) THEN
          IF(elem(i)%t.LT.t1) THEN
              chep=.true.
              newep=t1
          ELSEIF(elem(i)%t.GT.t2) THEN
              chep=.true.
              newep=t2
          END IF
      END IF

! Propagation to new epoch
! ELEM0 are the orbital elements to be used as a starting point for differential correction
      IF(chep .AND. ABS(newep-elem(i)%t).GT.epst) THEN
!         CALL coocha(elem(i)%coord,elem(i)%coo,gma1,elemn,'EQU',enne)
!         CALL proele('EQU',elem(i)%t,elemn,newep,elemt)
          CALL pro_ele(elem(i),newep,elem0)
      ELSE
!         CALL coocha(elem(i)%coord,elem(i)%coo,gma1,elemt,'EQU',enne)
          elem0=elem(i)
      END IF

      ln=lench(name(i))
      WRITE(unirep,206) name(i)(1:ln)
  206 FORMAT('Differential correction for object ',A,':')
      WRITE(unirep,123)
  123 FORMAT(5X,'Starting orbital elements:')
      CALL outele(unirep,elem0%coord,elem0%coo,elem0%t,' ',.true.,.false.)
      IF(magopt.NE.0 .AND. elem(i)%mag_set) THEN
          IF(elem0%mag_set .AND. (elem0%h_mag > -100.D0)) WRITE(unirep,140) elem0%h_mag
      END IF
  140 FORMAT(10X,'Absolute mag. H    =',F9.2)
      CALL diff_cor(n(i),obs(ip1(i):ip1(i)+n(i)-1),obsw(ip1(i):ip1(i)+n(i)-1),  &
     &              elem0,icor,unidif,elem1,elem_unc(i),csinor,delnor,ok)

      IF(ok) THEN
         elem(i)=elem1
         defcn(i)=.true.
         comele(i)='computed with least squares fit'
! Magnitude estimation
         IF(magopt.NE.0) THEN
!           CALL mag_est(m,obs,obsw,h0,rmsh)
            CALL mag_est(n(i),obs(ip1(i):ip1(i)+n(i)-1),obsw(ip1(i):ip1(i)+n(i)-1),hnew,rmsh)
            IF(rmsh.GT.0.d0) elem(i)%h_mag=hnew
         ELSE
            rmsh=0.d0
         END IF
         IF(opele) THEN
            OPEN(uniele,FILE=eleout,STATUS='UNKNOWN')
            opele=.false.
            WRITE(uniele,301) comcha
            CALL wromlh(uniele,'ECLM','J2000')
         END IF

         CALL coo_cha(elem(i),oetype,elem_print,fail,dxde)
!        IF(fail /= 0) THEN
!           elem_print=elem(i)
!           elem_unc_print=elem_unc(i)
!        ELSE
            CALL convertunc(elem_unc(i),dxde,elem_unc_print)
!        END IF
         WRITE(unirep,124)
         CALL outele(unirep,elem_print%coord,elem_print%coo,elem_print%t,' ',.true.,.false.)
         IF(magopt.NE.0 .AND. elem(i)%mag_set .AND. (elem(i)%h_mag > -100.D0)) WRITE(unirep,140) elem(i)%h_mag
         WRITE(unirep,230) csinor,delnor
         CALL wromlr(uniele,name(i),elem_print%coord,elem_print%coo,elem_print%t,elem_unc_print%g,.true.,elem_unc_print%c,   &
                     .true.,elem(i)%h_mag,elem(i)%g_mag,mass(i))
         ld=lench(dir(i))
         ln=lench(namof(i))
         IF(ld.GT.0) THEN
            file=dir(i)(1:ld)//namof(i)(1:ln)//'.rwo'
         ELSE
            file=namof(i)(1:ln)//'.rwo'
         END IF

!        CALL write_rwo(rwofi0,obs,obsw,m,error_model_priv,csino0,rmsh)
         CALL write_rwo(file,obs(ip1(i):ip1(i)+n(i)-1), obsw(ip1(i):ip1(i)+n(i)-1),n(i),error_model,csinor,rmsh)
      ELSE
         WRITE(unirep,202)
      END IF

    1 END DO
  124 FORMAT(5X,'Corrected orbital elements:')
  301 FORMAT(A,' Orbits computed with differential corrections')
  202 FORMAT(5X,'FAILED')
  230 FORMAT(5X,'Residual norm   =',1P,E11.3/                           &
     &       5X,'Correction norm =',1P,E11.3)

      END SUBROUTINE ofofit

!
!  *****************************************************************
!  *                                                               *
!  *                 O B S E R V _ R M S                           *
!  *                                                               *
!  *         Computation of a-priori RMS of observations           *
!  *    from their accuracy (= number of digits in input file)     *
!  *                                                               *
!  *****************************************************************
!
! INPUT:    OBS          - observations
!           ERROR_MODEL  - name of weighing scheme to be used
!           N            - number of observations
!           INIT         - reinitialize all the structure (when init=.false.,
!                          information different from weight/bias is preserved)
!
! OUTPUT:   OBSW         - observation weights
!
      SUBROUTINE observ_rms(obs,error_model,init,obsw,n)
      USE astrometric_observations
      USe fund_const
      IMPLICIT NONE

      INTEGER,                      INTENT(IN)    ::  n
      LOGICAL,                      INTENT(IN)    :: init
      TYPE(ast_obs), DIMENSION(n),  INTENT(IN)    :: obs
      CHARACTER*(*),                INTENT(IN)    :: error_model
      TYPE(ast_wbsr), DIMENSION(n), INTENT(INOUT) :: obsw ! when init=.false.,
                   ! information different from weight/bias is preserved

      INTEGER i,decstep(2)
      DOUBLE PRECISION rmsrmi,rmsvmi
! Get default magnitude weights from a function
      DOUBLE PRECISION magrms

      DO 1 i=1,n
         IF(init) obsw(i)=undefined_ast_wbsr
         IF(obs(i)%type.eq.'O'.or.obs(i)%type.eq.'S')THEN
! astrometric (telescope) observations; get weight
            CALL astrow(error_model,obs(i)%tech,obs(i)%time_tdt,obs(i)%obscod_i,      &
     &               obs(i)%acc_coord(1),obs(i)%acc_coord(2),            &
     &               obsw(i)%rms_coord(1),obsw(i)%rms_coord(2),          &
     &               obsw(i)%bias_coord(1),obsw(i)%bias_coord(2),decstep)
! magnitude weigthing
            obsw(i)%rms_mag=magrms(obs(i)%mag_str,obs(i)%time_tdt,obs(i)%obscod_i,obs(i)%tech)
         ELSEIF(obs(i)%type.eq.'R')THEN
! weighting of radar data is given with observations, but there is
! minimum credible (for a given acccuracy of the models)
            rmsrmi=1.d-3/aukm
            obsw(i)%rms_coord(1)=MAX(obs(i)%acc_coord(1),rmsrmi)
            obsw(i)%rms_coord(2)=9.d9
            obsw(i)%rms_mag=9.d9
            obsw(i)%bias_coord(2)=0.d0
         ELSEIF(obs(i)%type.eq.'V')THEN
            rmsvmi=1.d-3/aukm
            obsw(i)%rms_coord(2)=MAX(obs(i)%acc_coord(2),rmsvmi)
            obsw(i)%rms_coord(1)=9.d9
            obsw(i)%rms_mag=9.d9
            obsw(i)%bias_coord(1)=0.d0
         ELSE
            WRITE(*,*)'observ_rms: obs. type', obs(i)%type,' not known ', &
     &           'at time ',obs(i)%time_tdt,' number ',i
            STOP
        ENDIF
    1 END DO
      END SUBROUTINE observ_rms

! Copyright (C) 1997-1999 by Mario Carpino (carpino@brera.mi.astro.it)
! Version: November 26, 1999
! ---------------------------------------------------------------------
!
!  *****************************************************************
!  *                                                               *
!  *                         A S T R O W                           *
!  *                                                               *
!  *     Computation of a-priori RMS of optical observations       *
!  *         from results of a statistical analysis or             *
!  *    from their accuracy (= number of digits in input file)     *
!  *                                                               *
!  *****************************************************************
!
! INPUT:    MPCTYP    -  Observation type (column 15 of MPC record)
!           TDT       -  Time of observation (MJD, TDT)
!           IDSTA     -  Observatory code 9numeric)
!           ACCA      -  Accuracy of right ascension (rad)
!           ACCD      -  Accuracy of declination (rad)
!
! OUTPUT:   RMSA      -  A-priori RMS of right ascension (rad)
!           RMSD      -  A-priori RMS of declination (rad)
!           BIASA     -  Bias in  right ascension (rad)
!           BIASD     -  Bias in  declination (rad)
!           DECSTEP   -  Decision steps of RMS assignment
!
      SUBROUTINE astrow(error_model,mpctyp,tdt,idsta,acca,accd,        &
     &        rmsa,rmsd,biasa,biasd,decstep)
      USE fund_const
      IMPLICIT NONE
      CHARACTER*(*),INTENT(IN)    :: error_model
      INTEGER, INTENT(IN):: idsta
      DOUBLE PRECISION, INTENT(IN) :: tdt,acca,accd
      CHARACTER*(*), INTENT(IN) :: mpctyp

      DOUBLE PRECISION, INTENT(OUT) :: rmsa,rmsd,biasa,biasd
      INTEGER, DIMENSION(2), INTENT(OUT) :: decstep

      INTEGER le
      DOUBLE PRECISION rmsmin
      CHARACTER*5 ads(2)
      CHARACTER*80 filea,filed,ermnam
      LOGICAL ermuse,error

! Additional info on how RMS are obtained
! ---------------------------------------------------------------------
! former ASTROW.H
! RMS classes of astrometric observations: additional information
!
! idcl        -  Class progressive number in index
! orstep      -  Decision step (RA/DEC) :
!                    0 = no information
!                    1 = single-station class
!                    2 = multi-station class
!                    3 = default
! tdtlim      -  Class limits (MJD, TDT)
!
      INTEGER idcl(2),orstep(2)
!      COMMON/cmaow1/idcl,orstep
      DOUBLE PRECISION tdtlim(2,2)
!      COMMON/cmaow2/tdtlim

      INTEGER lench
      EXTERNAL lench
      LOGICAL first
      DATA first/.true./
      SAVE first
! decide on use of error model file
      ermuse=(error_model.ne.' ')
! Input of RMS class definitions
      IF(first.and.ermuse) THEN
         le=lench(error_model)
         filea=error_model(1:le)//'.cla'
         filed=error_model(1:le)//'.cld'
         CALL rrmscl(filea,filed,.false.)
         first=.false.
      ENDIF
! Negative values means: not yet assigned
      rmsa=-1.D0
      rmsd=-1.D0
      biasa=0.d0
      biasd=0.d0
      orstep(1)=0
      orstep(2)=0
! STEPs 1/2: look for a specific RMS class for the observatory
      IF(ermuse) THEN
          CALL accstr(acca,accd,ads(1),ads(2),error)
          IF(error) GOTO 2
          CALL crmscl(idsta,ads,mpctyp,tdt,rmsa,rmsd,biasa,biasd,idcl,orstep,tdtlim)
          IF(orstep(1).GT.0)THEN
               rmsa=rmsa*radsec
               biasa=biasa*radsec
          ENDIF
          IF(orstep(2).GT.0)THEN
               rmsd=rmsd*radsec
               biasd=biasd*radsec
          ENDIF
      ENDIF

! STEP 3: default, rough rule-of-thumb error model
    2 CONTINUE
      IF(rmsa.LT.0.D0 .OR. rmsd.LT.0.D0) THEN
! Deafault value of RMS (time dependent)
! Before 1890
          IF(tdt.LT.11368.d0) THEN
              rmsmin=3.d0
! From 1890 to 1950
          ELSE IF(tdt.LT.33282.d0) THEN
              rmsmin=2.d0
! After 1950
          ELSE
              rmsmin=1.d0
          END IF
          rmsmin=rmsmin*radsec
          IF(rmsa.LT.0.D0) THEN
              rmsa=MAX(rmsmin,acca)
              orstep(1)=3
          END IF
          IF(rmsd.LT.0.D0) THEN
              rmsd=MAX(rmsmin,accd)
              orstep(2)=3
          END IF
      END IF
      decstep=orstep

      END SUBROUTINE astrow
!***********************************************************************
! 'magrms' returns the default magnitude rms based on the MPC obs string
!***********************************************************************
      DOUBLE PRECISION FUNCTION magrms(magstr,tdt,idsta,typ)
      IMPLICIT NONE

      CHARACTER*6 magstr
! obs. type (from column 15 of .obs format), station code, time MJD
      CHARACTER*1 typ
      INTEGER idsta
      DOUBLE PRECISION tdt
      INTEGER ll,lench
! magnitude weighting for now is simple;
! should be based on digits given,color
      ll=lench(magstr)
      IF(ll.le.0)THEN
         magrms=-1.d0
         RETURN
      ENDIF
      IF(magstr(3:5).eq.'   ')THEN
         magrms=1.0
      ELSEIF(magstr(5:5).eq. ' ')THEN
         magrms=0.7
      ELSE
         magrms=0.5
      ENDIF

      RETURN
      END  FUNCTION magrms
! Copyright (C) 1999-2000 by Mario Carpino (carpino@brera.mi.astro.it)
! Version: November 28, 2000
! ---------------------------------------------------------------------
!
!  *****************************************************************
!  *                                                               *
!  *                         R R M S C L                           *
!  *                                                               *
!  *            Read from a file and store RMS classes             *
!  *                                                               *
!  *****************************************************************
!
! INPUT:    FILRA     -  File name (RA classes)
!           FILDEC    -  File name (DEC classes)
!           CLONLY    -  Load only class list, not accuracy description
!                        (this option is enabled when the routine is
!                        used by stat2.f)
!
! WARNING: the correct behaviour of the routine depends strictly on the
! sorting of input data. More precisely:
!   - "per-observatory" RMS classes must be sorted in increasing order
!      of observatory code (field 3);
!   -  "mixed" (all observatories together) classes must be sorted
!      according to the following sorting keys:
!         key 1: accuracy class (field 1)
!         key 2: type of observation (field 2)
!         key 3: starting time (field 4)
!
      SUBROUTINE rrmscl(filra,fildec,clonly)
      USE station_coordinates
      IMPLICIT NONE

      CHARACTER*(*) filra,fildec
      LOGICAL clonly

      INCLUDE 'parobc.h90'
      INCLUDE 'parrms.h90'

! Common blocks to be initialized:
      INCLUDE 'rmscl.h90'

      INTEGER ic,i,unit,n,obsc,obscp,nin,lf
      INTEGER n1,n2,n3
      DOUBLE PRECISION mjde1,mjde2,rms1,ave1,rms1a
      CHARACTER crmad1*5,crmod1*1,cobsc*3,file*100,rec*120
      LOGICAL mixobs,new1,new2

      INTEGER lench
      DOUBLE PRECISION tjm1
      EXTERNAL lench,tjm1

      DO 10 ic=1,2
      IF(ic.EQ.1) THEN
          file=filra
      ELSE
          file=fildec
      END IF
      lf=lench(file)
      CALL filopl(unit,file)

! Initializations (section 1)
      ncrm(ic)=0
      DO 1 i=obsc1,obsc2
      crmobp(1,i,ic)=0
      crmobp(2,i,ic)=0
    1 END DO
      nin=0
      n=1
      obscp=obsc1-999

! Initializations (section 2)
      n1=0
      n2=0
      n3=0

! Check file format (in order to avoid using files written according
! to old format, not containing bias information)
      IF(.NOT.clonly) THEN
          READ(unit,102,END=3) rec
          IF(lench(rec).LT.110) THEN
              WRITE(*,200) file(1:lf)
              STOP '**** rrmscl: abnormal end ****'
          END IF
          REWIND(unit)
      END IF
  102 FORMAT(A)
  200 FORMAT('**** ERROR: file "',A,'" is written in an obsolete ',     &
     &       'format: please use a more recente version ****')

! Start reading loop
    2 CONTINUE
      IF(clonly) THEN
          READ(unit,100,END=3) crmad1,crmod1,cobsc,mjde1,mjde2
          rms1=0
          ave1=0
          rms1a=0
      ELSE
          READ(unit,100,END=3) crmad1,crmod1,cobsc,mjde1,mjde2,         &
     &                         rms1,ave1,rms1a
      END IF
  100 FORMAT(A5,1X,A1,1X,A3,1X,F8.1,1X,F8.1,1X,F9.3,F11.5,F9.3,         &
     &       E11.3,F9.3,2I8,F7.2,I9)
      nin=nin+1
      IF(crmod1.EQ.' ') crmod1='P'

      mixobs=(cobsc.EQ.'ALL')

      IF(mixobs) THEN
          n3=n3+1
          IF(n3.GT.crx3nx) STOP '**** rrmscl: n3 > crx3nx ****'

! Check for change in accuracy descriptor
          IF(n1.LE.0) THEN
              new1=.true.
          ELSE
              new1=(crmad1.NE.crx1ad(n1,ic))
          END IF
          IF(new1) THEN
              n1=n1+1
              IF(n1.GT.crx1nx) STOP '**** rrmscl: n1 > crx1nx ****'
              n2=n2+1
              IF(n2.GT.crx2nx) STOP '**** rrmscl: n2 > crx2nx ****'
              crx1ad(n1,ic)=crmad1
              crx1pt(1,n1,ic)=n2
              crx1pt(2,n1,ic)=n2
              crx2ty(n2,ic)=crmod1
              crx2pt(1,n2,ic)=n3
          ELSE

! Check for change in observation type
              IF(n2.LE.0) THEN
                  new2=.true.
              ELSE
                  new2=(crmod1.NE.crx2ty(n2,ic))
              END IF
              IF(new2) THEN
                  n2=n2+1
                  IF(n2.GT.crx2nx) STOP '**** rrmscl: n2 > crx2nx ****'
                  crx1pt(2,n1,ic)=n2
                  crx2ty(n2,ic)=crmod1
                  crx2pt(1,n2,ic)=n3
              END IF
          END IF
          crx2pt(2,n2,ic)=n3
          crx3t(1,n3,ic)=mjde1
          crx3t(2,n3,ic)=mjde2
          crx3r(n3,ic)=rms1
          crx3a(n3,ic)=ave1
          crx3ra(n3,ic)=rms1a
      ELSE
          IF(n.GT.ncrmx) STOP '**** rrmscl: ncrm > ncrmx ****'
          CALL statcode(cobsc,obsc)
!         READ(cobsc,101) obsc
          IF(obsc.LT.obsc1 .OR. obsc.GT.obsc2)                          &
     &        STOP '**** rrmscl: input error (01) ****'
          crmad(n,ic)=crmad1
          crmotd(n,ic)=crmod1

! Integer codification of type descriptor
          crmoti(n,ic)=1000+ICHAR(crmotd(n,ic))
          crmot1(n,ic)=mjde1
          crmot2(n,ic)=mjde2
          crmrms(n,ic)=rms1
          crmave(n,ic)=ave1
          crmrma(n,ic)=rms1a

          IF(obsc.NE.obscp) THEN
              IF(crmobp(1,obsc,ic).NE.0)                                &
     &            STOP '**** rrmscl: input error (02) ****'
              crmobp(1,obsc,ic)=n
              obscp=obsc
          END IF
          crmobp(2,obsc,ic)=n
          n=n+1
      END IF
  101 FORMAT(I3)

      GOTO 2

    3 CONTINUE
      ncrm(ic)=n-1
      crx1n(ic)=n1
      crx2n(ic)=n2
      crx3n(ic)=n3
      CALL filclo(unit,' ')
   10 END DO

      iiccrm=36

      END SUBROUTINE rrmscl

! Copyright (C) 1999 by Mario Carpino (carpino@brera.mi.astro.it)
! Version: March 3, 1999
! ---------------------------------------------------------------------
!
!  *****************************************************************
!  *                                                               *
!  *                         A C C S T R                           *
!  *                                                               *
!  *               Accuracy description (string)                   *
!  *                                                               *
!  *****************************************************************
!
! INPUT:    ACCA      -  Accuracy of right ascension (rad)
!           ACCD      -  Accuracy of declination (rad)
!
! OUTPUT:   ADSA      -  Accuracy description string (RA)
!           ADSD      -  Accuracy description string (DEC)
!           ERROR     -  Error flag
!
      SUBROUTINE accstr(acca,accd,adsa,adsd,error)
      USE fund_const
      IMPLICIT NONE

      DOUBLE PRECISION acca,accd
      CHARACTER*(*) adsa,adsd
      LOGICAL error

      DOUBLE PRECISION accs
      INTEGER i,k,lads
      CHARACTER*10 ads

      error=.false.

      DO 1 i=1,2
      IF(i.EQ.1) THEN
          accs=acca*secrad/15.d0
      ELSE
          accs=accd*secrad
      END IF
      WRITE(ads,101) accs
  101 FORMAT(F10.4)
      CALL rmsp(ads,lads)
      DO 2 k=lads,1,-1
      IF(ads(k:k).EQ.'0') THEN
          ads(k:k)=' '
          lads=lads-1
      ELSE
          GOTO 3
      END IF
    2 END DO
    3 CONTINUE
      IF(ads(k:k).EQ.'.') THEN
          ads(k:k)=' '
          lads=lads-1
      END IF
      IF(i.EQ.1) THEN
          IF(lads.GT.LEN(adsa)) THEN
              error=.true.
              WRITE(*,200) 'ADSA',ads(1:lads)
          ELSE
              adsa=ads(1:lads)
          END IF
      ELSE
          IF(lads.GT.LEN(adsd)) THEN
              error=.true.
              WRITE(*,200) 'ADSD',ads(1:lads)
          ELSE
              adsd=ads(1:lads)
          END IF
      END IF
    1 END DO
  200 FORMAT('ERROR (accstr): ',A,' = "',A,'"')

      END SUBROUTINE accstr
! Copyright (C) 1999 by Mario Carpino (carpino@brera.mi.astro.it)
! Version: November 24, 1999
! ---------------------------------------------------------------------
!
!  *****************************************************************
!  *                                                               *
!  *                         C R M S C L                           *
!  *                                                               *
!  *   Computes a-priori observation RMS based on known classes    *
!  *                                                               *
!  *****************************************************************
!
! INPUT:    OBSCOD    -  Observatory code
!           ADS       -  Accuracy description string (RA/DEC)
!           MPCTYP    -  Observation type (column 15 of MPC record)
!           TDT       -  Time (MJD, TDT)
!
! OUTPUT:   RMSA      -  A-priori RMS of RA (arcsec)
!           RMSD      -  A-priori RMS of DEC (arcsec)
!           BIASA     -  Bias in RA (arcsec)
!           BIASD     -  Bias in DEC (arcsec)
!           IDCL      -  Class progressive number in index
!           STEP      -  RMS assignation steps:
!                           0 = not assigned
!                           1 = specific station class
!                           2 = mixed class
!           TDTLIM    -  Class limits (MJD, TDT)
!
! WARNING: if no valid class for the observation is found, the
!          corresponding output rms is set to a negative value
!
      SUBROUTINE crmscl(obscod,ads,mpctyp,tdt,rmsa,rmsd,biasa,biasd,idcl,step,tdtlim)
      IMPLICIT NONE

      INTEGER obscod,idcl(2),step(2)
      DOUBLE PRECISION rmsa,rmsd,biasa,biasd,tdt,tdtlim(2,2)
      CHARACTER*(*) mpctyp,ads(2)

      INCLUDE 'parobc.h90'
      INCLUDE 'parrms.h90'

! NEEDED common blocks:
      INCLUDE 'rmscl.h90'

      INTEGER ic,i,i1,i2,i3
      DOUBLE PRECISION rms1,bias1

      IF(iiccrm.NE.36) STOP '**** crmscl: internal error (01) ****'
      IF(obscod.LT.obsc1 .OR. obscod.GT.obsc2)                          &
     &    STOP '**** crmscl: input error (01) ****'

      DO 20 ic=1,2
      rms1=-1.d0
      bias1=0.d0
      step(ic)=0
      idcl(ic)=0
      tdtlim(1,ic)=-9.D9
      tdtlim(2,ic)=-9.D9
      IF(crmobp(1,obscod,ic).GT.0) THEN
          DO 1 i=crmobp(1,obscod,ic),crmobp(2,obscod,ic)
          IF(crmad(i,ic).NE.ads(ic)) GOTO 1
          IF(crmotd(i,ic).NE.mpctyp) GOTO 1
          IF(tdt.LT.crmot1(i,ic) .OR. tdt.GE.crmot2(i,ic)) GOTO 1
          rms1=crmrms(i,ic)
          bias1=crmave(i,ic)
          idcl(ic)=i
          step(ic)=1
          tdtlim(1,ic)=crmot1(i,ic)
          tdtlim(2,ic)=crmot2(i,ic)
          GOTO 2
    1     CONTINUE
    2     CONTINUE
      END IF
      IF(step(ic).GT.0) GOTO 10
      DO 3 i1=1,crx1n(ic)
      IF(ads(ic).EQ.crx1ad(i1,ic)) THEN
          DO 4 i2=crx1pt(1,i1,ic),crx1pt(2,i1,ic)
          IF(mpctyp.EQ.crx2ty(i2,ic)) THEN
              DO 5 i3=crx2pt(1,i2,ic),crx2pt(2,i2,ic)
              IF(tdt.GE.crx3t(1,i3,ic) .AND. tdt.LT.crx3t(2,i3,ic)) THEN
                  rms1=crx3r(i3,ic)
                  bias1=crx3a(i3,ic)
                  idcl(ic)=i3
                  step(ic)=2
                  tdtlim(1,ic)=crx3t(1,i3,ic)
                  tdtlim(2,ic)=crx3t(2,i3,ic)
                  GOTO 10
              END IF
    5         CONTINUE
          END IF
    4     CONTINUE
      END IF
    3 END DO

   10 CONTINUE
      IF(ic.EQ.1) THEN
          rmsa=rms1
          biasa=bias1
      ELSE
          rmsd=rms1
          biasd=bias1
      END IF
   20 END DO

      END SUBROUTINE crmscl
!
!  *****************************************************************
!  *                                                               *
!  *                         C R M S C N                           *
!  *                                                               *
!  *   Computes a-priori observation RMS based on known classes    *
!  *                                                               *
!  *****************************************************************
!
! INPUT:    OBSCOD    -  Observatory code
!           ADS       -  Accuracy description string (RA/DEC)
!           MPCTYP    -  Observation type (column 15 of MPC record)
!           TDT       -  Time (MJD, TDT)
!
! OUTPUT:   RMSA      -  A-priori RMS of RA (arcsec) before sub. ave.
!           RMSD      -  A-priori RMS of DEC (arcsec) before sub. ave.
!           AVEA      -  Average residual (bias) in RA (arcsec)
!           AVED      -  Average residual (bias) in DEC (arcsec)
!           RMSAA     -  A-priori RMS of RA (arcsec) after sub. ave.
!           RMSDA     -  A-priori RMS of DEC (arcsec) after sub. ave.
!           IDCL      -  Class progressive number in index
!           STEP      -  RMS assignation steps:
!                           0 = not assigned
!                           1 = specific station class
!                           2 = mixed class
!           TDTLIM    -  Class limits (MJD, TDT)
!
! WARNING: if no valid class for the observation is found, the
!          corresponding output rms is set to a negative value
!
      SUBROUTINE crmscn(obscod,ads,mpctyp,tdt,rmsa,rmsd,avea,aved,      &
     &                  rmsaa,rmsda,idcl,step,tdtlim)
      IMPLICIT NONE

      INTEGER obscod,idcl(2),step(2)
      DOUBLE PRECISION rmsa,rmsd,avea,aved,rmsaa,rmsda,tdt,tdtlim(2,2)
      CHARACTER*(*) mpctyp,ads(2)

      INCLUDE 'parobc.h90'
      INCLUDE 'parrms.h90'

! NEEDED common blocks:
      INCLUDE 'rmscl.h90'

      INTEGER ic,i,i1,i2,i3
      DOUBLE PRECISION rms1,ave1,rms1a

      IF(iiccrm.NE.36) STOP '**** crmscn: internal error (01) ****'
      IF(obscod.LT.obsc1 .OR. obscod.GT.obsc2)                          &
     &    STOP '**** crmscn: input error (01) ****'

      DO 20 ic=1,2
      rms1=-1
      ave1=0
      rms1a=-1
      step(ic)=0
      idcl(ic)=0
      tdtlim(1,ic)=-9.D9
      tdtlim(2,ic)=-9.D9
      IF(crmobp(1,obscod,ic).GT.0) THEN
          DO 1 i=crmobp(1,obscod,ic),crmobp(2,obscod,ic)
          IF(crmad(i,ic).NE.ads(ic)) GOTO 1
          IF(crmotd(i,ic).NE.mpctyp) GOTO 1
          IF(tdt.LT.crmot1(i,ic) .OR. tdt.GE.crmot2(i,ic)) GOTO 1
          rms1=crmrms(i,ic)
          ave1=crmave(i,ic)
          rms1a=crmrma(i,ic)
          idcl(ic)=i
          step(ic)=1
          tdtlim(1,ic)=crmot1(i,ic)
          tdtlim(2,ic)=crmot2(i,ic)
          GOTO 2
    1     CONTINUE
    2     CONTINUE
      END IF
      IF(step(ic).GT.0) GOTO 10
      DO 3 i1=1,crx1n(ic)
      IF(ads(ic).EQ.crx1ad(i1,ic)) THEN
          DO 4 i2=crx1pt(1,i1,ic),crx1pt(2,i1,ic)
          IF(mpctyp.EQ.crx2ty(i2,ic)) THEN
              DO 5 i3=crx2pt(1,i2,ic),crx2pt(2,i2,ic)
              IF(tdt.GE.crx3t(1,i3,ic) .AND. tdt.LT.crx3t(2,i3,ic)) THEN
                  rms1=crx3r(i3,ic)
                  ave1=crx3a(i3,ic)
                  rms1a=crx3ra(i3,ic)
                  idcl(ic)=i3
                  step(ic)=2
                  tdtlim(1,ic)=crx3t(1,i3,ic)
                  tdtlim(2,ic)=crx3t(2,i3,ic)
                  GOTO 10
              END IF
    5         CONTINUE
          END IF
    4     CONTINUE
      END IF
    3 END DO

   10 CONTINUE
      IF(ic.EQ.1) THEN
          rmsa=rms1
          avea=ave1
          rmsaa=rms1a
      ELSE
          rmsd=rms1
          aved=ave1
          rmsda=rms1a
      END IF
   20 END DO

      END

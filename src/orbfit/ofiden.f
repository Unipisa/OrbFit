* Copyright (C) 1997-2000 by Mario Carpino (carpino@brera.mi.astro.it)
*                            Zoran Knezevic (zoran@aob.aob.bg.ac.yu)
* Version: June 7, 2000
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         O F I D E N                           *
*  *                                                               *
*  *                Auxiliary routine for ORBFIT:                  *
*  *                     orbit identification                      *
*  *                                                               *
*  *****************************************************************
*
* IN/OUT:   UNIREP    -  FORTRAN unit for report
*           UNIELE    -  FORTRAN unit for output of orbital elements
*           UNIDIF    -  FORTRAN unit for logfilr of DIFCOR
*           OPELE     -  Flag (need to open UNIELE)
*           ELEOUT    -  Name of file for output of orbital elements
*           MAGOPT    -  Magnitude determination option
*           NAME      -  Object names
*           NAMOF     -  Names of observation files (w/o extension)
*           DEFORB    -  Orbit definition flag
*           DEFCN     -  Tells whether covariance/normal matrices
*                            are defined
*           ELEM      -  Orbital elements
*           TELEM     -  Epoch of orbital elements (MJD, TDT)
*           ELTYPE    -  Type of orbital elements (EQU/KEP/CAR)
*           COVE      -  Covariance matrix of orbital elements
*           NORE      -  Normal matrix of orbital elements
*           H,G       -  Magnitude parameters
*           MASS      -  Object masses
*           COMELE    -  Comment on orbital elements
*           DIR       -  Directory containing observation/residual file
*           NOBJ      -  Number of objects
*           OBJID     -  Object identifier
*           TDT       -  Times of observations (MJD, TDT)
*           TUTM      -  Time (MJD, UTM)
*           ALPHA     -  Right ascension (rad)
*           DELTA     -  Declination (rad)
*           IOBS      -  Observation type
*           OBSCOD    -  Observatory code
*           SMAG      -  Measured magnitude (string)
*           RMSA      -  A-priori RMS of RA (rad)
*           RMSD      -  A-priori RMS of DEC (rad)
*           RMSMAG    -  A-priori RMS of magnitude
*           SEL       -  Selection index (0=don't use; 1=use for fit;
*                        2=use for fit & Gauss method)
*           N         -  Number of observations for each object
*           NT        -  Total number of observations
*           IP1       -  Pointer to first observation for each object
*           IP2       -  Pointer to first weight for each object
*           RESA      -  Residuals in RA (rad)
*           RESD      -  Residuals in DEC (rad)
*           RESMAG    -  Residuals in magnitude
*           X2        -  Chi-square residual
*           GMSUN     -  G*Mass(Sun)
*           OETYPE    -  Type of output elements (CAR/EQU/KEP/EQP)
*
      SUBROUTINE ofiden(unirep,uniele,unidif,opele,eleout,magopt,name,
     +                  namof,deforb,defcn,elem,telem,eltype,cove,nore,
     +                  h,g,mass,comele,dir,nobj,objid,tdt,tutm,alpha,
     +                  delta,iobs,obscod,smag,rmsa,rmsd,rmsmag,sel,
     +                  n,nt,ip1,ip2,resa,resd,resmag,x2,gmsun,oetype)
      IMPLICIT NONE

      INCLUDE 'parnob.h'
      INCLUDE 'parobx.h'
      INCLUDE 'parcmc.h'

      INTEGER unirep,uniele,unidif,nobj,n(nobj),nt,magopt
c      INTEGER obscod(nt),ip1(nobj1x),ip2(nobj1x),sel(nt),iobs(nt)
      INTEGER obscod(nt),ip1(nobjx),ip2(nobjx),sel(nt),iobs(nt)
      DOUBLE PRECISION tdt(nt),tutm(nt),gmsun,h(nobj1x),g(nobj1x)
      DOUBLE PRECISION elem(6,nobj1x),telem(nobj1x),mass(nobj1x)
      DOUBLE PRECISION alpha(nt),delta(nt),rmsa(nt),rmsd(nt),rmsmag(nt)
      DOUBLE PRECISION cove(6,6,nobj1x),nore(6,6,nobj1x)
      DOUBLE PRECISION resa(nt),resd(nt),x2(nt),resmag(nt)
      LOGICAL opele,deforb(nobj1x),defcn(nobj1x)
      CHARACTER*(*) eleout,name(nobj1x),namof(nobj1x),eltype(nobj1x)
      CHARACTER*(*) dir(nobj1x),comele(nobj1x),oetype,objid(nt),smag(nt)

      INTEGER i,k,ng,lnt,ln1,icor(6),nsm
      DOUBLE PRECISION telemp(nobjx),elemp(6,nobjx),hnew,rmsh
      DOUBLE PRECISION csinor,delnor,elemo(6),telemo
      DOUBLE PRECISION telemt,elemt(6),masst,ht,gt,ennet,gma1,enne
      DOUBLE PRECISION lsfres(nob2x),w(nob2x),gtwg(6,6)
      CHARACTER namet*100,titnam*80,file*150
      LOGICAL ok,autrep

* NEEDED common blocks:
      INCLUDE 'comidn.h'
      INCLUDE 'comrej.h'

* TUNING PARAMETERS
* Small time (d)
      DOUBLE PRECISION epst
      PARAMETER (epst=1.0D-9)

      INCLUDE 'mag.h'

      INTEGER lench
      EXTERNAL lench

      IF(iicidn.NE.36) STOP '**** ofiden: internal error (01) ****'
      IF(iicrej.NE.36) STOP '**** ofiden: internal error (02) ****'

      IF(nobj.LT.2) RETURN
      DO 1 i=1,nobj
      IF(.NOT.(deforb(i).AND.defcn(i))) RETURN
 1    CONTINUE

      nobj=2
      IF(nobj.GT.nobjx) STOP '**** ofiden: internal error (03) ****'
      IF(nobj1x.LT.3) STOP '**** ofiden: internal error (04) ****'

* Weighting
      CALL fitwgt(rmsa,rmsd,delta,sel,iobs,w,nt,.false.)

      WRITE(unirep,207)
 207  FORMAT('Orbit identification:')

      DO 3 i=1,nobj
      gma1=gmsun*(1+mass(i))
      telemp(i)=telem(i)
      CALL coocha(elem(1,i),eltype(i),gma1,elemp(1,i),'EQU',enne)
      ln1=lench(name(i))
      WRITE(unirep,200) name(i)(1:ln1)
      CALL outele(unirep,elemp(1,i),'EQU',telemp(i),' ',.true.,.false.)
      IF(magopt.NE.0 .AND. h(i).GT.-100.D0) WRITE(unirep,140) h(i)
 3    CONTINUE
 200  FORMAT(5X,'Starting orbital elements for object ',A,':')
 140  FORMAT(10X,'Absolute mag. H    =',F9.2)

* Computation of starting elements (average of single objects)
      telemt=0
      masst=0
      DO 4 k=2,5
      elemt(k)=0
 4    CONTINUE
      ht=0
      gt=0
      nsm=0
      DO 5 i=1,nobj
      telemt=telemt+telemp(i)
      masst=masst+mass(i)
      DO 6 k=2,5
      elemt(k)=elemt(k)+elemp(k,i)
 6    CONTINUE
      IF(h(i).GT.-100.d0) THEN
          ht=ht+h(i)
          gt=gt+g(i)
          nsm=nsm+1
      END IF
 5    CONTINUE
      telemt=telemt/nobj
      masst=masst/nobj
      DO 7 k=2,5
      elemt(k)=elemt(k)/nobj
 7    CONTINUE
      IF(nsm.GT.0) THEN
          ht=ht/nsm
          gt=gt/nsm
      ELSE
          ht=-1.D9
          gt=0.d0
      END IF
      gma1=gmsun*(1+masst)
      CALL start(elemp(1,1),elemp(1,2),telemp(1),telemp(2),gma1,1,
     +           ng,ennet,elemt(1),elemt(6))
      CALL titast(3,namof(1),namof(2),titnam,namet,lnt)
      WRITE(unirep,200) namet(1:lnt)
      CALL outele(unirep,elemt,'EQU',telemt,' ',.true.,.false.)
      IF(magopt.NE.0 .AND. ht.GT.-100.D0) WRITE(unirep,140) ht
      WRITE(unirep,206) namet(1:lnt)
 206  FORMAT('Differential correction for object ',A,':')
      dir(3)=' '
      CALL fitwgt(rmsa,rmsd,delta,sel,iobs,w,nt,.false.)

* Preliminary 2-parameter fit (only semimajor axis and mean anomaly)
      IF(amfit) THEN
          DO 11 i=2,5
          icor(i)=0
 11       CONTINUE
          icor(1)=1
          icor(6)=1
* Disables temporarily outlier rejection
          autrep=autrej
          autrej=.false.
          IF(magopt.NE.0 .AND. ht.GT.-100.D0) gmagc=gt
          CALL difcor(nt,w,sel,telemt,iobs,tdt,obscod,elemt,alpha,
     +                delta,icor,2,unidif,delcr,elem(1,3),cove(1,1,3),
     +                gtwg,csinor,delnor,lsfres,x2,ok)
          autrej=autrep
          IF(ok) THEN
              elemt(1)=elem(1,3)
              elemt(6)=elem(6,3)
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
 10   CONTINUE
      IF(magopt.NE.0 .AND. ht.GT.-100.D0) gmagc=gt
      CALL difcor(nt,w,sel,telemt,iobs,tdt,obscod,elemt,alpha,
     +            delta,icor,2,unidif,delcr,elem(1,3),cove(1,1,3),
     +            nore(1,1,3),csinor,delnor,lsfres,x2,ok)
      IF(ok) THEN
          deforb(3)=.true.
          defcn(3)=.true.
          eltype(3)='EQU'
          comele(3)='computed with least squares fit'
          telem(3)=telemt
* Magnitude estimation
          IF(magopt.NE.0) THEN
              CALL magest(smag,rmsmag,sel,nt,hnew,resmag,rmsh)
              IF(rmsh.GT.0.d0) ht=hnew
          END IF
          h(3)=ht
          g(3)=gt
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
          CALL coocha(elem(1,3),eltype(3),gma1,elemo,oetype,enne)
          telemo=telem(3)
          WRITE(unirep,124)
          CALL outele(unirep,elemo,oetype,telemo,' ',.true.,.false.)
          IF(magopt.NE.0 .AND. h(3).GT.-100.D0) WRITE(unirep,140) h(3)
          WRITE(unirep,230) csinor,delnor
          CALL wromlr(uniele,name(3),elem(1,3),eltype(3),telem(3),
     +                cove(1,1,3),defcn(3),nore(1,1,3),defcn(3),
     +                h(3),g(3),mass(3))
          DO 8 k=1,nt
          resa(k)=lsfres(2*k-1)
          resd(k)=lsfres(2*k)
 8        CONTINUE
          file=namet(1:lnt)//'.rwo'
          CALL wrirwo(file,objid,iobs,tutm,obscod,alpha,rmsa,resa,
     +                delta,rmsd,resd,smag,rmsmag,resmag,rmsh,
     +                sel,x2,nt,csinor)
      ELSE
          WRITE(unirep,202)
      END IF
 202  FORMAT(5X,'FAILED')
 301  FORMAT(A,' Orbits computed with differential corrections')
 124  FORMAT(5X,'Corrected orbital elements:')
 125  FORMAT(5X,'Intermediate orbital elements (a-M correction):')
 230  FORMAT(5X,'Residual norm   =',1P,E11.3/
     +       5X,'Correction norm =',1P,E11.3)

      END

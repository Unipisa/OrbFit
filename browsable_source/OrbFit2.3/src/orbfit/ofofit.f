* Copyright (C) 1997-2000 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: June 7, 2000
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         O F O F I T                           *
*  *                                                               *
*  *                Auxiliary routine for ORBFIT:                  *
*  *                  least squares orbital fit                    *
*  *                                                               *
*  *****************************************************************
*
* IN/OUT:   UNIREP    -  FORTRAN unit for report
*           UNIELE    -  FORTRAN unit for output of orbital elements
*           UNIDIF    -  FORTRAN unit for logfilr of DIFCOR
*           OPELE     -  Flag (need to open UNIELE)
*           ELEOUT    -  Name of file for output of orbital elements
*           ORBOPT    -  Orbit determination option
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
*           OEPTIM    -  Epoch of output elements (MJD, TDT)
*           OEPSET    -  Flag stating that an output epoch is requested
*           OETYPE    -  Type of output elements (CAR/EQU/KEP/EQP)
*
      SUBROUTINE ofofit(unirep,uniele,unidif,opele,eleout,orbopt,magopt,
     +                  name,namof,deforb,defcn,elem,telem,eltype,cove,
     +                  nore,h,g,mass,comele,dir,nobj,objid,tdt,tutm,
     +                  alpha,delta,iobs,obscod,smag,rmsa,rmsd,rmsmag,
     +                  sel,n,nt,ip1,ip2,resa,resd,resmag,x2,
     +                  gmsun,oeptim,oepset,oetype)
      IMPLICIT NONE

      INTEGER unirep,uniele,unidif,orbopt,magopt,nobj,n(nobj),nt
      INTEGER ip1(nobj),ip2(nobj),sel(nt),iobs(nt),obscod(nt)
      DOUBLE PRECISION mass(nobj),gmsun,elem(6,nobj),telem(nobj)
      DOUBLE PRECISION tdt(nt),tutm(nt),alpha(nt),delta(nt)
      DOUBLE PRECISION cove(6,6,nobj),nore(6,6,nobj)
      DOUBLE PRECISION h(nobj),g(nobj),oeptim,x2(nt)
      DOUBLE PRECISION rmsa(nt),rmsd(nt),rmsmag(nt)
      DOUBLE PRECISION resa(nt),resd(nt),resmag(nt)
      LOGICAL deforb(nobj),defcn(nobj),opele,oepset
      CHARACTER*(*) eltype(nobj),name(nobj),namof(nobj),eleout,dir(nobj)
      CHARACTER*(*) comele(nobj),oetype,objid(nt),smag(nt)

      INCLUDE 'parobx.h'
      INCLUDE 'parcmc.h'
      INCLUDE 'trig.h'

* NEEDED common blocks:
      INCLUDE 'comlsf.h'

* TUNING PARAMETERS
* Max RMS of a "good" observation (arcsec)
      DOUBLE PRECISION maxrms
      PARAMETER (maxrms=5.d0)
* Small time (d)
      DOUBLE PRECISION epst
      PARAMETER (epst=1.0D-9)

      INTEGER i,ln,icor(6),j1,j2,ld
      DOUBLE PRECISION gma1,elemt(6),enne,csinor,delnor,rms1,t1,t2
      DOUBLE PRECISION elemn(6),coven(6,6),gtwg(6,6),xea(6),telemp
      DOUBLE PRECISION lsfres(nob2x),w(nob2x),newep,elemp(6),rmsh
      DOUBLE PRECISION dxde(6,6),ddxde(3,6,6),covep(6,6),norep(6,6)
      DOUBLE PRECISION hnew
      CHARACTER file*150
      LOGICAL doit,ok,chep,error,defcov,defnor

      INCLUDE 'mag.h'

      INTEGER lench
      EXTERNAL lench

      IF(iiclsf.NE.36) STOP '**** ofofit: internal error (01) ****'
      rms1=maxrms*radsec

      DO 10 i=1,6
      icor(i)=1
 10   CONTINUE

* Weighting
      CALL fitwgt(rmsa,rmsd,delta,sel,iobs,w,nt,.false.)

      DO 1 i=1,nobj
      IF(.NOT.deforb(i)) GOTO 1
* This check is not enough (many things could have changed since the last
* orbit determination)
***   IF(orbopt.EQ.1) THEN
***       doit=(.NOT.defcn(i))
***   ELSE
          doit=.true.
***   END IF
      IF(.NOT.doit) GOTO 1
      gma1=gmsun*(1+mass(i))

* Selection of a suitable time of initial conditions for the fit
* Useful timespan of observations (discarding observation of low
* accuracy)
      CALL ustsp(tdt(ip1(i)),rmsa(ip1(i)),rmsd(ip1(i)),sel(ip1(i)),
     +           n(i),rms1,t1,t2)
      chep=.false.
      newep=telem(i)
* Use the time requested for output of orbital elements, but only
* if it is within the useful timespan of observations
      IF(oepset) THEN
          IF(oeptim.GE.t1 .AND. oeptim.LE.t2) THEN
              chep=.true.
              newep=oeptim
          END IF
      END IF
* Otherwise, use the end-point of the observation timespan nearest
* to the requested time
      IF(.NOT.chep) THEN
          IF(telem(i).LT.t1) THEN
              chep=.true.
              newep=t1
          ELSEIF(telem(i).GT.t2) THEN
              chep=.true.
              newep=t2
          END IF
      END IF

* Propagation to new epoch
      IF(chep .AND. ABS(newep-telem(i)).GT.epst) THEN
          CALL coocha(elem(1,i),eltype(i),gma1,elemt,'EQU',enne)
          CALL propag(telem(i),elemt,newep,elemn,xea,0,dxde,ddxde)
          CALL coocha(elemn,'CAR',gma1,elemt,'EQU',enne)
          telem(i)=newep
      ELSE
          CALL coocha(elem(1,i),eltype(i),gma1,elemt,'EQU',enne)
      END IF

      ln=lench(name(i))
      WRITE(unirep,206) name(i)(1:ln)
 206  FORMAT('Differential correction for object ',A,':')
      WRITE(unirep,123)
 123  FORMAT(5X,'Starting orbital elements:')
      CALL outele(unirep,elemt,'EQU',telem(i),' ',.true.,.false.)
      IF(magopt.NE.0) THEN
          gmagc=g(i)
          IF(h(i).GT.-100.D0) WRITE(unirep,140) h(i)
      END IF
 140  FORMAT(10X,'Absolute mag. H    =',F9.2)
      CALL difcor(n(i),w(ip2(i)),sel(i),telem(i),iobs(ip1(i)),
     +            tdt(ip1(i)),obscod(ip1(i)),elemt,alpha(ip1(i)),
     +            delta(ip1(i)),icor,2,unidif,delcr,elemn,coven,gtwg,
     +            csinor,delnor,lsfres,x2(ip1(i)),ok)
      IF(ok) THEN
* Magnitude estimation
          IF(magopt.NE.0) THEN
              CALL magest(smag(ip1(i)),rmsmag(ip1(i)),sel(ip1(i)),n(i),
     +                    hnew,resmag(ip1(i)),rmsh)
              IF(rmsh.GT.0.d0) h(i)=hnew
          END IF
          DO 2 j1=1,6
          elem(j1,i)=elemn(j1)
          DO 3 j2=1,6
          cove(j1,j2,i)=coven(j1,j2)
          nore(j1,j2,i)=gtwg(j1,j2)
 3        CONTINUE
 2        CONTINUE
          defcn(i)=.true.
          eltype(i)='EQU'
          comele(i)='computed with least squares fit'
          IF(opele) THEN
              OPEN(uniele,FILE=eleout,STATUS='UNKNOWN')
              opele=.false.
              WRITE(uniele,301) comcha
              CALL wromlh(uniele,'ECLM','J2000')
          END IF
          CALL cooder(elem(1,i),eltype(i),gma1,elemp,oetype,enne,dxde)
          CALL covprs(cove(1,1,i),dxde,6,covep)
          defcov=.true.
          CALL norprs(nore(1,1,i),dxde,6,norep,error)
          defnor=(.NOT.error)
          telemp=telem(i)
          WRITE(unirep,124)
          CALL outele(unirep,elemp,oetype,telem(i),' ',.true.,.false.)
          IF(magopt.NE.0 .AND. h(i).GT.-100.D0) WRITE(unirep,140) h(i)
          WRITE(unirep,230) csinor,delnor
          DO 4 j1=ip1(i),ip1(i)+n(i)-1
          resa(j1)=lsfres(2*j1-1)
          resd(j1)=lsfres(2*j1)
 4        CONTINUE
          CALL wromlr(uniele,name(i),elemp,oetype,telemp,
     +                covep,defcov,norep,defnor,h(i),g(i),mass(i))
          ld=lench(dir(i))
          ln=lench(namof(i))
          IF(ld.GT.0) THEN
              file=dir(i)(1:ld)//namof(i)(1:ln)//'.rwo'
          ELSE
              file=namof(i)(1:ln)//'.rwo'
          END IF
          rmsh=0
          CALL wrirwo(file,objid(ip1(i)),iobs(ip1(i)),tutm(ip1(i)),
     +                obscod(ip1(i)),alpha(ip1(i)),rmsa(ip1(i)),
     +                resa(ip1(i)),delta(ip1(i)),rmsd(ip1(i)),
     +                resd(ip1(i)),smag(ip1(i)),rmsmag(ip1(i)),
     +                resmag(ip1(i)),rmsh,sel(ip1(i)),x2(ip1(i)),
     +                n(i),csinor)
      ELSE
          WRITE(unirep,202)
      END IF

 1    CONTINUE
 124  FORMAT(5X,'Corrected orbital elements:')
 301  FORMAT(A,' Orbits computed with differential corrections')
 202  FORMAT(5X,'FAILED')
 230  FORMAT(5X,'Residual norm   =',1P,E11.3/
     +       5X,'Correction norm =',1P,E11.3)

      END

* Copyright (C) 1997-1999 Mario Carpino (carpino@brera.mi.astro.it)
* Version: May 31, 1999
* ---------------------------------------------------------------------
*
* Orbit determination
*
      PROGRAM orbfit
      IMPLICIT NONE

      INCLUDE 'parobx.h'

* =========================== DECLARATIONS ===========================
* Fixed FORTRAN I/O units
      INTEGER unirep,uniele,unidif,unieph
* File names
      CHARACTER*80 run,file,optfil,eleout

      INCLUDE 'parnob.h'

* Object names and input directories
      CHARACTER*80 name(nobj1x),nameo(nobj1x),dir(nobj1x)
      CHARACTER*80 namof(nobj1x)

* OBSERVATIONS
      INTEGER nt,n(nobjx),ip1(nobjx),ip2(nobjx)
      DOUBLE PRECISION tdt(nobx),tutm(nobx),alpha(nobx),delta(nobx)
      DOUBLE PRECISION x2(nobx)
      INTEGER obscod(nobx),sel(nobx),iobs(nobx)
      CHARACTER objid(nobx)*9,smag(nobx)*6
* A-priori RMS of observations (defined from accuracy OR read from file)
      DOUBLE PRECISION rmsa(nobx),rmsd(nobx),rmsmag(nobx)
* Fit residuals
      DOUBLE PRECISION resa(nobx),resd(nobx),resmag(nobx)

* Max number of orbital element files (for each object)
      INTEGER nifx
      PARAMETER (nifx=10)

* ORBITS
      INTEGER nelft,nelf1(nobjx)
      DOUBLE PRECISION telem(nobj1x),elem(6,nobj1x),mass(nobj1x)
      DOUBLE PRECISION cove(6,6,nobj1x),nore(6,6,nobj1x)
      DOUBLE PRECISION h(nobj1x),g(nobj1x)
      LOGICAL deforb(nobj1x),defcn(nobj1x)
      CHARACTER*10 eltype(nobj1x),oetype
      CHARACTER*120 comele(nobj1x),elft(nifx),elf1(nifx,nobjx)

      DOUBLE PRECISION gmsun,oeptim
      INTEGER lr,unit,nobj,op(10),lf,nfound
      LOGICAL found,opdif,opele,oepset,needrm
c modification to select preliminary orbit method

      INTEGER lench
      EXTERNAL lench

c verbosity control

      INCLUDE 'verbosity.h'
c verbosity levels for an interactive program, but some a bit less
      verb_pro=10
      verb_clo=10
      verb_dif=10
      verb_mul=5
      verb_rej=5

* Fixed FORTRAN I/O units:
* report file
      unirep=1
* output orbital elements
      uniele=2
* log of difcor
      unidif=3
* ephemerides
      unieph=7

* Flags indicating whether output orbital element file and .odc file
* are to be opened or not
      opdif=.true.
      opele=.true.

      WRITE(*,110)
 110  FORMAT(' Run name =')
      READ(*,100) run
 100  FORMAT(A)
      lr=lench(run)
      WRITE(*,111) run(1:lr)
 111  FORMAT(' Run name = ',A)

      file=run(1:lr)//'.olg'
      OPEN(unirep,FILE=file,STATUS='UNKNOWN')

* ========================== INITIALIZATION ==========================
      nobj=0
      gmsun=0.01720209895d0**2
      CALL namini
      CALL libini

* ========================= INPUT OF OPTIONS =========================
* Read option files
      CALL filopf(unit,'orbfit.def',found)
      IF(found) THEN
          CALL rdnam(unit)
          CALL filclo(unit,' ')
      ELSE
          WRITE(*,210)
      END IF
 210  FORMAT(' WARNING: no default option file (orbfit.def) found')
      optfil=run(1:lr)//'.oop'
      INQUIRE(FILE=optfil,EXIST=found)
      IF(found) THEN
          CALL filopn(unit,optfil,'OLD')
          CALL rdnam(unit)
          CALL filclo(unit,' ')
      ELSE
          lf=lench(optfil)
          WRITE(*,211) optfil(1:lf)
      END IF
 211  FORMAT(' WARNING: no run-specific option file (',A,') found')

* Check of keywords
      CALL rdklst('orbfit.key')
      CALL chkkey

* Input of options
      CALL rdopto(run,op,name,nameo,namof,dir,nobj,elft,nelft,
     +            elf1,nelf1,eleout,oeptim,oepset,oetype,
     +            nifx,'nifx')
      IF(nobj.GT.nobjx) STOP '**** orbfit: nobj > nobjx ****'
* Check list of perturbing asteroids
      IF(nameo(2).ne.' ') THEN
          CALL selpert2(nameo(1),nameo(2),nfound)
      ELSE
          CALL selpert(nameo(1),found)
      END IF
* Additional options (not always required)
      needrm=.false.
* Options for initial orbit determination
      IF(op(2).GT.0) CALL iodini
* Options for differential correction
      IF(op(3).GT.0) THEN
          CALL rdoptf
          needrm=.true.
      END IF
* Options for identifications
      IF(op(4).GT.0) THEN
          CALL rdopti
          needrm=.true.
      END IF
* Options for ephemerides
      IF(op(5).GT.0) THEN
          CALL rdopte
          needrm=.true.
      END IF
      IF(oepset) needrm=.true.
      IF(needrm) THEN
          CALL rmodel(run)
          CALL ofinip(run)
      END IF

* ====================== INPUT OF OBSERVATIONS =======================
      CALL ofiobs(unirep,name,namof,dir,nobj,nobx,op(1),n,nt,
     +            ip1,ip2,objid,tdt,tutm,alpha,delta,iobs,smag,
     +            obscod,sel,rmsa,rmsd,rmsmag)

* ==================== INPUT OF ORBITAL ELEMENTS =====================
      CALL ofiorb(unirep,elft,nelft,elf1,nelf1,name,nameo,
     +            nobj,nobj1x,deforb,defcn,eltype,elem,
     +            telem,cove,nore,h,g,mass,comele,nifx)

* =================== INITIAL ORBIT DETERMINATION ====================
      IF(op(2).GT.0) THEN
          CALL ofinod(unirep,eleout,op(2),name,namof,deforb,defcn,h,g,
     +                dir,nobj,objid,tdt,tutm,alpha,delta,iobs,
     +                obscod,smag,rmsa,rmsd,rmsmag,sel,n,nt,ip1,elem,
     +                telem,eltype,comele,resa,resd)
      END IF

* =========== LEAST SQUARES ORBITAL FIT (SEPARATE ORBITS) ============
      IF(op(3).GT.0) THEN
          IF(opdif) THEN
              file=run(1:lr)//'.odc'
              OPEN(unidif,FILE=file,STATUS='UNKNOWN')
              opdif=.false.
          END IF
          CALL ofofit(unirep,uniele,unidif,opele,eleout,op(3),op(6),
     +                name,namof,deforb,defcn,elem,telem,eltype,cove,
     +                nore,h,g,mass,comele,dir,nobj,objid,tdt,tutm,
     +                alpha,delta,iobs,obscod,smag,rmsa,rmsd,rmsmag,
     +                sel,n,nt,ip1,ip2,resa,resd,resmag,x2,
     +                gmsun,oeptim,oepset,oetype)
      END IF

* ======================= ORBIT IDENTIFICATION =======================
      IF(op(4).GT.0) THEN
          IF(opdif) THEN
              file=run(1:lr)//'.odc'
              OPEN(unidif,FILE=file,STATUS='UNKNOWN')
              opdif=.false.
          END IF
          CALL ofiden(unirep,uniele,unidif,opele,eleout,op(6),name,
     +                namof,deforb,defcn,elem,telem,eltype,cove,nore,
     +                h,g,mass,comele,dir,nobj,objid,tdt,tutm,alpha,
     +                delta,iobs,obscod,smag,rmsa,rmsd,rmsmag,sel,
     +                n,nt,ip1,ip2,resa,resd,resmag,x2,gmsun,oetype)
      END IF

* ================= PROPAGATION OF ORBITAL ELEMENTS ==================
      IF(oepset) THEN
          IF(.NOT.opele) CLOSE(uniele)
          opele=.true.
          CALL ofprop(unirep,uniele,opele,eleout,name,deforb,defcn,
     +                elem,telem,eltype,cove,nore,h,g,mass,
     +                comele,nobj,gmsun,oeptim,oepset,oetype)
      END IF

* ============================ EPHEMERIDES ===========================
      IF(op(5).GT.0) THEN
          file=run(1:lr)//'.oep'
          OPEN(unieph,FILE=file,STATUS='UNKNOWN')
          CALL ofephe(unieph,name,deforb,defcn,elem,telem,eltype,cove,
     +                h,g,mass,comele,nobj)
      END IF

      CALL ofclrf

      END

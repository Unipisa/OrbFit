* Copyright (C) 1998-2000 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: June 9, 2000
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                          I O D E T                            *
*  *                                                               *
*  *      Initial orbit determination (Gauss, Vaisala, etc.)       *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    UNIREP    -  FORTRAN unit for report (0 for no report)
*           RWOFIL    -  File for residual output (or blank)
*           NAME      -  Object name
*           OBJID     -  Object identifier
*           TDT       -  Times of observations (MJD, TDT)
*           TUTM      -  Time (MJD, UTM)
*           ALPHA     -  Right ascension (rad)
*           DELTA     -  Declination (rad)
*           IOBS      -  Observation type (1=RA,DEC; 2=R,RDOT)
*           OBSCOD    -  Observatory code
*           SMAG      -  Measured magnitude (string)
*           RMSA      -  A-priori RMS of RA (rad)
*           RMSD      -  A-priori RMS of DEC (rad)
*           RMSMAG    -  A-priori RMS of magnitude
*           SEL       -  Selection index (0=don't use; 1=use for fit;
*                                         2=use for fit & Gauss method)
*           N         -  Number of observations
*           IFO       -  Vaisala interactive level:
*                           1=interactive (fitobs)
*                           2=noninteractive (orbfit)
*
* OUTPUT:   ELEM      -  Orbital elements
*           TELEM     -  Epoch of orbital elements (MJD, TDT)
*           RESA      -  Residuals in RA (rad)
*           RESD      -  Residuals in DEC (rad)
*           ELTYPE    -  Type of orbital elements (EQU/KEP/CAR)
*           COMELE    -  Comment on orbital elements
*           FAIL      -  Error flag
*
* CONCEPT: the routines tries to compute a preliminary orbit using
* different triplets of observations and different method, possibly
* applying also small random increments to the astrometric measurements
* in order to counteract observation errors, but returns only one
* solution (the first elliptic solution encountered)
*
      SUBROUTINE iodet(unirep,rwofil,name,objid,tdt,tutm,alpha,delta,
     +                 iobs,obscod,smag,rmsa,rmsd,rmsmag,sel,n,ifo,elem,
     +                 telem,resa,resd,eltype,comele,fail)
      IMPLICIT NONE

      INTEGER unirep,n,ifo
      INTEGER sel(n),obscod(n),iobs(n)
      DOUBLE PRECISION tdt(n),rmsa(n),rmsd(n),rmsmag(n)
      DOUBLE PRECISION alpha(n),delta(n)
      DOUBLE PRECISION elem(6),telem,resa(n),resd(n),tutm(n)
      CHARACTER*(*) name,eltype,comele,rwofil,objid(n),smag(n)
      LOGICAL fail

      INCLUDE 'pariod.h'
      INCLUDE 'parobx.h'
      INCLUDE 'trig.h'

* NEEDED common blocks:
      INCLUDE 'comiod.h'

      INTEGER ncx
      PARAMETER (ncx=7)

      DOUBLE PRECISION elemv(6,8),t0v(8),tdt3(3),alpha3(3),delta3(3)
      DOUBLE PRECISION rmsb,rms1,noisea(3),noised(3),dt1,dt2,dtc
      INTEGER ln,lm,selipt(3),itry,nsol,i,lmsg,ngr,is1,is2,is3,is
      INTEGER ir(2),irb(2),obsc3(3),nkep,in,selb(3)
      INTEGER nc(ncx)
      CHARACTER msg*100,methb*20,label*50
      CHARACTER*4 ieltyv(8),ieltyb
      LOGICAL eot,existb

      INTEGER lench
      DOUBLE PRECISION rvnorm
      EXTERNAL lench,rvnorm

      IF(iiciod.NE.36) STOP '**** iodet: internal error (01) ****'

      fail=.true.
      existb=.false.
      rmsb=0
      methb=' '

      ln=lench(name)
      WRITE(*,202) name(1:ln)
      IF(unirep.GT.0) WRITE(unirep,202) name(1:ln)
 202  FORMAT('Preliminary orbit for object ',A,':')

* Computation of best triplets of observations
      CALL sel3mc(tdt,iobs,rmsa,rmsd,sel,n,iodntr)

* Statistics:
* NC(1) = total number of triplets tried
* NC(2) = total number of noise sample added (including no noise)
* NC(3) = total number of roots (Gauss' method)
* NC(4) = total number of acceptable solutions (Gauss' method)
* NC(5) = total number of elliptic solutions (Gauss' method)
* NC(6) = total number of trials (Vaisala method)
* NC(7) = total number of acceptable solutions (Vaisala method)
      DO 7 i=1,ncx
      nc(i)=0
 7    CONTINUE

* LOOP 1 (on triplets of observations)
 10   CONTINUE
      CALL sel3mg(selipt,eot)
      IF(eot) GOTO 20
* Total number of triplets tried
      nc(1)=nc(1)+1

* Selected points
      is1=selipt(1)
      is2=selipt(2)
      is3=selipt(3)
      DO 4 i=1,3
      is=selipt(i)
      tdt3(i)=tdt(is)
      obsc3(i)=obscod(is)
 4    CONTINUE

      IF(iodvrb.GE.2) THEN
          WRITE(*,203) 'Trying',selipt,
     +                 tdt(is2)-tdt(is1),tdt(is3)-tdt(is2)
          IF(unirep.GT.0) WRITE(unirep,203) 'Trying',selipt,
     +                    tdt(is2)-tdt(is1),tdt(is3)-tdt(is2)
      END IF
 203  FORMAT(4X,A,' observations:',3I6,' (DT=',2F10.2,' d)')
      CALL iodsdt(selipt,tdt,n,iodexp,ioddtm,ir)
      IF(iodvrb.GE.2) THEN
          WRITE(*,204) ir,tdt(ir(2))-tdt(ir(1))
          IF(unirep.GT.0) WRITE(unirep,204) ir,tdt(ir(2))-tdt(ir(1))
      END IF
 204  FORMAT(4X,'RMS check observations:', I6,' -',I6,
     +          ' (DT=',F10.2,' d)')

* LOOP 2 (on different realizations of noise)
      DO 11 in=0,iodnit
      IF(in.EQ.0) THEN
          DO 14 i=1,3
          is=selipt(i)
          noisea(i)=0
          noised(i)=0
          alpha3(i)=alpha(is)
          delta3(i)=delta(is)
 14       CONTINUE
      ELSE
          DO 15 i=1,3
          is=selipt(i)
          noisea(i)=rmsa(is)*iodksi*rvnorm()
          noised(i)=rmsd(is)*iodksi*rvnorm()
          alpha3(i)=alpha(is)+noisea(i)
          delta3(i)=delta(is)+noised(i)
 15       CONTINUE
      END IF
      nc(2)=nc(2)+1

* LOOP 3 (on different methods)
      DO 1 itry=1,iodnm
      lm=lench(iodmen(itry))
      msg=' '
      IF((iodvrb.GE.2) .AND. iodmul) THEN
          WRITE(*,210) iodmen(itry)(1:lm)
          IF(unirep.GT.0) WRITE(unirep,210) iodmen(itry)(1:lm)
      END IF
 210  FORMAT(8X,'Trying ',A,' method')
      IF(iodvrb.GE.3) THEN
          WRITE(*,240) iodmen(itry)(1:lm),name(1:ln),in,
     +                 (noisea(i)*secrad,noised(i)*secrad,i=1,3)
          IF(unirep.GT.0) WRITE(unirep,240) iodmen(itry)(1:lm),
     +                    name(1:ln),in,
     +                    (noisea(i)*secrad,noised(i)*secrad,i=1,3)
      END IF
      IF(iodmet(itry).EQ.1) THEN
          CALL gaussn(tdt3,alpha3,delta3,obsc3,elemv,ieltyv,t0v,
     +                ngr,nsol,fail,msg,(iodvrb.GE.3),iodmul)
          nc(3)=nc(3)+ngr
          nc(4)=nc(4)+nsol
          nkep=0
          DO 2 is=1,nsol
          IF(ieltyv(is).EQ.'KEP') nkep=nkep+1
 2        CONTINUE
          nc(5)=nc(5)+nkep

          IF(iodvrb.GE.2) THEN
              lmsg=lench(msg)
              IF(lmsg.GT.0) THEN
                  WRITE(*,220) 'Gauss',name(1:ln),ngr,nsol,nkep,
     +                         msg(1:lmsg)
                  IF(unirep.GT.0) WRITE(unirep,220) 'Gauss',name(1:ln),
     +                            ngr,nsol,nkep,msg(1:lmsg)
              ELSE
                  WRITE(*,221) 'Gauss',name(1:ln),ngr,nsol,nkep
                  IF(unirep.GT.0) WRITE(unirep,221)'Gauss', name(1:ln),
     +                            ngr,nsol,nkep
              END IF
          END IF
          DO 3 is=1,nsol
          CALL iodrms(elemv(1,is),ieltyv(is),t0v(is),tdt,alpha,delta,
     +                iobs,obscod,sel,ir(1),ir(2),resa,resd,rms1)
          IF(iodvrb.GE.2) THEN
              IF(ieltyv(is).EQ.'KEP') THEN
                  WRITE(*,222) 'Gauss',name(1:ln),is,'a',elemv(1,is),
     +                          elemv(2,is),rms1*secrad
                  IF(unirep.GT.0) WRITE(unirep,222)'Gauss', name(1:ln),
     +                    is,'a',elemv(1,is),elemv(2,is),rms1*secrad
              ELSE
                  WRITE(*,222) 'Gauss',name(1:ln),is,'q',elemv(1,is),
     +                          elemv(2,is),rms1*secrad
                  IF(unirep.GT.0) WRITE(unirep,222) 'Gauss',name(1:ln),
     +                     is,'q',elemv(1,is),elemv(2,is),rms1*secrad
              END IF
          END IF
          IF((ieltyv(is).EQ.'KEP') .AND. (rms1.LE.iodrmx))
     +        CALL iodsbs(elem,ieltyb,telem,rmsb,methb,selb,irb,
     +                    elemv(1,is),ieltyv(is),t0v(is),rms1,'Gauss',
     +                    selipt,ir,existb)
 3        CONTINUE
          IF(existb .AND. (rmsb.LE.iodrok)) GOTO 20

      ELSEIF(iodmet(itry).EQ.2) THEN
          nc(6)=nc(6)+1
          CALL vaisala(tdt3,alpha3,delta3,obsc3,ifo,
     +         elemv(1,1),t0v(1),fail)
          nsol=1
          IF(fail) nsol=0
          nc(7)=nc(7)+nsol
          IF(iodvrb.GE.2) THEN
              WRITE(*,221) 'Vaisala',name(1:ln),nsol,nsol,nsol
              IF(unirep.GT.0) WRITE(unirep,221) 'Vaisala',name(1:ln),
     +                                           nsol,nsol,nsol
          END IF
          ieltyv(1)='KEP'
          DO 5 is=1,nsol
          CALL iodrms(elemv(1,is),ieltyv(is),t0v(is),tdt,alpha,delta,
     +                iobs,obscod,sel,ir(1),ir(2),resa,resd,rms1)
          IF(iodvrb.GE.2) THEN
              IF(ieltyv(is).EQ.'KEP') THEN
                  WRITE(*,222) 'Vaisala',name(1:ln),is,'a',elemv(1,is),
     +                         elemv(2,is),rms1*secrad
                  IF(unirep.GT.0) WRITE(unirep,222) 'Vaisala',
     +                            name(1:ln),is,'a',elemv(1,is),
     +                            elemv(2,is),rms1*secrad
              ELSE
                  WRITE(*,222) 'Vaisala',name(1:ln),is,'q',elemv(1,is),
     +                          elemv(2,is),rms1*secrad
                  IF(unirep.GT.0) WRITE(unirep,222) 'Vaisala',
     +                            name(1:ln),is,'q',elemv(1,is),
     +                            elemv(2,is),rms1*secrad
              END IF
          END IF
          IF((ieltyv(is).EQ.'KEP') .AND. (rms1.LE.iodrmx))
     +        CALL iodsbs(elem,ieltyb,telem,rmsb,methb,selb,irb,
     +                    elemv(1,is),ieltyv(is),t0v(is),rms1,'Vaisala',
     +                    selipt,ir,existb)
 5        CONTINUE
          IF(existb .AND. (rmsb.LE.iodrok)) GOTO 20
 6        CONTINUE
      END IF
 220  FORMAT(8X,A,'(',A,'): NR=',I1,' NC=',I1,' NK=',I1,' (',A,')')
 221  FORMAT(8X,A,'(',A,'): NR=',I1,' NC=',I1,' NK=',I1)
 222  FORMAT(12X,A,'(',A,'#',I1,'): ',A,'=',F9.4,' e=',F9.4,' RMS=',
     +       1P,E10.2,' arcsec')
 240  FORMAT(8X,A,'(',A,'): Iter =',I4,'; Noise =',6F9.3)

* END OF LOOP 3
 1    CONTINUE

* END OF LOOP 1
 11   CONTINUE

* END OF LOOP 1
      GOTO 10
 20   CONTINUE
      IF(existb) THEN
          dt1=tdt(selb(2))-tdt(selb(1))
          dt2=tdt(selb(3))-tdt(selb(2))
          dtc=tdt(irb(2))-tdt(irb(1))
          WRITE(*,230) name(1:ln),nc,selb,irb,dt1,dt2,dtc
          IF(unirep.GT.0) WRITE(unirep,230)
     +                        name(1:ln),nc,selb,ir,dt1,dt2,dtc
      ELSE
          WRITE(*,235) name(1:ln),nc
          IF(unirep.GT.0) WRITE(unirep,235) name(1:ln),nc
      END IF
 230  FORMAT(4X,'IOD(',A,') STAT:',7I6,' SEL:',3I6,' IR:',2I6,' DT:',
     +       3F10.2)
 235  FORMAT(4X,'IOD(',A,') STAT:',7I6)

      fail=(.NOT.existb)
      IF(existb) THEN
          lm=lench(methb)
          IF(in.LE.0) THEN
              WRITE(*,232) name(1:ln),methb(1:lm),rmsb*secrad
              IF(unirep.GT.0) WRITE(unirep,232) name(1:ln),methb(1:lm),
     +                                          rmsb*secrad
              comele=methb(1:lm)
          ELSE
              WRITE(*,233) name(1:ln),methb(1:lm),rmsb*secrad
              IF(unirep.GT.0) WRITE(unirep,233) name(1:ln),methb(1:lm),
     +                                          rmsb*secrad
              comele=methb(1:lm)//'(wN)'
          END IF
          label=name(1:ln)//'/'//methb(1:lm)
          CALL outele(unirep,elem,ieltyb,telem,label,iodmul,.true.)
          IF(ieltyb.NE.'KEP')
     +        STOP '**** iodet: internal error (02) ****'
          eltype='KEP'
          CALL iodrms(elem,ieltyb,telem,tdt,alpha,delta,
     +                iobs,obscod,sel,1,n,resa,resd,rms1)
          sel(selb(1))=2
          sel(selb(2))=2
          sel(selb(3))=2
          IF(rwofil.NE.' ')
     +        CALL wrirwg(rwofil,objid,iobs,tutm,obscod,alpha,rmsa,
     +             resa,delta,rmsd,resd,smag,rmsmag,sel,n,rmsb*secrad)
      ELSE
          WRITE(*,231) name(1:ln)
          IF(unirep.GT.0) WRITE(unirep,231) name(1:ln)
          comele='FAIL'
      END IF
 231  FORMAT(4X,'IOD(',A,'): FAILED')
 232  FORMAT(4X,'IOD(',A,'): solution with ',A,' method (no noise): ',
     +           'RMS=',1P,E10.2,' arcsec')
 233  FORMAT(4X,'IOD(',A,'): solution with ',A,' method (with noise): ',
     +           'RMS=',1P,E10.2,' arcsec')

      END

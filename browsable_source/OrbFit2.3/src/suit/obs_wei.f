c ====================(part of) MODULE obs_wei========================
c CONTAINS
c SUBROUTINES
c OBSERVATIONS:
c    ERROR MODEL
c ***rdorms	read a-priori standard deviation of observations
c ***wrorms	writes a-priori standard deviation of observations
c ***wrores	writes a-priori standard deviation of observations and fit residuals
c inobs         high level observations input
c wrirwg        writes residuals (for Gauss method only) 
c obsrms	computation of a-priori RMS of observations
c astrow        computation of astrometry RMS based upon past performances
c *rrmscl       read from a file and store RMS classes
c *crmscl	computes a-priori observation RMS based on known classes
c magrms        computation of photometry RMS based upon past performances
c rearwo        reads residuals, weights and observations file .rwo
c wrirwo        writes residuals, weights and observations file .rwo
c wrirms        writes .rwo file, but without residuals
c statcod       fix for alphanumeric station code
c codestat       " 
c addobs        updates .rwo files using input .obs, .rad files
c fitwgt        computes weights as a function of standard deviations
c   INPUT
c mpcin		input of astrometric observations (MPC format)
c mpctr		transformation of an astrometric observation (MPC format)
c mpcrad        transformation of a radar observation (MPC format)
c jplin		input of radar observations (JPL format)
c jplrad        transformation of a radar observation (JPL format)
c sessag	transform an angle into sessagesimal notation
c astrad        gets radius of asteroid for radar
c
c
* Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: November 3, 1997
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         O B S R M S                           *
*  *                                                               *
*  *         Computation of a-priori RMS of observations           *
*  *    from their accuracy (= number of digits in input file)     *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    IOBS      -  Observation type: 1000+x=astrometry 
*                        2000+x radar 3000+x satellite 
*           IDSTA     -  Observatory code
*           TDT       -  Time (MJD, TDT)
*           ACCT      -  Accuracy of time (d)
*           ACCA      -  Accuracy of right ascension (rad)
*           ACCD      -  Accuracy of declination (rad)
*           SMAG      -  Apparent magnitude and color (string)
*           N         -  Number of observations
*
* OUTPUT:   RMSA      -  A-priori RMS of right ascension (rad)
*           RMSD      -  A-priori RMS of declination (rad)
*           RMSMAG    -  A-priori RMS of magnitude (magn)
*
      SUBROUTINE obsrms(iobs,idsta,tdt,acct,acca,accd,smag,
     +           rmsa,rmsd,rmsmag,n)
      IMPLICIT NONE

      INTEGER n,idsta(n),iobs(n)
      CHARACTER*6 smag(n)
      DOUBLE PRECISION tdt(n),acct(n),acca(n),accd(n),
     +              rmsa(n),rmsd(n),rmsmag(n)

      INTEGER i
      DOUBLE PRECISION rmsrmi,rmsvmi
c Get default magnitude weights from a function
      DOUBLE PRECISION magrms
      CHARACTER*1 typ

c      INCLUDE 'trig.h'
      INCLUDE 'jplhdr.h'
      DO 1 i=1,n
         IF(iobs(i)/1000.eq.1)THEN
c astrometric (telescope) observations; get weight
            typ=char(iobs(i)-1000)
            CALL astrow(typ,tdt(i),idsta(i),acca(i),accd(i),
     +               rmsa(i),rmsd(i))
c magnitude weigthing
            rmsmag(i)=magrms(smag(i),tdt(i),idsta(i),typ)
         ELSEIF(iobs(i)/1000.eq.2)THEN
c weighting of radar data is given with observations, but there is 
c minimum credible (for a given acccuracy of the models)
            rmsrmi=0.0d0/au
            rmsvmi=0.0d0/au
            rmsa(i)=MAX(acca(i),rmsrmi)
            rmsd(i)=MAX(accd(i),rmsvmi)
c            rmsa(i)=acca(i)
c            rmsd(i)=accd(i)
            rmsmag(i)=-1.d0
         ELSE
            WRITE(*,*)'obsrms: obs. type not known ', 
     +           tdt(i),' ',iobs(i),' ',i
            STOP
        ENDIF 
 1    CONTINUE
      END
* Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: December 5, 1997
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         R D O R M S                           *
*  *                                                               *
*  *      Read a-priori standard deviation of observations         *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    FILE      -  Input file name
* OUPUT:    TUTM      -  Time (MJD, UTM)
*           OBSCOD    -  Observatory code
*           RMSA      -  A-priori RMS of right ascension (rad)
*           RMSD      -  A-priori RMS of declination (rad)
*           SEL       -  Selection index (0=don't use; 1=use for fit;
*                        2=use for fit & Gauss method)
*           N         -  Number of observations
*
* On output, RMSA, RMSD and SEL of observations found in the file
* are updated; if the input file does not exists, returns without
* any change
*
      SUBROUTINE rdorms(file,tutm,obscod,rmsa,rmsd,sel,n)
      IMPLICIT NONE

      CHARACTER*(*) file
      INTEGER n,obscod(n),sel(n)
      DOUBLE PRECISION tutm(n),rmsa(n),rmsd(n)

      INCLUDE 'trig.h'
      INCLUDE 'parcmc.h'

      INTEGER unit,nr,lf,year,month,obsc1,iday,i,sel1
      DOUBLE PRECISION day,ra1,rd1,tutm1,eps
      CHARACTER*200 rec
      LOGICAL found

      INTEGER lench
      DOUBLE PRECISION tjm1
      EXTERNAL lench,tjm1

      DATA eps/5.D-6/

      INQUIRE(FILE=file,EXIST=found)
      IF(.NOT.found) RETURN
      CALL filopn(unit,file,'OLD')
      nr=0

 1    CONTINUE
      READ(unit,100,END=10) rec
 100  FORMAT(A)
      nr=nr+1
      IF(rec(1:1).EQ.comcha) GOTO 1
      READ(rec,*,ERR=20) year,month,day,obsc1,ra1,rd1,sel1
      iday=day
      tutm1=tjm1(iday,month,year,0.d0)+(day-iday)
      DO 3 i=1,n
      IF(obsc1.NE.obscod(i)) GOTO 3
      IF(ABS(tutm1-tutm(i)).GT.eps) GOTO 3
      rmsa(i)=ra1*radsec
      rmsd(i)=rd1*radsec
      sel(i)=sel1
      GOTO 1
 3    CONTINUE
      GOTO 1

 10   CONTINUE
      CALL filclo(unit,' ')
      RETURN

 20   CONTINUE
      lf=lench(file)
      WRITE(*,101) file(1:lf),nr
 101  FORMAT(' Input error (file "',A,'", record',I5,')')
      STOP '**** rdorms: abnormal end ****'

      END
* Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: November 7, 1997
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         W R O R M S                           *
*  *                                                               *
*  *       Writes a-priori standard deviation of observations      *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    FILE      -  Output file name
*           TUTM      -  Time (MJD, UTM)
*           OBSCOD    -  Observatory code
*           RMSA      -  A-priori RMS of right ascension (rad)
*           RMSD      -  A-priori RMS of declination (rad)
*           SEL       -  Selection indicator (0=don't use; 1=use for fit;
*                        2=use for fit & Gauss method)
*           N         -  Number of observations
*
      SUBROUTINE wrorms(file,tutm,obscod,rmsa,rmsd,sel,n)
      IMPLICIT NONE

      CHARACTER*(*) file
      INTEGER n,obscod(n),sel(n)
      DOUBLE PRECISION tutm(n),rmsa(n),rmsd(n)

      INCLUDE 'trig.h'
      INCLUDE 'parcmc.h'

      INTEGER unit,i,iday,month,year
      DOUBLE PRECISION day,hour

      CALL filopn(unit,file,'UNKNOWN')
      WRITE(unit,101) comcha
 101  FORMAT(A1,'YYYY MM DD.dddddd  OBS   rms(RA)"    rms(DEC)" S')

      DO 1 i=1,n
      CALL mjddat(tutm(i),iday,month,year,hour)
      day=iday+hour/24
      IF(sel(i).GT.9) STOP '**** wrorms: sel > 9 ****'
      WRITE(unit,100) year,month,day,obscod(i),
     +                rmsa(i)*secrad,rmsd(i)*secrad,sel(i)
 100  FORMAT(I5,I3,F10.6,I5,1P,2E12.4,I2)
 1    CONTINUE

      CALL filclo(unit,' ')

      END
* Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: November 7, 1997
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         W R O R E S                           *
*  *                                                               *
*  *       Writes a-priori standard deviation of observations      *
*  *                       and fit residuals                       *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    FILE      -  Output file name
*           TUTM      -  Time (MJD, UTM)
*           OBSCOD    -  Observatory code
*           RMSA      -  A-priori RMS of right ascension (rad)
*           RMSD      -  A-priori RMS of declination (rad)
*           SEL       -  Selection indicator (0=don't use; 1=use for fit;
*                        2=use for fit & Gauss method)
*           RESA      -  Residuals in right ascension (rad)
*           RESD      -  Residuals in declination (rad)
*           N         -  Number of observations
*
      SUBROUTINE wrores(file,tutm,obscod,rmsa,rmsd,sel,resa,resd,n)
      IMPLICIT NONE

      CHARACTER*(*) file
      INTEGER n,obscod(n),sel(n)
      DOUBLE PRECISION tutm(n),rmsa(n),rmsd(n),resa(n),resd(n)

      INCLUDE 'trig.h'
      INCLUDE 'parcmc.h'

      INTEGER unit,i,iday,month,year
      DOUBLE PRECISION day,hour,ra,rd
      CALL filopn(unit,file,'UNKNOWN')
      WRITE(unit,110) comcha
 110  FORMAT(A1,'YYYY MM DD.dddddd  OBS   rms(RA)"    rms(DEC)" S',
     +          '    res(RA)"  res(DEC)"')

      DO 1 i=1,n
      CALL mjddat(tutm(i),iday,month,year,hour)
      day=iday+hour/24
      ra=resa(i)*secrad
      rd=resd(i)*secrad
      IF(sel(i).GT.9) STOP '**** wrores: sel > 9 ****'
      IF(MAX(ABS(ra),ABS(rd)).GE.999.D0) THEN
          WRITE(unit,100) year,month,day,obscod(i),
     +                    rmsa(i)*secrad,rmsd(i)*secrad,sel(i),ra,rd
      ELSE
          WRITE(unit,101) year,month,day,obscod(i),
     +                    rmsa(i)*secrad,rmsd(i)*secrad,sel(i),ra,rd
      END IF
 100  FORMAT(I5,I3,F10.6,I5,1P,2E12.4,I2,1P,2E11.3)
 101  FORMAT(I5,I3,F10.6,I5,1P,2E12.4,I2,0P,2F11.3)
 1    CONTINUE

      CALL filclo(unit,' ')

      END
* Copyright (C) 1998 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: January 15, 1998; corrected AM/ZK March 16, 1998
* Modified 2 Dec. 1998 by Steven Chesley
* Changes: make optional to set w=0, correct for spherical metric
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         F I T W G T                           *
*  *                                                               *
*  *      Computation of observation weights for orbital fit       *
*  *                   from their a-priori RMS                     *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    RMSA      -  A-priori RMS of right ascension (rad)
*           RMSD      -  A-priori RMS of declination (rad)
*           DEL       -  Declination
*           SEL       -  Selection index
*           IOBS      -  Observation type: 2002,2003 are 1-dim
*           N         -  Number of observations
*           DISCARD   -  TRUE to set w=0 if sel=0
*
* OUTPUT:   W         -  Weights
*
      SUBROUTINE fitwgt(rmsa,rmsd,del,sel,iobs,w,n,discard)
      IMPLICIT NONE

      INTEGER n,sel(n),iobs(n)
      DOUBLE PRECISION rmsa(n),rmsd(n),del(n),w(2*n)
      LOGICAL discard

      INTEGER i

      DO i=1,n
         IF(sel(i).EQ.0 .AND. discard) THEN
            w(2*i-1)=0
            w(2*i)=0
         ELSEIF(iobs(i)/1000.eq.1)THEN
            w(2*i-1)=(cos(del(i))/rmsa(i))**2
            w(2*i)=1.d0/rmsd(i)**2
         ELSEIF(iobs(i).eq.2001.or.iobs(i).eq.2101)THEN
            w(2*i-1)=1.d0/rmsa(i)**2
            w(2*i)=1.d0/rmsd(i)**2
         ELSEIF(iobs(i).eq.2002.or.iobs(i).eq.2102)THEN
            w(2*i-1)=1.d0/rmsa(i)**2
            w(2*i)=0.d0
         ELSEIF(iobs(i).eq.2003.or.iobs(i).eq.2103)THEN
            w(2*i-1)=0.d0
            w(2*i)=1.d0/rmsd(i)**2
         END IF
      ENDDO

      END
* Copyright (C) 1998 by OrbFit Consortium
* Version: December 15, 1997 Steven Chesley
* Revised by Genny in Jun 20, 2001 to handle designations as K00Sa3P
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         I A U C O D                           *
*  *                                                               *
*  * Computes official IAU code from MPC-style packed designation  *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    MPCCOD    -  MPC-style code as used in observation archives
*                        5 characters for number, 7 characters for
*                        packed provisional designation 
*
* OUTPUT:   IAUDES    -  IAU code
*           ERROR     -  Error flag (cannot understand input code)
*
      SUBROUTINE iaucod(mpccod,iaudes,error)
      IMPLICIT NONE

      CHARACTER*(*) iaudes,mpccod
      LOGICAL error

      INTEGER ln,i,temp
      CHARACTER*2 head
      CHARACTER*3 tail

      INTEGER lench
      LOGICAL isnum
      EXTERNAL lench,isnum

      CHARACTER*12 numfield, desfield

      error=.false.
      iaudes=' '

      ln=lench(mpccod)
      IF(ln.LE.0) GOTO 10
      
      numfield=mpccod(1:5)
      call rmsp(numfield,ln)
* Numbered asteroids
      if(ln.ne.0) then
         iaudes=numfield
         do i=1,ln-1
            if (iaudes(i:i).eq.'0') then
               iaudes(i:i)=' '
            else
               goto 123
            endif
         enddo  
 123     call rmsp(iaudes,ln)
         return 
      endif
      
* Unnumbered asteroids
      desfield=mpccod(6:12)

      if(desfield(3:3).eq.'S')then
c Survey Asteroid
         iaudes=desfield(4:7)//desfield(1:1)//'-'//desfield(2:2)

      elseif (desfield(1:1).eq.'I' .or. 
     +        desfield(1:1).eq.'J' .or. 
     +        desfield(1:1).eq.'K')then
c Temporary designation
         if(desfield(5:6).eq.'00')then
c           1999AA = J99A00A
            tail=''
         elseif(desfield(5:5).eq.'0')then
c           1999AA1 = J99A01A
            tail=desfield(6:6)
         elseif(isnum(desfield(5:5)))then
c           1999AA12 = J99A12A
            tail=desfield(5:6)
         else
c           1999AA103 = J99AA3A
            temp=ichar(desfield(5:5))-55
c           1999AA363 = J99Aa3A
            if (temp.gt.35) temp = temp - 6
            write(head,103) temp
            tail=head//desfield(6:6)
         endif
         temp=ichar(desfield(1:1))-55
         write(head,103) temp
 103     format(I2)
         iaudes=head//desfield(2:3)//desfield(4:4)//desfield(7:7)//tail
      else
c Unknown type
         write(*,*)'cannot understand MPC designation: ',mpccod
         goto 10
      endif
      return

 10   CONTINUE
      iaudes=mpccod
      error=.true.

      END
* Copyright (C) 1998 by OrbFit Consortium
* Version: May 2000 AM MES
* Revised by Genny in Jun 20, 2001 to handle designations as K00Sa3P
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         I A U C O D 2                         *
*  *                                                               *
*  * Computes official IAU code from MPC-style packed designation  *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    MPCCOD    -  MPC-style packed 7 character code as used in orbit archives
*                        5 digits for numbered, 7 characters for unnumbered
*
* OUTPUT:   IAUDES    -  IAU code
*           ERROR     -  Error flag (cannot understand input code)
*
      SUBROUTINE iaucod2(mpccod,iaudes,error)
      IMPLICIT NONE

      CHARACTER*7 mpccod
      CHARACTER*9 iaudes
      LOGICAL error

      INTEGER ln,i,temp
      CHARACTER*2 head
      CHARACTER*3 tail

      INTEGER lench
      LOGICAL isnum
      EXTERNAL lench,isnum

      CHARACTER*5 numfield
      CHARACTER*7 desfield

      error=.false.
      iaudes=' '

      ln=lench(mpccod)
      IF(ln.LE.0) GOTO 10
      IF(ln.eq.5)THEN
* Numbered asteroids
         numfield=mpccod(1:5)
         call rmsp(numfield,ln)
         iaudes=numfield
         do i=1,ln-1
            if (iaudes(i:i).eq.'0') then
               iaudes(i:i)=' '
            else
               goto 123
            endif
         enddo  
 123     call rmsp(iaudes,ln)
         DO i=ln+1,9
           iaudes(i:i)='w'
         ENDDO
         return 
       ELSEIF(ln.eq.7)THEN      
* Unnumbered asteroids
          desfield=mpccod(1:7)
          if(desfield(3:3).eq.'S')then
c Survey Asteroid
          iaudes=desfield(4:7)//desfield(1:1)//'-'//desfield(2:2)//'ww'

          elseif (desfield(1:1).eq.'I' .or. 
     +            desfield(1:1).eq.'J' .or. 
     +            desfield(1:1).eq.'K')then
c Temporary designation
             if(desfield(5:6).eq.'00')then
c           1999AA = J99A00A
                tail='www'
             elseif(desfield(5:5).eq.'0')then
c           1999AA1 = J99A01A
                tail=desfield(6:6)//'ww'
             elseif(isnum(desfield(5:5)))then
c           1999AA12 = J99A12A
                tail=desfield(5:6)//'w'
             else
c           1999AA103 = J99AA3A
                temp=ichar(desfield(5:5))-55
c           1999AA363 = J99Aa3A
                if (temp.gt.35) temp = temp - 6
                write(head,103) temp
                tail=head//desfield(6:6)
             endif
             temp=ichar(desfield(1:1))-55
             write(head,103) temp
 103         format(I2)
         iaudes=head//desfield(2:3)//desfield(4:4)//desfield(7:7)//tail
          else
c Unknown type
             write(*,*)'cannot understand MPC designation: ',mpccod
             goto 10
          endif
          return
      ELSE
c wrong length of designator
          write(*,*) 'designation ',mpccod,' not understood'          
      ENDIF
 10   CONTINUE
      iaudes=mpccod
      error=.true.

      END






* Copyright (C) 1998 by Steven Chesley (chesley@dm.unipi.it)
* Version: Dec. 15, 1998
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         W R I R W O                           *
*  *                                                               *
*  *       Writes a-priori standard deviation of observations      *
*  *                       and fit residuals                       *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    FILE      -  Output file name
*           OBJID     -  IAU Identifier for each observation
*           IOBS      -  Observation type 1=a,d 2=r,rdot 3=r 4=rdot 5=sat
*           TUTM      -  Time (MJD, UTM)
*           OBSCOD    -  Observatory code
*           ALPHA     -  Right Ascension (radians)
*           RMSA      -  A-priori RMS of right ascension (rad)
*           RESA      -  Residuals in right ascension (rad) (aka O-C)
*           DELTA     -  Declination (radians)
*           RMSD      -  A-priori RMS of declination (rad)
*           RESD      -  Residuals in declination (rad)
*           SMAG      -  Magnitude observations (string)
*           RMSMAG       A-priori RMS of magnitude
*           RESMAG    -  Residuals in magnitude
*           RMSH      -  RMS of residuals in magnitude
*           SEL       -  Selection indicator (0=don't use; 1=use for fit;
*                              2=use for fit & Gauss method)
*           CHI2      -  CHI**2 value for each observation
*           N         -  Number of observations
*           RMSRES    -  RMS of the residuals
*
      SUBROUTINE wrirwo(file,objid,iobs,tutm,obscod,alpha,rmsa,resa,
     +     delta,rmsd,resd,
     +     smag,rmsmag,resmag,rmsh,
     +     sel,chi2,n,rmsres)
      IMPLICIT NONE

      CHARACTER*(*) file
      INTEGER n,obscod(n),sel(n),iobs(n)
      CHARACTER*(*) objid(n)      
      DOUBLE PRECISION tutm(n),alpha(n),delta(n),rmsa(n),rmsd(n)
      DOUBLE PRECISION resa(n),resd(n),chi2(n),rmsres
c magnitudes (string), a priori rms'
      CHARACTER*6 smag(n)
      DOUBLE PRECISION rmsmag(n)
c fit residuals, rms
      DOUBLE PRECISION resmag(n),rmsh
c =================END INTERFACE=============================
      INCLUDE 'trig.h'
      INCLUDE 'parcmc.h'
      INCLUDE 'jplhdr.h'
      INCLUDE 'parobx.h'

      INTEGER unit,i,is(nobx),ln
      INTEGER iyear,imonth,iday,ihour,imin,isec,ideg
      DOUBLE PRECISION day,hour,minu,sec,resnor
      CHARACTER*19 magstri,tmpstr
      CHARACTER*30 rastri,rdstri
      CHARACTER*37 rstri
      CHARACTER*3 radtyp
      CHARACTER*1 obstr,signo
      LOGICAL radar
      INTEGER truncat,iobcur,iotr,iore
c fix for exadecimal station code      
      CHARACTER*3 obsstr

      CALL filopn(unit,file,'UNKNOWN')

c check data set first
      radar=.false.
      DO i=1,n
         IF(iobs(i)/1000.eq.2)THEN
            radar=.true.
         ELSEIF(iobs(i)/1000.ne.1)THEN
            WRITE(*,*)'wrirwo: obs.type ',iobs(i),' unknown, rec.no=',i
            STOP
         ENDIF
      ENDDO
c ========= SORT OBSERVATIONS BY TIME ==============
      call heapsort(tutm,n,is)
c ========= HANDLE ASTROMETRY OBSERVATIONS =========
c astrometry header
      IF(rmsres.GT.0.D0) WRITE(unit,110) comcha,rmsres
 110  FORMAT(A1,'RMS of orbit residuals = ',F8.3)

      IF(rmsh.gt.0.d0)WRITE(unit,111) comcha,rmsh
 111  FORMAT(A1,'RMS of mag residuals = ',F5.2)

      WRITE(unit,120) comcha
 120  FORMAT(A1,'++ OBJ ++ OBS +++++ DATE ++++++',
     +     '  +++++ RIGHT ASCENSION ++++++',
     +     '  ++++++ DECLINATION ++++++++',
     +     '  +++++ APP MAG +++++',
     +     '  ++++ QUALITY +++')
      WRITE(unit,221) comcha
 221  FORMAT(A1,'+ DESIG + TYP YYYY MM DD.dddddd',
     +     '  HH MM SS.sss   rms     resid',
     +     '  DD MM SS.ss   rms     resid',
     +     '  MAG COL rms   resid',
     +     '  OBS     CHI  SEL')
c ========= ASTROMETRY LOOP==========================
      DO   i=1,n
      IF(iobs(is(i))/1000.eq.1)THEN
         obstr=char(iobs(is(i))-1000)
c convert time
         CALL mjddat(tutm(is(i)),iday,imonth,iyear,hour)
         day=iday+hour/24.d0
c convert RA
         CALL sessag(alpha(is(i))*degrad/15.d0,signo,ihour,imin,sec)
         IF(signo.eq.'-')STOP 'wrirwo error: negative right ascension.'
c prepare RA string
         resnor=resa(is(i))*secrad*cos(delta(is(i)))
         WRITE(tmpstr,FMT='(F6.3)') sec
         CALL rmsp(tmpstr,ln)
         IF(ln.lt.6)tmpstr='0'//tmpstr
         IF(abs(resnor).gt.999.d0)THEN
            WRITE(rastri,131)ihour,imin,tmpstr, 
     +           rmsa(is(i))*secrad,resnor
 131        FORMAT(2x,I2.2,1x,I2.2,1x,A6,F7.2,1P,E9.1)
         ELSE
            WRITE(rastri,130)ihour,imin,tmpstr,
     +           rmsa(is(i))*secrad,resnor
 130        FORMAT(2x,I2.2,1x,I2.2,1x,A6,F7.2,F9.3)
         ENDIF
c convert DEC
         CALL sessag(delta(is(i))*degrad,signo,ideg,imin,sec)
c prepare DEC string
         WRITE(tmpstr,FMT='(F5.2)') sec
         CALL rmsp(tmpstr,ln)
         IF(ln.lt.5)tmpstr='0'//tmpstr
         IF(abs(resd(is(i))*secrad).gt.999.d0)THEN
            WRITE(rdstri,171)signo,ideg,imin,tmpstr,
     +           rmsd(is(i))*secrad,resd(is(i))*secrad
 171        FORMAT(1x,A1,I2.2,1x,I2.2,1x,A5,F7.2,1P,E9.1)
         ELSE
            WRITE(rdstri,170)signo,ideg,imin,tmpstr,
     +           rmsd(is(i))*secrad,resd(is(i))*secrad
 170        FORMAT(1x,A1,I2.2,1x,I2.2,1x,A5,F7.2,F9.3)
         ENDIF
c prepare MAG string
         IF(rmsmag(is(i)).lt.0)THEN
            WRITE(magstri,121)smag(is(i))
 121        FORMAT(a6,13x)
         ELSEIF(resmag(is(i)).gt.1.d6)THEN
            WRITE(magstri,122)smag(is(i)),rmsmag(is(i))
 122        FORMAT(a6,1x,f5.2,7x)
         ELSE
            WRITE(magstri,123)smag(is(i)),rmsmag(is(i)),resmag(is(i))
 123        FORMAT(a6,1x,f5.2,1x,f6.2)
         ENDIF
c output  astrometry   
         WRITE(tmpstr,FMT='(F9.6)') day
         CALL rmsp(tmpstr,ln)
         IF(ln.lt.9)tmpstr='0'//tmpstr
c fix for exadecimal station code      
         CALL codestat(obscod(is(i)),obsstr)
         WRITE(unit,101) objid(is(i)),obstr,iyear,imonth,tmpstr,
     +        rastri,rdstri,magstri,
     +        obsstr,sqrt(chi2(is(i))),sel(is(i))
 101     FORMAT(1x,A9,2x,a1,2x,I4,1x,I2.2,1x,A9,
     +        A30,A29,2x,A19,
     +        2x,A3,f9.2,2X,I1,3x)
      ENDIF
      ENDDO
c ========= HANDLE RADAR OBSERVATIONS ==========
      If(.not.radar) GOTO 99
c radar header      
      WRITE(unit,128) comcha
 128  FORMAT(A1,'++ OBJ ++ OBS ++++++ DATE +++++++ ',
     +     '++++++++ RADAR RANGE/RANGE RATE +++++++++ ',
     +     '++++++ QUALITY +++++')
      WRITE(unit,228) comcha
 228  FORMAT(A1,'+ DESIG + TYP YYYY MM DD hh:mm:ss ',
     +     'TYP   KM or KM/DAY  a priori rms residual ',
     +     'TRX REC     CHI  SEL')
c ========= RADAR LOOP==========================
      DO i=1,n
         IF(iobs(is(i))/1000.eq.2)THEN
            IF(iobs(is(i))-2000.ge.100)THEN
c surface bounce
               iobcur=iobs(is(i))-2100
               obstr='r'
            ELSE
c already corrected to center of mass
               iobcur=iobs(is(i))-2000
               obstr='R'
            ENDIF
c convert time
            CALL mjddat(tutm(is(i)),iday,imonth,iyear,hour)
c convert hour to 12:12:12
            ihour=truncat(hour,1d-7)
            minu=(hour-ihour)*60.d0
            imin=truncat(minu,1d-5)
            sec=(minu-imin)*60.d0
            isec=truncat(sec,1d-3)
            IF(iobcur.eq.2)THEN
c range observation
               radtyp='DEL'               
               IF(abs(resa(is(i))*au).gt.999.d0)THEN
                  WRITE(rstri,143)alpha(is(i))*au,rmsa(is(i))*au,
     +                 resa(is(i))*au
               ELSE
                  WRITE(rstri,141)alpha(is(i))*au,rmsa(is(i))*au,
     +                 resa(is(i))*au
               ENDIF
            ELSEIF(iobcur.eq.3)THEN
c range-rate observation
               radtyp='DOP'               
               IF(abs(resd(is(i))*au).gt.999.d0)THEN
                  WRITE(rstri,143)delta(is(i))*au,rmsd(is(i))*au,
     +                 resd(is(i))*au
               ELSE
                  WRITE(rstri,141)delta(is(i))*au,rmsd(is(i))*au,
     +                 resd(is(i))*au
               ENDIF
            ELSE
               STOP '*** wrirwo.f: internal error(1) ***'
            ENDIF
 143        FORMAT(f16.5,1x,f9.5,1x,1p,e10.4)
 141        FORMAT(f16.5,1x,f9.5,1x,f10.5)
c find codes of two observatories
            iotr= obscod(is(i))/10000
            iore= obscod(is(i))-iotr*10000
c output radar data              
            WRITE(unit,102) objid(is(i)),obstr,iyear,imonth,iday,
     +           ihour,imin,isec,radtyp,rstri,
     +           iotr,iore,sqrt(chi2(is(i))),sel(is(i))
 102        FORMAT(1x,A9,2x,a1,2x,I4,1x,I2.2,1x,I2.2,1x,
     +           I2.2,':',I2.2,':',I2.2,1x,A3,1x,A37,
     +           1X,I3.3,1X,I3.3,f9.2,2X,I1,3x)
         ENDIF
      ENDDO
c ============================================
 99   CALL filclo(unit,' ')
      END
* Copyright (C) 1998 by Steven Chesley (chesley@dm.unipi.it)
* Version: Dec. 15, 1998
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         W R I R M S                           *
*  *                                                               *
*  *       Writes a-priori standard deviation of observations      *
*  *                       and NOT fit residuals                   *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    FILE      -  Output file name
*           OBJID     -  IAU Identifier for each observation
*           IOBS      -  Observation type 
*           TUTM      -  Time (MJD, UTM)
*           OBSCOD    -  Observatory code
*           ALPHA     -  Right Ascension (radians)
*           RMSA      -  A-priori RMS of right ascension (rad)
*           DELTA     -  Declination (radians)
*           RMSD      -  A-priori RMS of declination (rad)
*           SMAG      -  Magnitude observations (string)
*           RMSMAG       A-priori RMS of magnitude
*           SEL       -  Selection indicator (0=don't use; 1=use for fit;
*                              2=use for fit & Gauss method)
*           N         -  Number of observations
*
      SUBROUTINE wrirms(file,objid,iobs,tutm,obscod,alpha,rmsa,
     +     delta,rmsd,smag,rmsmag,sel,n)
      IMPLICIT NONE

      CHARACTER*(*) file
      INTEGER n,obscod(n),sel(n),iobs(n)
      CHARACTER*(*) objid(n)
      DOUBLE PRECISION tutm(n),alpha(n),delta(n),rmsa(n),rmsd(n)
c magnitudes (string), a priori rms'
      CHARACTER*6 smag(n)
      DOUBLE PRECISION rmsmag(n)
c =================END INTERFACE=============================
      INCLUDE 'trig.h'
      INCLUDE 'parcmc.h'
      INCLUDE 'jplhdr.h'
      INCLUDE 'parobx.h'

      INTEGER unit,i,is(nobx),ln
      INTEGER iyear,imonth,iday,ihour,imin,isec,ideg
      DOUBLE PRECISION day,hour,minu,sec
      CHARACTER*19 magstri,tmpstr
      CHARACTER*30 rastri,rdstri
      CHARACTER*37 rstri
      CHARACTER*3 radtyp
      CHARACTER*1 obstr,signo
      LOGICAL radar
      INTEGER truncat,iobcur,iotr,iore
c fix for exadecimal station code      
      CHARACTER*3 obsstr

      CALL filopn(unit,file,'UNKNOWN')

c check data set first
      radar=.false.
      DO i=1,n
         IF(iobs(i)/1000.eq.2)THEN
            radar=.true.
         ELSEIF(iobs(i)/1000.ne.1)THEN
            WRITE(*,*)'wrirms: obs.type ',iobs(i),' unknown, rec.no=',i
            STOP
         ENDIF
      ENDDO
c ========= SORT OBSERVATIONS BY TIME ==============
      call heapsort(tutm,n,is)
c ========= HANDLE ASTROMETRY OBSERVATIONS =========
      WRITE(unit,120) comcha
 120  FORMAT(A1,'++ OBJ ++ OBS +++++ DATE ++++++',
     +     '  +++++ RIGHT ASCENSION ++++++',
     +     '  ++++++ DECLINATION ++++++++',
     +     '  +++++ APP MAG +++++',
     +     '  ++++ QUALITY +++')
      WRITE(unit,221) comcha
 221  FORMAT(A1,'+ DESIG + TYP YYYY MM DD.dddddd',
     +     '  HH MM SS.sss   rms     resid',
     +     '  DD MM SS.ss   rms     resid',
     +     '  MAG COL rms   resid',
     +     '  OBS     CHI  SEL')
c ========= ASTROMETRY LOOP==========================
      DO   i=1,n
      IF(iobs(is(i))/1000.eq.1)THEN
         obstr=char(iobs(is(i))-1000)
c convert time
         CALL mjddat(tutm(is(i)),iday,imonth,iyear,hour)
         day=iday+hour/24.d0
******************
c convert RA
         CALL sessag(alpha(is(i))*degrad/15.d0,signo,ihour,imin,sec)
         IF(signo.eq.'-')STOP 'wrirms error: negative right ascension.'
c prepare RA string
         WRITE(tmpstr,FMT='(F6.3)') sec
         CALL rmsp(tmpstr,ln)
         IF(ln.lt.6)tmpstr='0'//tmpstr
         WRITE(rastri,130)ihour,imin,tmpstr,rmsa(is(i))*secrad
 130     FORMAT(2x,I2.2,1x,I2.2,1x,a6,F7.2,9x)
c convert DEC
         CALL sessag(delta(is(i))*degrad,signo,ideg,imin,sec)
c prepare DEC string
         WRITE(tmpstr,FMT='(F5.2)') sec
         CALL rmsp(tmpstr,ln)
         IF(ln.lt.5)tmpstr='0'//tmpstr
         WRITE(rdstri,170)signo,ideg,imin,tmpstr,rmsd(is(i))*secrad
 170     FORMAT(1x,A1,I2.2,1x,I2.2,1x,a5,F7.2,9x)
c prepare MAG string
         IF(rmsmag(is(i)).lt.0)THEN
            magstri=smag(is(i))//'            '
         ELSE
            WRITE(magstri,121)smag(is(i)),rmsmag(is(i))
 121        FORMAT(a6,1x,f5.2,1x,6x)
         ENDIF
c output  astrometry             
         WRITE(tmpstr,FMT='(F9.6)') day
         CALL rmsp(tmpstr,ln)
         IF(ln.lt.9)tmpstr='0'//tmpstr
c fix for exadecimal station code      
         CALL codestat(obscod(is(i)),obsstr)
         WRITE(unit,101) objid(is(i)),obstr,iyear,imonth,tmpstr,
     +        rastri,rdstri,magstri,
     +        obsstr,sel(is(i))
 101     FORMAT(1x,A9,2x,a1,2x,I4,1x,I2.2,1x,A9,
     +        A30,A29,2x,A19,
     +        2X,A3,9x,2X,I1,3x)
      ENDIF
      ENDDO
c ========= HANDLE RADAR OBSERVATIONS ==========
      If(.not.radar) GOTO 99
c radar header      
      WRITE(unit,128) comcha
 128  FORMAT(A1,'++ OBJ ++ OBS ++++++ DATE +++++++ ',
     +     '++++++++ RADAR RANGE/RANGE RATE +++++++++ ',
     +     '++++++ QUALITY +++++')
      WRITE(unit,228) comcha
 228  FORMAT(A1,'+ DESIG + TYP YYYY MM DD hh:mm:ss ',
     +     'TYP   KM or KM/DAY  a priori rms residual ',
     +     'TRX REC     CHI  SEL')
c ========= RADAR LOOP==========================
      DO i=1,n
         IF(iobs(is(i))/1000.eq.2)THEN
            IF(iobs(is(i))-2000.ge.100)THEN
c surface bounce
               iobcur=iobs(is(i))-2100
               obstr='r'
            ELSE
c already corrected to center of mass
               iobcur=iobs(is(i))-2000
               obstr='R'
            ENDIF
c convert time
            CALL mjddat(tutm(is(i)),iday,imonth,iyear,hour)
c convert hour to 12:12:12
            ihour=truncat(hour,1d-7)
            minu=(hour-ihour)*60.d0
            imin=truncat(minu,1d-5)
            sec=(minu-imin)*60.d0
            isec=truncat(sec,1d-3)
            IF(iobcur.eq.2)THEN
c range observation
               radtyp='DEL'               
               WRITE(rstri,141)alpha(is(i))*au,rmsa(is(i))*au
 141           FORMAT(f16.5,1x,f9.5,1x,10x)
            ELSEIF(iobcur.eq.3)THEN
c range-rate observation
               radtyp='DOP'               
               WRITE(rstri,141)delta(is(i))*au,rmsd(is(i))*au
            ELSE
               STOP '*** wrirms.f: internal error(1) ***'
            ENDIF
c find codes of two observatories
            iotr= obscod(is(i))/10000
            iore= obscod(is(i))-iotr*10000
c output radar data              
            WRITE(unit,102) objid(is(i)),obstr,iyear,imonth,iday,
     +           ihour,imin,isec,radtyp,rstri,
     +           iotr,iore,sel(is(i))
 102        FORMAT(1x,A9,2x,a1,2x,I4,1x,I2.2,1x,I2.2,1x,
     +           I2.2,':',I2.2,':',I2.2,1x,A3,1x,A37,
     +           1X,I3.3,1X,I3.3,9x,2X,I1,3x)
         ENDIF
      ENDDO
c ============================================
 99   CALL filclo(unit,' ')
      END
c ==============================================
c  statcode, codestat
c conversion from/to exadecimal to/from numeric code for observing stations
c answer to mess done by MPC in April 2002
c ==============================================
      SUBROUTINE statcode(obsstr,iobs)
      IMPLICIT NONE
c input
      CHARACTER*3 obsstr
c output 
      INTEGER iobs
c end interface
      CHARACTER*1 alfanum
      INTEGER hundreds, temp
      LOGICAL isnum
      READ(obsstr,100)iobs
 100  FORMAT(1x,i2)
      READ(obsstr,'(A1)')alfanum
      IF(isnum(alfanum))THEN
         READ(alfanum,'(I1)')hundreds
         iobs=iobs+100*hundreds
      ELSEIF(alfanum.eq.' ')THEN
         iobs=iobs
      ELSE
         temp=ichar(alfanum) - 55
         IF(temp.gt.35) temp = temp - 6
         iobs=iobs+100*temp
      ENDIF
      RETURN
      END
c================================================  
      SUBROUTINE codestat(iobs,obsstr)
      IMPLICIT NONE
c input
      INTEGER iobs
c output 
      CHARACTER*3 obsstr
c end interface
      CHARACTER*1 alfanum
      INTEGER hundreds, units
      units=mod(iobs,100)
      hundreds=(iobs-units)/100
      IF(hundreds.le.9)THEN
         WRITE(obsstr,'(I3.3)')iobs
      ELSE
         alfanum=char(55+hundreds)
         WRITE(obsstr,101)alfanum,units
 101     FORMAT(A1,I2.2)
      ENDIF
      RETURN
      END  
* Copyright (C) 1998 by Steven Chesley (chesley@dm.unipi.it)
* Version: Jan. 22, 1998
* Version: September 10, 1999
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         R E A R W O                           *
*  *                                                               *
*  *                 Reads observations, apriori rms,              *
*  *                   and post fit redisuals                      *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    FILE      -  Input file name
*           NLEF      -  max dimension for arrays
* OUTPUT:
*           OBJID     -  IAU Identifier for each observation
*           IOBS      -  Observation type 1000+x astrometry 2000+x radar
*           TAU       -  Time(MJD, TDT)
*           TUTM      -  Time (MJD, UTM)
*           OBSCOD    -  Observatory code
*           ALPHA     -  Right Ascension (radians)
*           RMSA      -  A-priori RMS of right ascension (rad)
*           RESA      -  Residuals in right ascension (rad) (aka O-C)
*           DELTA     -  Declination (radians)
*           RMSD      -  A-priori RMS of declination (rad)
*           RESD      -  Residuals in declination (rad)
*           SMAG      -  Magnitude observations (string)
*           RMSMAG       A-priori RMS of magnitude
*           RESMAG    -  Residuals in magnitude
*           SEL       -  Selection indicator (0=don't use; 1=use for fit;
*                              2=use for fit & Gauss method)
*           CHI2      -  CHI**2 value for each observation
*           N         -  Number of observations
*
      SUBROUTINE rearwo(file,objid,iobs,tau,tutm,obscod,
     +     alpha,rmsa,resa,
     +     delta,rmsd,resd,
     +     smag,rmsmag,resmag,
     +     sel,chi2,n,nlef)
      IMPLICIT NONE
c input file, space left in arrays
      CHARACTER*(*) file
      INTEGER nlef
c number of obs, observatory, selection flag, obs. type, object identifier
      INTEGER n,obscod(nlef),sel(nlef),iobs(nlef)
      CHARACTER*(*) objid(nlef)
c time of observation      
      DOUBLE PRECISION tutm(nlef),tau(nlef)
c observations and residuals
      DOUBLE PRECISION alpha(nlef),delta(nlef),rmsa(nlef),rmsd(nlef)
      DOUBLE PRECISION resa(nlef),resd(nlef),chi2(nlef)
c magnitudes (string), a priori rms'
      CHARACTER*6 smag(nlef)
      DOUBLE PRECISION rmsmag(nlef)
c fit residuals, rms
      DOUBLE PRECISION resmag(nlef)
c =================END INTERFACE=============================
      INCLUDE 'trig.h'
      INCLUDE 'parcmc.h'
      INCLUDE 'jplhdr.h'
c ===========units for err,pro,clo files ========
      INCLUDE 'proout.h'

      INTEGER unit,i,iday,month,year,deg,mindec,hr,minra,isec,imin,ihour
      DOUBLE PRECISION day,secra,secdec

      CHARACTER*140 rec,tmprec
      CHARACTER*1 obstyp,sign
      CHARACTER*3 scale
      CHARACTER*37 rstri
      CHARACTER*3 radtyp
      CHARACTER*9 chistr
      DOUBLE PRECISION chi,sec,sect
      INTEGER mjd,mjdt,ll,iotr,iore

      DOUBLE PRECISION tjm1
      EXTERNAL tjm1
c fix for exadecimal station code      
      CHARACTER*3 obsstr

      CALL filopn(unit,file,'UNKNOWN')

c ========= PROCESS RECORDS SEQUENTIALLY ===========
      n=0
      DO 1  i=1,nlef+20
         READ(unit,100,END=99,ERR=10)rec
 100     FORMAT(A)
c skip comments
         tmprec=rec
         CALL rmsp(tmprec,ll)
         IF(tmprec(1:1).eq.comcha) GOTO 1
c otherwise get observation code
         obstyp=rec(13:13)
         IF(obstyp.ne.'R'.and.obstyp.ne.'r'.and.obstyp.ne.'S')THEN
            n=n+1
            IF(n.gt.nlef) STOP 'rearwo: nobs > nobx.'
            iobs(n)=1000+ichar(obstyp)
c Read astrometry observation
            READ(rec,101,ERR=10) objid(n),year,month,day,
     +           hr, minra, secra, rmsa(n),resa(n),
     +           sign,deg,mindec,secdec,rmsd(n),resd(n),
     +           smag(n),rmsmag(n),resmag(n),
     +           obsstr,chistr,sel(n)
 101        FORMAT(1x,A9,5x,I4,I3,F10.6,
     +           2x,I2,1x,I2,F7.3,F7.2,G9.3,
     +           1x,A1,I2,1x,I2,F6.2,F7.2,G9.1,2x,
     +           a6,1x,f5.2,1x,f6.2,
     +           2X,A3,A9,2X,I1,3x)
            READ(chistr,FMT='(F9.2)',ERR=11)chi
            chi2(n)=chi**2
c fix for exadecimal station code      
            CALL statcode(obsstr,obscod(n))
            GOTO 12
 11         WRITE(ierrou,*)'rearwo: error in chi ',chistr
            WRITE(*,*)'rearwo: error in chi ',chistr
            chi2(n)=0.d0
            numerr=numerr+1
 12         CONTINUE
            IF(rec(100:112).eq.'             ') rmsmag(n)=-1.d0
c            IF(rmsmag(n).eq.0.d0) rmsmag(n)=-1.d0
c convert time
            IF(year.LT.1972) THEN
               scale='UT1'
            ELSE
               scale='UTC'
            ENDIF
            iday=day
            sec=(day-iday)*86400.d0
            mjd=nint(tjm1(iday,month,year,0.d0))
            CALL cnvtim(mjd,sec,scale,mjdt,sect,'TDT')
            tutm(n)=mjd+sec/86400.d0
            tau(n)=mjdt+sect/86400.d0 !test
c convert DEC
            delta(n)=(deg*3600.d0+mindec*60.d0+secdec)/secrad
            IF(sign.eq.'-')delta(n)=-delta(n)
            rmsd(n)=rmsd(n)/secrad
            resd(n)=resd(n)/secrad
c convert RA (residuals conversion depends upon delta)
            alpha(n)=15.d0*(hr*3600.d0+minra*60.d0+secra)/secrad
            rmsa(n)=rmsa(n)/secrad
            resa(n)=resa(n)/cos(delta(n))/secrad
c Radar Observation
         ELSEIF(obstyp.eq.'R'.or.obstyp.eq.'r')THEN
            n=n+1
            IF(n.gt.nlef) STOP 'rearwo: nobs > nobx.'
c Read radar observation
            READ(rec,102,ERR=10)objid(n),year,month,iday,ihour,
     +           imin,isec,radtyp,rstri,
     +           iotr,iore,chi,sel(n)
 102        FORMAT(1x,A9,2x,1x,2x,I4,1x,I2,1x,I2,1x,I2,':',I2,':',I2,1x,
     +           A3,1x,A37,
     +           1X,I3.3,1X,I3.3,f9.2,2X,I1,3x)
c warning: if a radar station has  a non-numeric code, we are in a mess.
            obscod(n)=iotr*10000+iore
            chi2(n)=chi**2
c           WRITE(*,102) objid(n),year,month,iday,ihour,imin,isec,
c    +           rstri,vstri,
c    +           iotr,iore,chi,sel(n)
c convert time
            sec=isec+60.d0*imin+3600.d0*ihour
            mjd=nint(tjm1(iday,month,year,0.d0))
            tutm(n)=mjd+sec/86400.d0
            IF(year.LT.1972) THEN
               scale='UT1'
            ELSE
               scale='UTC'
            ENDIF
            CALL cnvtim(mjd,sec,scale,mjdt,sect,'TDT')
            tau(n)=mjdt+sect/86400.d0
            IF(radtyp.eq.'DEL')THEN
c Prcess range string
               READ(rstri,141,ERR=10)alpha(n),rmsa(n),resa(n)
 141           FORMAT(f16.5,1x,f9.5,1x,f10.5)
               alpha(n)=alpha(n)/au
               rmsa(n)=rmsa(n)/au
               resa(n)=resa(n)/au
               iobs(n)=2002
               delta(n)=0.d0
               rmsd(n)=-1.d0
            ELSEIF(radtyp.eq.'DOP')THEN
c Prcess range rate string
               READ(rstri,141,ERR=10)delta(n),rmsd(n),resd(n)
               delta(n)=delta(n)/au
               rmsd(n)=rmsd(n)/au
               resd(n)=resd(n)/au
               iobs(n)=2003
               alpha(n)=0.d0
               rmsa(n)=-1.d0
            ELSE
               STOP'*** rearwo: internal error (1) ***'
            ENDIF
c photometry has no meaning
            smag(n)='      '
            rmsmag(n)=-1.d0
c surface bounce correction required
            IF(obstyp.eq.'r')iobs(n)=iobs(n)+100
         ELSE
c           Skip record if unknown type          
            WRITE(*,*) 'Unknown obs type ',obstyp,' at line ',
     +           i,' in ',file
         ENDIF
c spaghetti code: only go to 10 on error, else skip line 10.
         GOTO 1
 10      WRITE(ierrou,*) 'ERROR while reading line ',i,' of ',file
         WRITE(ierrou,*) 'skipping record: ',rec
         WRITE(*,*) 'ERROR while reading line ',i,' of ',file
         WRITE(*,*) 'skipping record: ',rec
         n=n-1
         numerr=numerr+1
 1    ENDDO
c ============================================
 99   CALL filclo(unit,' ')
      RETURN
      END
c =========================================
c  A D D O B S
c
c add information from .obs and .rad
c The number of observations can increase, from m to a maximum
c which is m+mt; however, the dimensioning of the vectors
c is nlef
      SUBROUTINE addobs(     
     +           objid,iobs,tau,tut,idsta,
     +           aln,rmsa,den,rmsd,smag,rmsmag,sel,m,
     +           objidt,iobst,taut,tutt,idstat,
     +           alnt,rmsat,dent,rmsdt,smagt,rmsmagt,mt,
     +           nlef,mnew,change)
      IMPLICIT NONE
c dummy variables length
      INTEGER nlef
c ===============================================
c change flag, new obs. number
      LOGICAL change
      INTEGER mnew
c ===== observational data ===========================
c observation number
      INTEGER m
c observations: alpha, delta, time (ET and UT), station code, type
      DOUBLE PRECISION aln(nlef),den(nlef),tau(nlef),tut(nlef)
      INTEGER idsta(nlef),iobs(nlef)
c identifier, app. magnitude
      CHARACTER*9 objid(nlef)
      CHARACTER*6 smag(nlef) 
c selection flag 0=discard 1=select 2=prelim
      INTEGER sel(nlef)
c RMS of observation error, of magnitude
      DOUBLE PRECISION rmsa(nlef),rmsd(nlef),rmsmag(nlef),rmseps
c ===== observational data: temporary copy===========
c observation number
      INTEGER mt
c observations: alpha, delta, time (ET and UT), station code, type
      DOUBLE PRECISION alnt(mt),dent(mt),taut(mt),tutt(mt)
      INTEGER idstat(mt),iobst(mt)
c identifier, app. magnitude
      CHARACTER*9 objidt(mt)
      CHARACTER*6 smagt(mt) 
c RMS of observation error, of magnitude
      DOUBLE PRECISION rmsat(mt),rmsdt(mt),rmsmagt(mt)
c ===========================================
      INTEGER j,mj,findob,double
      LOGICAl chaobs
c error file, number of errors
      INCLUDE 'proout.h'
      INCLUDE 'trig.h'
      INCLUDE 'jplhdr.h'
c ==================================================
c monitor changes
      mnew=m
      change=.false. 
c scan supposedly new observations
      DO 1 j=1,mt
         mj=findob(tutt(j),tut,idstat(j),idsta,iobst(j),iobs,double,m)
         IF(mj.ne.0.and.double.eq.0)THEN
c     the observation was already there
            IF(chaobs(alnt(j),aln(mj),dent(j),den(mj),
     +           smagt(j),smag(mj),iobst(j),iobs(mj)) )THEN 
c ... but it is changed 
c               write(*,*)'change',j,mj,alnt(j),aln(mj),dent(j),den(mj),
c     +           smagt(j),smag(mj),iobst(j),iobs(mj)
               change=.true.
               aln(mj)=alnt(j)
               den(mj)=dent(j)
               smag(mj)=smagt(j)
               rmsa(mj)=rmsat(j)
               rmsd(mj)=rmsdt(j)
               rmsmag(mj)=rmsmagt(j)
               sel(mj)=1
               iobs(mj)=iobst(j)
               IF(ierrou.gt.0)THEN
                  WRITE(ierrou,*)'addobs: changed observ. at record ',mj
                  numerr=numerr+1
               ELSE
                  WRITE(*,*)'addobs: changed observation at record ',mj
               ENDIF
            ELSE
c do not give precedence to rwo weighting unless rwo weights are negative
c cutoff is 5 milliarcsec or 1 meter (radar) for new weights
               if(iobs(mj)/1000.eq.2)then
                  rmseps=0.0005/au
               elseif(iobs(mj)/1000.eq.1)then
                  rmseps=0.0055*radsec
               endif
               if(rmsa(mj).ge.0d0.and.
     +              abs(rmsa(mj)-rmsat(j)).gt.rmseps)then
                  change=.true.
c                  write(*,'(a,4g16.8)')'alf',rmseps,
c     +                 rmsa(mj)-rmsat(j),rmsa(mj),rmsat(j)
                  rmsa(mj)=rmsat(j) 
               endif  
               if(rmsd(mj).ge.0d0.and.
     +              abs(rmsd(mj)-rmsdt(j)).gt.rmseps)then
c                  write(*,'(a,4g16.8)')'del',rmseps,
c     +                 rmsd(mj)-rmsdt(j),rmsd(mj),rmsdt(j)
                   change=.true.
                  rmsd(mj)=rmsdt(j) 
               endif  
            ENDIF
         ELSEIF(mj.ne.0.and.double.ne.0)THEN
c the observation was already there, in double copy!
            IF(chaobs(alnt(j),aln(mj),dent(j),den(mj),smagt(j),smag(mj)
     +           ,iobst(j),iobs(mj)) ) THEN 
               IF(chaobs(alnt(j),aln(double),dent(j),den(double),
     +              smagt(j),smag(double),iobst(j),iobs(double)))THEN 
                  change=.true.
c double, and changed! human intervention required
                  WRITE(*,*)'addobs: double and changed'
                  WRITE(*,*)' records ',mj,' and ',double,' in .rwo'
                  WRITE(*,*)' record ',j,' in .obs'
c                STOP
               ELSE
c OK, it is the double
               ENDIF
            ELSE
c OK, it is the first one
            ENDIF  
         ELSEIF(mj.eq.0)THEN
c the observation is new: add it
           change=.true.
c           WRITE(*,*)'addobs: new observation at record ',j
           mnew=mnew+1
           aln(mnew)=alnt(j)
           den(mnew)=dent(j)
           smag(mnew)=smagt(j)
           rmsa(mnew)=rmsat(j)
           rmsd(mnew)=rmsdt(j)
           rmsmag(mnew)=rmsmagt(j)
           tau(mnew)=taut(j)
           tut(mnew)=tutt(j)
           idsta(mnew)=idstat(j)
           iobs(mnew)=iobst(j)
           sel(mnew)=1
           objid(mnew)=objidt(j)
        ENDIF
 1    ENDDO
      RETURN
      END
c =======================================================
c F I N D O B
c find an observation from a list, matching the given one
      INTEGER FUNCTION findob(tutt,tut,idstat,idsta,iobst,iobs,double,m)
      IMPLICIT NONE
c ============= INPUT =====================
c number of obs. record to be scanned, time, time to be found
      INTEGER m
      DOUBLE PRECISION tut(m),tutt
c station code vector, of the one to be found
      INTEGER idsta(m),idstat,iobs(m),iobst
c ============OUTPUT (impure function!) ===========
      INTEGER double
c error number, file
      INCLUDE 'proout.h'
c =========END INTERFACE=====================
      INTEGER j
      DOUBLE PRECISION epst
      PARAMETER (epst=1.d-8)
      findob=0
      double=0
      DO 1 j=1,m
        IF(abs(tutt-tut(j)).lt.epst.and.
     +      idsta(j).eq.idstat.and.iobs(j).eq.iobst)THEN
           IF(findob.ne.0)THEN
              IF(double.eq.0)THEN
                 double=j
                 IF(ierrou.gt.0)THEN
                 WRITE(ierrou,*)'findob: two same time',
     +                findob,j,tutt,idsta(j)
                 numerr=numerr+1
                 ELSE
                    WRITE(*,*)'findob: two same time',
     +                findob,j,tutt,idsta(j)
                 ENDIF
              ELSE
                 IF(ierrou.gt.0)THEN
                    WRITE(ierrou,*)'findob: three same time',
     +                   findob,double,j,tutt,idsta(j)
                 ELSE
                    WRITE(*,*)'findob: three same time',
     +                   findob,double,j,tutt,idsta(j)
                 ENDIF
c                STOP
              ENDIF
           ELSE
              findob=j
           ENDIF
        ENDIF
 1    ENDDO
      RETURN
      END
c =======================================================
c C H A O B S
c is an observation changed?
      LOGICAL FUNCTION chaobs(aln,alnt,dent,den,smagt,smag,iobs,iobst)
      IMPLICIT NONE
      DOUBLE PRECISION aln,alnt,den,dent
      CHARACTER*6 smag,smagt
      INTEGER iobst,iobs
      DOUBLE PRECISION epsa
      PARAMETER (epsa=1.d-9)
      chaobs=.false.
      IF(abs(aln-alnt).gt.epsa*abs(aln).or.
     +      abs(den-dent).gt.epsa*abs(den))chaobs=.true.
      IF(smagt.ne.smag.or.iobs.ne.iobst)chaobs=.true.
      RETURN
      END
c =========================================
c  A D D R W O
c
c add information from .rwo file (weights, selection flags)
c to the observations from .obs and .rad
c The number of observations cannot increase
      SUBROUTINE addrwo(
     +           objid,iobs,tut,idsta,
     +           aln,rmsa,den,rmsd,smag,rmsmag,sel,m,
     +           objidt,iobst,tutt,idstat,
     +           alnt,rmsat,dent,rmsdt,smagt,rmsmagt,selt,mt,
     +           change)
      IMPLICIT NONE
c logical change flag
      LOGICAl change
c ===== observational data ===========================
c observation number
      INTEGER m
c observations: alpha, delta, time (ET and UT), station code, type
      DOUBLE PRECISION aln(m),den(m),tut(m)
      INTEGER idsta(m),iobs(m)
c identifier, app. magnitude
      CHARACTER*9 objid(m)
      CHARACTER*6 smag(m) 
c selection flag 0=discard 1=select 2=prelim
      INTEGER sel(m)
c RMS of observation error, of magnitude
      DOUBLE PRECISION rmsa(m),rmsd(m),rmsmag(m),rmseps
c ===== observational data: temporary copy===========
c observation number
      INTEGER mt
c observations: alpha, delta, time (ET and UT), station code, type
      DOUBLE PRECISION alnt(mt),dent(mt),tutt(mt)
      INTEGER idstat(mt),iobst(mt)
c identifier, app. magnitude
      CHARACTER*9 objidt(mt)
      CHARACTER*6 smagt(mt) 
c selection flag 0=discard 1=select 2=prelim
      INTEGER selt(mt)
c RMS of observation error, of magnitude
      DOUBLE PRECISION rmsat(mt),rmsdt(mt),rmsmagt(mt)
c ===========================================
      INTEGER j,mj,findob,double
      LOGICAl chaobs
      INCLUDE 'trig.h'
      INCLUDE 'jplhdr.h'
c if all the same...
      change=.false.
c scan old weigts and selection flags, see if they match  observations
      DO 1 j=1,mt
        mj=findob(tutt(j),tut,idstat(j),idsta,iobst(j),iobs,double,m)
        IF(mj.ne.0.and.double.eq.0)THEN
c this observation is still present in .obs, .rad files
           IF(chaobs(alnt(j),aln(mj),dent(j),den(mj),smagt(j),smag(mj)
     +        ,iobst(j),iobs(mj)))THEN 
c ... but it is changed, thus leave selection flag=1 and default weight
              change=.true.
           ELSE
c do not give precedence to rwo weighting unless rwo weights are negative
c cutoff is 5 milliarcsec or 1 meter for new weights
c here we just set the change flag to accept the computed weights
              if(iobs(mj)/1000.eq.2)then
                 rmseps=0.0005/au
              elseif(iobs(mj)/1000.eq.1)then
                 rmseps=0.0055*radsec
              endif
              if(rmsat(j).ge.0d0.and.
     +             abs(rmsat(j)-rmsa(mj)).gt.rmseps)then
                 change=.true.
              endif  
              if(rmsdt(j).ge.0d0.and.
     +             abs(rmsdt(j)-rmsd(mj)).gt.rmseps)then
                 change=.true.
              endif  
c if no change to observation then preserve the 
c manually fixed weights and selection flags
              IF(rmsat(j).lt.0.d0)rmsa(mj)=rmsat(j)
              IF(rmsdt(j).lt.0.d0)rmsd(mj)=rmsdt(j)
              IF(rmsmagt(j).lt.0.d0)rmsmag(mj)=rmsmagt(j)
c selection flags are preserved anyway
              sel(mj)=selt(j)
           ENDIF
        ELSEIF(mj.ne.0.and.double.ne.0)THEN
c this observation is present in .obs, .rad files, in duplicate!
           IF(chaobs(alnt(j),aln(mj),dent(j),den(mj),smagt(j),smag(mj)
     +        ,iobst(j),iobs(mj)))
     +     THEN 
              IF(chaobs(alnt(j),aln(double),dent(j),den(double),
     +             smagt(j),smag(double),iobst(j),iobs(double)))THEN
                 change=.true.
c double, and changed! human intervention required
                 WRITE(*,*)'addobs: double and changed'
                 WRITE(*,*)' records ',mj,' and ',double,' in .obs'
                 WRITE(*,*)' record ',j,' in .rwo'
c                STOP
              ELSE
c OK, it is the duplicate
c it is the same, so preserve the weights and selection flags
                 rmsa(double)=rmsat(j)
                 rmsd(double)=rmsdt(j)
                 rmsmag(double)=rmsmagt(j)
                 sel(double)=selt(j)
              ENDIF
           ELSE
c OK, it is the first one
c it is the same, so preserve the weights and selection flags
              rmsa(mj)=rmsat(j)
              rmsd(mj)=rmsdt(j)
              rmsmag(mj)=rmsmagt(j)
              sel(mj)=selt(j)
           ENDIF
        ELSEIF(mj.eq.0)THEN
c if it is not found in .rwo, leave the default weights and selection flags
           change=.true.
        ENDIF
 1    ENDDO
c check if there are extra (added) observations in .obs file
c it might be better to loop on the .obs rather than the .rwo data
      IF (mt.lt.m) change=.true.
      RETURN
      END
c ===========================================
c INOBS
c observation input control routine
c reads sequentially the .rwo, .obs, .rad file and combines
c the data according to the precedence rule specified by precob
c
      SUBROUTINE inobs(obsdir,astna0,precob,objid,obs0,m,iobs,tau,
     +    aln,den,tut,idsta,sel,rmsa,rmsd,rmsmag,smag,nlef,iun20,change)
      IMPLICIT NONE
c ==============INPUT==================
c input directory (all files astna0.rwo, astna0.obs, astna0.rad must be there)
      CHARACTER*60 obsdir
c asteroid name
      CHARACTER*(*) astna0
c messages unit
      INTEGER iun20
c observation numbers: maximum, space left
      INCLUDE 'parobx.h'
      INTEGER nlef
c logical flag: .true. for overwrite .rwo, .false. for update .rwo
      LOGICAL precob
c =============OUTPUT================================
c successful input flag
      LOGICAL obs0,change
c ===== observational data ===========================
c observation number
      INTEGER m
c observations: alpha, delta, time (ET and UT), station code, type
      DOUBLE PRECISION aln(nlef),den(nlef),tau(nlef),tut(nlef)
      INTEGER idsta(nlef),iobs(nlef)
c identifier, app. magnitude
      CHARACTER*9 objid(nlef)
      CHARACTER*6 smag(nlef) 
c selection flag 0=discard 1=select 2=prelim
      INTEGER sel(nlef)
c RMS of observation error, of magnitude
      DOUBLE PRECISION rmsa(nlef),rmsd(nlef),rmsmag(nlef)
c ===========END INTERFACE=========================
c file names
      CHARACTER*77 file
      INTEGER lfile
      LOGICAL rwo,mpc,rad
c accuracy (?) of observational data from mpcin
      double precision acct(nobx),acca(nobx),accd(nobx)
c new obs. number
      INTEGER mnew,mr,nlefm
c ===== observational data: temporary copy===========
c observation number
      INTEGER mt
c observations: alpha, delta, time (ET and UT), station code, type
      DOUBLE PRECISION alnt(nobx),dent(nobx),taut(nobx),tutt(nobx)
      INTEGER idstat(nobx),iobst(nobx)
c identifier, app. magnitude
      CHARACTER*9 objidt(nobx)
      CHARACTER*6 smagt(nobx) 
c selection flag 0=discard 1=select 2=prelim
      INTEGER selt(nobx)
c RMS of observation error, of magnitude
      DOUBLE PRECISION rmsat(nobx),rmsdt(nobx),rmsmagt(nobx)
c ============residuals and other data from rearwo============
      DOUBLE PRECISION resa(nobx),resd(nobx),chi2(nobx)
      DOUBLE PRECISION resmag(nobx)
c sorting 
      INTEGER iperm(nobx)
c directory char
      INCLUDE 'sysdep.h'
c loop index
      INTEGER i,j
      INTEGER ld
      INTEGER lench
      EXTERNAL lench
c =============EXECUTION BEGINS======================
c  compute file name
      ld=lench(obsdir)
      IF(ld.GT.0) THEN
          IF(obsdir(ld:ld).EQ.dircha) THEN
              file=obsdir(1:ld)//astna0
          ELSE
              file=obsdir(1:ld)//dircha//astna0
          END IF
      ELSE
          file=astna0
      END IF
      CALL rmsp(file,lfile)
c existence of .rwo, .obs, .rad
      INQUIRE(file=file(1:lfile)//'.rwo',exist=rwo)
      INQUIRE(file=file(1:lfile)//'.obs',exist=mpc)
      INQUIRE(file=file(1:lfile)//'.rad',exist=rad)
c For now we will not do radar only orbits
      IF(.not.rwo .and. .not. mpc)THEN
         WRITE(*,*)'You must provide either a .obs or .rwo file',
     +        'in directory ',obsdir
         obs0=.false.
         RETURN
      ENDIF
c select operations mode
      IF(.not.rwo)THEN
         WRITE(*,*) 'No .rwo file, so reading .obs and/or .rad files.'
c there is no .rwo, so read .obs and/or .rad
         IF(mpc)THEN
c Input of astrometric observations from a file (MPC format)
            CALL mpcin(mpc,file(1:lfile)//'.obs',objid,iobs,tau,tut,
     +           aln,den,idsta,acct,acca,accd,smag,m,nlef)
            WRITE(*,*)'mpcin: ',m,' obs in ',file(1:lfile)//'.obs'
            WRITE(iun20,*)'mpcin: ',m,' obs in ',file(1:lfile)//'.obs'
         ENDIF
c read radar jpl data
         IF(rad)THEN
            nlefm=nlef-m
            CALL jplin(rad,file(1:lfile)//'.rad',objid(m+1),iobs(m+1),
     +           tau(m+1),tut(m+1),aln(m+1),den(m+1),idsta(m+1),
     +           acct(m+1),acca(m+1),accd(m+1),smag(m+1),rmsmag(m+1),mr,
     +           nlefm)
            m=mr+m
         ENDIF
         obs0=mpc.or.rad
         IF(.not.obs0)RETURN
c find weights for these; a priori RMS of astrometric observations
         CALL obsrms(iobs,idsta,tau,acct,acca,accd,smag,
     +             rmsa,rmsd,rmsmag,m)
c give default selection flag of 1
         DO j=1,m
           sel(j)=1
         ENDDO
c output data for possible manual fixing: create weights file
         CALL wrirms(file(1:lfile)//'.rwo',objid,iobs,tut,idsta,
     +           aln,rmsa,den,rmsd,smag,rmsmag,sel,m)
         change=.true.
c select between update and overwrite of .rwo
      ELSEIF(precob)THEN
         WRITE(*,*)file(1:lfile),
     +        '.rwo found but ALL obs will come from .obs/.rad files.'
c give the precedence to the observation files .obs and .rad 
c with respect to .rwo, which is overwritten
         IF(mpc)THEN
c Input of astrometric observations from a file (MPC format)
            CALL mpcin(mpc,file(1:lfile)//'.obs',objid,iobs,tau,tut,
     +           aln,den,idsta,acct,acca,accd,smag,m,nlef)
            WRITE(*,*)'mpcin: ',m,' obs in ',file(1:lfile)//'.obs'
            WRITE(iun20,*)'mpcin: ',m,' obs in ',file(1:lfile)//'.obs'
         ELSE
            m=0
         ENDIF
c read radar jpl data
         IF(rad)THEN
            nlefm=nlef-m
            CALL jplin(rad,file(1:lfile)//'.rad',objid(m+1),iobs(m+1),
     +           tau(m+1),tut(m+1),aln(m+1),den(m+1),idsta(m+1),
     +           acct(m+1),acca(m+1),accd(m+1),smag(m+1),rmsmag(m+1),mr,
     +           nlefm)
            m=mr+m
         ENDIF
         obs0=rad.or.mpc
c If no obs then object does not "exist", so rwo should mnot be read:
         if(.not.obs0)return
c find weights for these; a priori RMS of astrometric observations
         CALL obsrms(iobs,idsta,tau,acct,acca,accd,smag,
     +             rmsa,rmsd,rmsmag,m)
c give default selection flag of 1
         DO j=1,m
           sel(j)=1
         ENDDO
c read .rwo  anyway, but store in temporary array the data
         CALL rearwo(file(1:lfile)//'.rwo',
     +     objidt,iobst,taut,tutt,idstat,
     +     alnt,rmsat,resa,dent,rmsdt,resd,
     +     smagt,rmsmagt,resmag,
     +     selt,chi2,mt,nlef)
c recover informations from .rwo (weights, selection flags)
         CALL addrwo(
     +           objid,iobs,tut,idsta,
     +           aln,rmsa,den,rmsd,smag,rmsmag,sel,m,
     +           objidt,iobst,tutt,idstat,
     +           alnt,rmsat,dent,rmsdt,smagt,rmsmagt,selt,mt,change)
         IF(change)THEN
            WRITE(*,*)'There are new/changed obs. New numobs=',m
c output updated .rwo file, if there are new observations (erasing residuals)
            CALL wrirms(file(1:lfile)//'.rwo',objid,iobs,tut,idsta,
     +           aln,rmsa,den,rmsd,smag,rmsmag,sel,m)
         ELSE
            WRITE(*,*)'There are no updates in .obs or .rad files.'
         ENDIF
      ELSE
c give the precedence to .rwo, the .obs and .rad files are intended
c as additional observations only; data in .rwo are not erased, can only be
c changed
c read .rwo, and store data in final array
         WRITE(*,*)'Using .rwo file, but checking .obs,.rad for update.'
         CALL rearwo(file(1:lfile)//'.rwo',
     +     objid,iobs,tau,tut,idsta,
     +     aln,rmsa,resa,den,rmsd,resd,
     +     smag,rmsmag,resmag,
     +     sel,chi2,m,nlef)
         WRITE(*,*)'rearwo: ',m,' obs from  ',file(1:lfile)//'.rwo'
         WRITE(iun20,*)'rearwo: ',m,' obs from ',file(1:lfile)//'.rwo'
         IF(m.eq.0)THEN
            obs0=.false.
         ELSE
            obs0=.true.
         ENDIF
c if there are input data
         IF(mpc)THEN
c Input of astrometric observations into temporary from a file (MPC format)
            CALL mpcin(mpc,file(1:lfile)//'.obs',
     +           objidt,iobst,taut,tutt,
     +           alnt,dent,idstat,acct,acca,accd,smagt,mt,nobx)
            IF(.not.mpc)THEN
               WRITE(*,*) file(1:lfile)//'.obs is possibly corrupt. ',
     +              'Not using any data from this file.'
            ELSE
               WRITE(*,*)'mpcin:',mt,' obs from  ',file(1:lfile)//'.obs'
               WRITE(iun20,*)'mpcin:',mt,' from ',file(1:lfile)//'.obs'
            ENDIF
         ENDIF
c read radar jpl data
         IF(rad)THEN
            nlefm=nobx-mt
            CALL jplin(rad,file(1:lfile)//'.rad',objidt(mt+1),
     +           iobst(mt+1),taut(mt+1),tutt(mt+1),alnt(mt+1),
     +           dent(mt+1),idstat(mt+1),acct(mt+1),acca(mt+1),
     +           accd(mt+1),smagt(mt+1),rmsmag(mt+1),mr,nlefm)
            IF(.not.rad)THEN
               WRITE(*,*) file(1:lfile)//'.rad is possibly corrupt. ',
     +              'Not using any data from this file.'
            ELSE
c               WRITE(*,*)'mpcin:',mr,' obs from  ',file(1:lfile)//'.rad'
c               WRITE(iun20,*)'mpcin:',mr,' obs  ',file(1:lfile)//'.rad'
               mt=mr+mt
            ENDIF
         ENDIF
c add information from .obs and .rad
         IF(mpc.or.rad)THEN
            obs0=.true.
c find weights for these; a priori RMS of astrometric observations
            CALL obsrms(iobst,idstat,taut,acct,acca,accd,smagt,
     +           rmsat,rmsdt,rmsmagt,mt)
            CALL addobs(     
     +           objid,iobs,tau,tut,idsta,
     +           aln,rmsa,den,rmsd,smag,rmsmag,sel,m,
     +           objidt,iobst,taut,tutt,idstat,
     +           alnt,rmsat,dent,rmsdt,smagt,rmsmagt,mt,
     +           nlef,mnew,change)
            IF(change)THEN
               WRITE(*,*)'There are new/changed obs. New numobs=',mnew
c output updated .rwo file, if there are new observations (erasing residuals)
               CALL wrirms(file(1:lfile)//'.rwo',objid,iobs,tut,idsta,
     +              aln,rmsa,den,rmsd,smag,rmsmag,sel,mnew)
               m=mnew
            ELSE
               WRITE(*,*)'There are no updates in .obs or .rad files.'
            ENDIF
         else
c check for new weights anyway???
         ENDIF
      ENDIF
c get asteroid radius (if necessary) before returning 
      call astrad(objid,iobs,m)
c =======  sort data before returning =========
      call heapsort(tut,m,iperm)
c copy output into temp vectors
      do i=1,m
         objidt(i)=objid(i)
         iobst(i)=iobs(i)
         taut(i)=tau(i)
         alnt(i)=aln(i)
         dent(i)=den(i)
         tutt(i)=tut(i)
         idstat(i)=idsta(i)
         selt(i)=sel(i)
         rmsat(i)=rmsa(i)
         rmsdt(i)=rmsd(i)
         rmsmagt(i)=rmsmag(i)
         smagt(i)=smag(i)
      enddo
c copy back input output vectors in sorted order
      do i=1,m
         objid(i)=objidt(iperm(i))
         iobs(i)=iobst(iperm(i))
         tau(i)=taut(iperm(i))
         aln(i)=alnt(iperm(i))
         den(i)=dent(iperm(i))
         tut(i)=tutt(iperm(i))
         idsta(i)=idstat(iperm(i))
         sel(i)=selt(iperm(i))
         rmsa(i)=rmsat(iperm(i))
         rmsd(i)=rmsdt(iperm(i))
         rmsmag(i)=rmsmagt(iperm(i))
         smag(i)=smagt(iperm(i))
      enddo

      RETURN
      END
* Copyright (C) 1999 by Steven Chesley (chesley@dm.unipi.it)
*                    and Mario Carpino (carpino@brera.mi.astro.it)
* Version: February 10, 1999
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         W R I R W G                           *
*  *                                                               *
*  *       Writes a-priori standard deviation of observations      *
*  *         and fit residuals (for Gauss' method only)            *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    FILE      -  Output file name
*           OBJID     -  IAU Identifier for each observation
*           IOBS      -  Observation type 1=a,d 2=r,rdot 3=r 4=rdot 5=sat
*           TUTM      -  Time (MJD, UTM)
*           OBSCOD    -  Observatory code
*           ALPHA     -  Right Ascension (radians)
*           RMSA      -  A-priori RMS of right ascension (rad)
*           RESA      -  Residuals in right ascension (rad) (aka O-C)
*           DELTA     -  Declination (radians)
*           RMSD      -  A-priori RMS of declination (rad)
*           RESD      -  Residuals in declination (rad)
*           SMAG      -  Magnitude observations (string)
*           RMSMAG       A-priori RMS of magnitude
*           SEL       -  Selection indicator (0=don't use; 1=use for fit;
*                              2=use for fit & Gauss method)
*           N         -  Number of observations
*           RMSRES    -  RMS of the residuals
*
      SUBROUTINE wrirwg(file,objid,iobs,tutm,obscod,alpha,rmsa,resa,
     +                  delta,rmsd,resd,smag,rmsmag,sel,n,rmsres)
      IMPLICIT NONE

      CHARACTER*(*) file
      INTEGER n,obscod(n),sel(n),iobs(n)
      CHARACTER*(*) objid(n)      
      DOUBLE PRECISION tutm(n),alpha(n),delta(n),rmsa(n),rmsd(n)
      DOUBLE PRECISION resa(n),resd(n),rmsres
c magnitudes (string), a priori rms'
      CHARACTER*6 smag(n)
      DOUBLE PRECISION rmsmag(n)
c =================END INTERFACE=============================
      INCLUDE 'trig.h'
      INCLUDE 'parcmc.h'
      INCLUDE 'jplhdr.h'
      INCLUDE 'parobx.h'

      INTEGER unit,i,is(nobx),ln
      INTEGER iyear,imonth,iday,ihour,imin,isec,ideg
      DOUBLE PRECISION day,hour,minu,sec,resnor
      CHARACTER*19 magstri,tmpstr
      CHARACTER*30 rastri,rdstri
      CHARACTER*37 rstri
      CHARACTER*3 radtyp
      CHARACTER*1 obstr,signo
      LOGICAL radar
      INTEGER truncat,iobcur,iotr,iore
c fix for exadecimal station code      
      CHARACTER*3 obsstr

      CALL filopn(unit,file,'UNKNOWN')

c check data set first
      radar=.false.
      DO i=1,n
         IF(iobs(i)/1000.eq.2)THEN
            radar=.true.
         ELSEIF(iobs(i)/1000.ne.1)THEN
            WRITE(*,*)'wrirwg: obs.type ',iobs(i),' unknown, rec.no=',i
            STOP
         ENDIF
      ENDDO
c ========= SORT OBSERVATIONS BY TIME ==============
      call heapsort(tutm,n,is)
c ========= HANDLE ASTROMETRY OBSERVATIONS =========
c astrometry header
      IF(rmsres.GT.0.D0) WRITE(unit,110) comcha,rmsres
 110  FORMAT(A1,'RMS of orbit residuals = ',F8.3)

      WRITE(unit,120) comcha
 120  FORMAT(A1,'++ OBJ ++ OBS +++++ DATE ++++++',
     +     '  +++++ RIGHT ASCENSION ++++++',
     +     '  ++++++ DECLINATION ++++++++',
     +     '  +++++ APP MAG +++++',
     +     '  ++++ QUALITY +++')
      WRITE(unit,221) comcha
 221  FORMAT(A1,'+ DESIG + TYP YYYY MM DD.dddddd',
     +     '  HH MM SS.sss   rms     resid',
     +     '  DD MM SS.ss   rms     resid',
     +     '  MAG COL rms   resid',
     +     '  OBS     CHI  SEL')
c ========= ASTROMETRY LOOP==========================
      DO   i=1,n
      IF(iobs(is(i))/1000.eq.1)THEN
         obstr=char(iobs(is(i))-1000)
c convert time
         CALL mjddat(tutm(is(i)),iday,imonth,iyear,hour)
         day=iday+hour/24.d0
c convert RA
         CALL sessag(alpha(is(i))*degrad/15.d0,signo,ihour,imin,sec)
         IF(signo.eq.'-')STOP 'wrirwg error: negative right ascension.'
c prepare RA string
         resnor=resa(is(i))*secrad*cos(delta(is(i)))
         WRITE(tmpstr,FMT='(F6.3)') sec
         CALL rmsp(tmpstr,ln)
         IF(ln.lt.6)tmpstr='0'//tmpstr
         IF(abs(resnor).gt.999.d0)THEN
            WRITE(rastri,131)ihour,imin,tmpstr, 
     +           rmsa(is(i))*secrad,resnor
 131        FORMAT(2x,I2.2,1x,I2.2,1x,A6,F7.2,1P,E9.1)
         ELSE
            WRITE(rastri,130)ihour,imin,tmpstr,
     +           rmsa(is(i))*secrad,resnor
 130        FORMAT(2x,I2.2,1x,I2.2,1x,A6,F7.2,F9.3)
         ENDIF
c convert DEC
         CALL sessag(delta(is(i))*degrad,signo,ideg,imin,sec)
c prepare DEC string
         WRITE(tmpstr,FMT='(F5.2)') sec
         CALL rmsp(tmpstr,ln)
         IF(ln.lt.5)tmpstr='0'//tmpstr
         IF(abs(resd(is(i))*secrad).gt.999.d0)THEN
            WRITE(rdstri,171)signo,ideg,imin,tmpstr,
     +           rmsd(is(i))*secrad,resd(is(i))*secrad
 171        FORMAT(1x,A1,I2.2,1x,I2.2,1x,A5,F7.2,1P,E9.1)
         ELSE
            WRITE(rdstri,170)signo,ideg,imin,tmpstr,
     +           rmsd(is(i))*secrad,resd(is(i))*secrad
 170        FORMAT(1x,A1,I2.2,1x,I2.2,1x,A5,F7.2,F9.3)
         ENDIF
c prepare MAG string
         magstri=' '
         IF(rmsmag(is(i)).lt.0)THEN
            WRITE(magstri,121)smag(is(i))
 121        FORMAT(a6,13x)
         ELSE
            WRITE(magstri,122)smag(is(i)),rmsmag(is(i))
 122        FORMAT(a6,1x,f5.2,7x)
         ENDIF
c output  astrometry   
         WRITE(tmpstr,FMT='(F9.6)') day
         CALL rmsp(tmpstr,ln)
         IF(ln.lt.9)tmpstr='0'//tmpstr
c fix for exadecimal station code      
         CALL codestat(obscod(is(i)),obsstr)
c write on .rwo file
         WRITE(unit,101) objid(is(i)),obstr,iyear,imonth,tmpstr,
     +        rastri,rdstri,magstri,
     +        obsstr,sel(is(i))
 101     FORMAT(1x,A9,2x,a1,2x,I4,1x,I2.2,1x,A9,
     +        A30,A29,2x,A19,
     +        2X,A3,9X,2X,I1,3x)
      ENDIF
      ENDDO
c ========= HANDLE RADAR OBSERVATIONS ==========
      If(.not.radar) GOTO 99
c radar header      
      WRITE(unit,128) comcha
 128  FORMAT(A1,'++ OBJ ++ OBS ++++++ DATE +++++++ ',
     +     '++++++++ RADAR RANGE/RANGE RATE +++++++++ ',
     +     '++++++ QUALITY +++++')
      WRITE(unit,228) comcha
 228  FORMAT(A1,'+ DESIG + TYP YYYY MM DD hh:mm:ss ',
     +     'TYP   KM or KM/DAY  a priori rms residual ',
     +     'TRX REC     CHI  SEL')
c ========= RADAR LOOP==========================
      DO i=1,n
         IF(iobs(is(i))/1000.eq.2)THEN
            IF(iobs(is(i))-2000.ge.100)THEN
c surface bounce
               iobcur=iobs(is(i))-2100
               obstr='r'
            ELSE
c already corrected to center of mass
               iobcur=iobs(is(i))-2000
               obstr='R'
            ENDIF
c convert time
            CALL mjddat(tutm(is(i)),iday,imonth,iyear,hour)
c convert hour to 12:12:12
            ihour=truncat(hour,1d-7)
            minu=(hour-ihour)*60.d0
            imin=truncat(minu,1d-5)
            sec=(minu-imin)*60.d0
            isec=truncat(sec,1d-3)
            IF(iobcur.eq.2)THEN
c range observation
               radtyp='DEL'               
               IF(abs(resa(is(i))*au).gt.999.d0)THEN
                  WRITE(rstri,143)alpha(is(i))*au,rmsa(is(i))*au,
     +                 resa(is(i))*au
               ELSE
                  WRITE(rstri,141)alpha(is(i))*au,rmsa(is(i))*au,
     +                 resa(is(i))*au
               ENDIF
            ELSEIF(iobcur.eq.3)THEN
c range-rate observation
               radtyp='DOP'               
               IF(abs(resd(is(i))*au).gt.999.d0)THEN
                  WRITE(rstri,143)delta(is(i))*au,rmsd(is(i))*au,
     +                 resd(is(i))*au
               ELSE
                  WRITE(rstri,141)delta(is(i))*au,rmsd(is(i))*au,
     +                 resd(is(i))*au
               ENDIF
            ELSE
               STOP '*** wrirwg.f: internal error(1) ***'
            ENDIF
 143        FORMAT(f16.5,1x,f9.5,1x,1p,e10.4)
 141        FORMAT(f16.5,1x,f9.5,1x,f10.5)
c find codes of two observatories
            iotr= obscod(is(i))/10000
            iore= obscod(is(i))-iotr*10000
c output radar data              
            WRITE(unit,102) objid(is(i)),obstr,iyear,imonth,iday,
     +           ihour,imin,isec,radtyp,rstri,
     +           iotr,iore,sel(is(i))
 102        FORMAT(1x,A9,2x,a1,2x,I4,1x,I2.2,1x,I2.2,1x,
     +           I2.2,':',I2.2,':',I2.2,1x,A3,1x,A37,
     +           1X,I3.3,1X,I3.3, 9X, 2X,I1,3x)
         ENDIF
      ENDDO
c ============================================
 99   CALL filclo(unit,' ')
      END
* Copyright (C) 1997-1999 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: November 26, 1999
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         A S T R O W                           *
*  *                                                               *
*  *     Computation of a-priori RMS of optical observations       *
*  *         from results of a statistical analysis or             *
*  *    from their accuracy (= number of digits in input file)     *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    MPCTYP    -  Observation type (column 15 of MPC record)
*           TDT       -  Time of observation (MJD, TDT)
*           IDSTA     -  Observatory code
*           ACCA      -  Accuracy of right ascension (rad)
*           ACCD      -  Accuracy of declination (rad)
*
* OUTPUT:   RMSA      -  A-priori RMS of right ascension (rad)
*           RMSD      -  A-priori RMS of declination (rad)
*
      SUBROUTINE astrow(mpctyp,tdt,idsta,acca,accd,rmsa,rmsd)
      IMPLICIT NONE

      INTEGER idsta
      DOUBLE PRECISION tdt,acca,accd,rmsa,rmsd
      CHARACTER*(*) mpctyp

      INCLUDE 'trig.h'

      INTEGER le
      DOUBLE PRECISION rmsmin
      CHARACTER*5 ads(2)
      CHARACTER*80 filea,filed,ermnam
      LOGICAL ermuse,fail,fail1,found,error

      SAVE ermuse

      LOGICAL first
      DATA first/.true./
      SAVE first

* Additional info on how RMS are obtained
      INCLUDE 'astrow.h'

      INTEGER lench
      EXTERNAL lench

* Input of RMS class definitions
      IF(first) THEN
          fail=.false.
          ermuse=.false.
          CALL rdnlog('errmod.','use',ermuse,.false.,found,
     +                fail1,fail)
          IF(ermuse) THEN
              ermnam='num'
              CALL rdncha('errmod.','name',ermnam,.false.,found,
     +                    fail1,fail)
              IF(fail) STOP '**** astrow: abnormal END ****'
              le=lench(ermnam)
              filea=ermnam(1:le)//'.cla'
              filed=ermnam(1:le)//'.cld'
              CALL rrmscl(filea,filed,.false.)
          END IF
          IF(fail) STOP '**** astrow: abnormal END ****'
          first=.false.
      END IF

* Negative values means: not yet assigned
      rmsa=-1.D0
      rmsd=-1.D0
      orstep(1)=0
      orstep(2)=0

* STEPs 1/2: look for a specific RMS class for the observatory
      IF(ermuse) THEN
          CALL accstr(acca,accd,ads(1),ads(2),error)
          IF(error) GOTO 2
          CALL crmscl(idsta,ads,mpctyp,tdt,rmsa,rmsd,idcl,orstep,tdtlim)
          IF(orstep(1).GT.0) rmsa=rmsa*radsec
          IF(orstep(2).GT.0) rmsd=rmsd*radsec
      END IF

* STEP 3: default, rough rule-of-thumb error model
 2    CONTINUE
      IF(rmsa.LT.0.D0 .OR. rmsd.LT.0.D0) THEN
* Minimum value of RMS (time dependent)
* Before 1890
          IF(tdt.LT.11368.d0) THEN
              rmsmin=3.d0
* From 1890 to 1950
          ELSE IF(tdt.LT.33282.d0) THEN
              rmsmin=2.d0
* After 1950
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

      END
************************************************************************
* 'magrms' returns the default magnitude rms based on the MPC obs string
************************************************************************
      DOUBLE PRECISION FUNCTION magrms(magstr,tdt,idsta,typ)
      IMPLICIT NONE
      
      CHARACTER*6 magstr
c obs. type (from column 15 of .obs format), station code, time MJD
      CHARACTER*1 typ
      INTEGER idsta
      DOUBLE PRECISION tdt
      INTEGER ll,lench
c magnitude weighting for now is simple; 
c should be based on digits given,color
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
      END
* Copyright (C) 1999-2000 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: November 28, 2000
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         R R M S C L                           *
*  *                                                               *
*  *            Read from a file and store RMS classes             *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    FILRA     -  File name (RA classes)
*           FILDEC    -  File name (DEC classes)
*           CLONLY    -  Load only class list, not accuracy description
*                        (this option is enabled when the routine is
*                        used by stat2.f)
*
* WARNING: the correct behaviour of the routine depends strictly on the
* sorting of input data. More precisely:
*   - "per-observatory" RMS classes must be sorted in increasing order
*      of observatory code (field 3);
*   -  "mixed" (all observatories together) classes must be sorted
*      according to the following sorting keys:
*         key 1: accuracy class (field 1)
*         key 2: type of observation (field 2)
*         key 3: starting time (field 4)
*
      SUBROUTINE rrmscl(filra,fildec,clonly)
      IMPLICIT NONE

      CHARACTER*(*) filra,fildec
      LOGICAL clonly

      INCLUDE 'parobc.h'
      INCLUDE 'parrms.h'

* Common blocks to be initialized:
      INCLUDE 'rmscl.h'

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

* Check file format (in order to avoid using files written according
* to old format, not containing bias information)
      READ(unit,102) rec
 102  FORMAT(A)
      IF(lench(rec).LT.110) THEN
          WRITE(*,200) file(1:lf)
          STOP '**** rrmscl: abnormal end ****'
      END IF
 200  FORMAT('**** ERROR: file "',A,'" is written in an obsolete ',
     +       'format: please use a more recente version ****')
      REWIND(unit)

* Initializations (section 1)
      ncrm(ic)=0
      DO 1 i=obsc1,obsc2
      crmobp(1,i,ic)=0
      crmobp(2,i,ic)=0
 1    CONTINUE
      nin=0
      n=1
      obscp=obsc1-999

* Initializations (section 2)
      n1=0
      n2=0
      n3=0

* Start reading loop
 2    CONTINUE
      IF(clonly) THEN
          READ(unit,100,END=3) crmad1,crmod1,cobsc,mjde1,mjde2
          rms1=0
          ave1=0
          rms1a=0
      ELSE
          READ(unit,100,END=3) crmad1,crmod1,cobsc,mjde1,mjde2,
     +                         rms1,ave1,rms1a
      END IF
 100  FORMAT(A5,1X,A1,1X,A3,1X,F8.1,1X,F8.1,1X,F9.3,F11.5,F9.3,
     +       E11.3,F9.3,2I8,F7.2,I9)
      nin=nin+1
      IF(crmod1.EQ.' ') crmod1='P'

      mixobs=(cobsc.EQ.'ALL')

      IF(mixobs) THEN
          n3=n3+1
          IF(n3.GT.crx3nx) STOP '**** rrmscl: n3 > crx3nx ****'

* Check for change in accuracy descriptor
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

* Check for change in observation type
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
          READ(cobsc,101) obsc
          IF(obsc.LT.obsc1 .OR. obsc.GT.obsc2)
     +        STOP '**** rrmscl: input error (01) ****'
          crmad(n,ic)=crmad1
          crmotd(n,ic)=crmod1

* Integer codification of type descriptor
          crmoti(n,ic)=1000+ICHAR(crmotd(n,ic))
          crmot1(n,ic)=mjde1
          crmot2(n,ic)=mjde2
          crmrms(n,ic)=rms1
          crmave(n,ic)=ave1
          crmrma(n,ic)=rms1a

          IF(obsc.NE.obscp) THEN
              IF(crmobp(1,obsc,ic).NE.0)
     +            STOP '**** rrmscl: input error (02) ****'
              crmobp(1,obsc,ic)=n
              obscp=obsc
          END IF
          crmobp(2,obsc,ic)=n
          n=n+1
      END IF
 101  FORMAT(I3)

      GOTO 2

 3    CONTINUE
      ncrm(ic)=n-1
      crx1n(ic)=n1
      crx2n(ic)=n2
      crx3n(ic)=n3
      CALL filclo(unit,' ')
 10   CONTINUE

      iiccrm=36

      END
* Copyright (C) 1999 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: March 3, 1999
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         A C C S T R                           *
*  *                                                               *
*  *               Accuracy description (string)                   *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    ACCA      -  Accuracy of right ascension (rad)
*           ACCD      -  Accuracy of declination (rad)
*
* OUTPUT:   ADSA      -  Accuracy description string (RA)
*           ADSD      -  Accuracy description string (DEC)
*           ERROR     -  Error flag
*
      SUBROUTINE accstr(acca,accd,adsa,adsd,error)
      IMPLICIT NONE

      DOUBLE PRECISION acca,accd
      CHARACTER*(*) adsa,adsd
      LOGICAL error

      INCLUDE 'trig.h'

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
 101  FORMAT(F10.4)
      CALL rmsp(ads,lads)
      DO 2 k=lads,1,-1
      IF(ads(k:k).EQ.'0') THEN
          ads(k:k)=' '
          lads=lads-1
      ELSE
          GOTO 3
      END IF
 2    CONTINUE
 3    CONTINUE
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
 1    CONTINUE
 200  FORMAT('ERROR (accstr): ',A,' = "',A,'"')

      END
* Copyright (C) 1999 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: November 24, 1999
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         C R M S C L                           *
*  *                                                               *
*  *   Computes a-priori observation RMS based on known classes    *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    OBSCOD    -  Observatory code
*           ADS       -  Accuracy description string (RA/DEC)
*           MPCTYP    -  Observation type (column 15 of MPC record)
*           TDT       -  Time (MJD, TDT)
*
* OUTPUT:   RMSA      -  A-priori RMS of RA (arcsec)
*           RMSD      -  A-priori RMS of DEC (arcsec)
*           IDCL      -  Class progressive number in index
*           STEP      -  RMS assignation steps:
*                           0 = not assigned
*                           1 = specific station class
*                           2 = mixed class
*           TDTLIM    -  Class limits (MJD, TDT)
*
* WARNING: if no valid class for the observation is found, the
*          corresponding output rms is set to a negative value
*
      SUBROUTINE crmscl(obscod,ads,mpctyp,tdt,rmsa,rmsd,idcl,step,
     +                  tdtlim)
      IMPLICIT NONE

      INTEGER obscod,idcl(2),step(2)
      DOUBLE PRECISION rmsa,rmsd,tdt,tdtlim(2,2)
      CHARACTER*(*) mpctyp,ads(2)

      INCLUDE 'parobc.h'
      INCLUDE 'parrms.h'

* NEEDED common blocks:
      INCLUDE 'rmscl.h'

      INTEGER ic,i,i1,i2,i3
      DOUBLE PRECISION rms1

      IF(iiccrm.NE.36) STOP '**** crmscl: internal error (01) ****'
      IF(obscod.LT.obsc1 .OR. obscod.GT.obsc2)
     +    STOP '**** crmscl: input error (01) ****'

      DO 20 ic=1,2
      rms1=-1.d0
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
          idcl(ic)=i
          step(ic)=1
          tdtlim(1,ic)=crmot1(i,ic)
          tdtlim(2,ic)=crmot2(i,ic)
          GOTO 2
 1        CONTINUE
 2        CONTINUE
      END IF
      IF(step(ic).GT.0) GOTO 10
      DO 3 i1=1,crx1n(ic)
      IF(ads(ic).EQ.crx1ad(i1,ic)) THEN
          DO 4 i2=crx1pt(1,i1,ic),crx1pt(2,i1,ic)
          IF(mpctyp.EQ.crx2ty(i2,ic)) THEN
              DO 5 i3=crx2pt(1,i2,ic),crx2pt(2,i2,ic)
              IF(tdt.GE.crx3t(1,i3,ic) .AND. tdt.LT.crx3t(2,i3,ic)) THEN
                  rms1=crx3r(i3,ic)
                  idcl(ic)=i3
                  step(ic)=2
                  tdtlim(1,ic)=crx3t(1,i3,ic)
                  tdtlim(2,ic)=crx3t(2,i3,ic)
                  GOTO 10
              END IF
 5            CONTINUE
          END IF
 4        CONTINUE
      END IF
 3    CONTINUE

 10   CONTINUE
      IF(ic.EQ.1) THEN
          rmsa=rms1
      ELSE
          rmsd=rms1
      END IF
 20   CONTINUE

      END
* Copyright (C) 1997-2000 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: November 28, 2000
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         A S T R W B                           *
*  *                                                               *
*  *     Computation of a-priori RMS of optical observations       *
*  *         from results of a statistical analysis or             *
*  *    from their accuracy (= number of digits in input file)     *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    MPCTYP    -  Observation type (column 15 of MPC record)
*           TDT       -  Time of observation (MJD, TDT)
*           IDSTA     -  Observatory code
*           ACCA      -  Accuracy of right ascension (rad)
*           ACCD      -  Accuracy of declination (rad)
*
* OUTPUT:   RMSA      -  A-priori RMS of right ascension (rad)
*           RMSD      -  A-priori RMS of declination (rad)
*           BIASA     -  Bias in right ascension (rad)
*           BIASD     -  Bias in declination (rad)
*           DSTEP     -  Decision step (RA/DEC):
*                           0 = ERROR
*                           1 = single-station class
*                           2 = multi-station class
*                           3 = default (1/2/3 rule)
*
      SUBROUTINE astrwb(mpctyp,tdt,idsta,acca,accd,rmsa,rmsd,
     +                  biasa,biasd,dstep)
      IMPLICIT NONE

      INTEGER idsta,dstep(2)
      DOUBLE PRECISION tdt,acca,accd,rmsa,rmsd,biasa,biasd
      CHARACTER*(*) mpctyp

      INCLUDE 'trig.h'

      INTEGER le
      DOUBLE PRECISION rmsmin,rmsap,rmsdp,avea,aved
      CHARACTER*5 ads(2)
      CHARACTER*80 filea,filed,ermnam
      LOGICAL ermuse,fail,fail1,found,error

      SAVE ermuse

      LOGICAL first
      DATA first/.true./
      SAVE first

* Additional info on how RMS are obtained
      INCLUDE 'astrow.h'

      INTEGER lench
      EXTERNAL lench

* Input of RMS class definitions
      IF(first) THEN
          fail=.false.
          ermuse=.false.
          CALL rdnlog('errmod.','use',ermuse,.false.,found,
     +                fail1,fail)
          IF(ermuse) THEN
              ermnam='num'
              CALL rdncha('errmod.','name',ermnam,.false.,found,
     +                    fail1,fail)
              IF(fail) STOP '**** astrwb: abnormal END ****'
              le=lench(ermnam)
              filea=ermnam(1:le)//'.cla'
              filed=ermnam(1:le)//'.cld'
              WRITE(*,210) ermnam(1:le)
              CALL rrmscl(filea,filed,.false.)
          ELSE
              WRITE(*,211)
          END IF
          IF(fail) STOP '**** astrwb: abnormal END ****'
          first=.false.
      END IF
 210  FORMAT('INFO(astrwb): reading error RMS model ',
     +       'from files "',A,'.cl[ad]"')
 211  FORMAT('INFO(astrwb): NO specific error RMS model is used')

* Negative values means: not yet assigned
      rmsa=-1.D0
      rmsd=-1.D0
      biasa=0.D0
      biasd=0.D0
      dstep(1)=0
      dstep(2)=0

* STEPs 1/2: look for a specific RMS class for the observatory
      IF(ermuse) THEN
          CALL accstr(acca,accd,ads(1),ads(2),error)
          IF(error) GOTO 2
          CALL crmscn(idsta,ads,mpctyp,tdt,rmsap,rmsdp,
     +                avea,aved,rmsa,rmsd,idcl,dstep,tdtlim)
          IF(dstep(1).GT.0) THEN
              rmsa=rmsa*radsec
              biasa=avea*radsec
          END IF
          IF(dstep(2).GT.0) THEN
              rmsd=rmsd*radsec
              biasd=aved*radsec
          END IF
      END IF

* STEP 3: default, rough rule-of-thumb error model
 2    CONTINUE
      IF(rmsa.LT.0.D0 .OR. rmsd.LT.0.D0) THEN
* Minimum value of RMS (time dependent)
* Before 1890
          IF(tdt.LT.11368.d0) THEN
              rmsmin=3.d0
* From 1890 to 1950
          ELSE IF(tdt.LT.33282.d0) THEN
              rmsmin=2.d0
* After 1950
          ELSE
              rmsmin=1.d0
          END IF
          rmsmin=rmsmin*radsec
          IF(rmsa.LT.0.D0) THEN
              rmsa=MAX(rmsmin,acca)
              dstep(1)=3
          END IF
          IF(rmsd.LT.0.D0) THEN
              rmsd=MAX(rmsmin,accd)
              dstep(2)=3
          END IF
      END IF
      orstep(1)=dstep(1)
      orstep(2)=dstep(2)

      END
* Copyright (C) 1999-2000 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: October 31, 2000
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         C R M S C N                           *
*  *                                                               *
*  *   Computes a-priori observation RMS based on known classes    *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    OBSCOD    -  Observatory code
*           ADS       -  Accuracy description string (RA/DEC)
*           MPCTYP    -  Observation type (column 15 of MPC record)
*           TDT       -  Time (MJD, TDT)
*
* OUTPUT:   RMSA      -  A-priori RMS of RA (arcsec) before sub. ave.
*           RMSD      -  A-priori RMS of DEC (arcsec) before sub. ave.
*           AVEA      -  Average residual (bias) in RA (arcsec)
*           AVED      -  Average residual (bias) in DEC (arcsec)
*           RMSAA     -  A-priori RMS of RA (arcsec) after sub. ave.
*           RMSDA     -  A-priori RMS of DEC (arcsec) after sub. ave.
*           IDCL      -  Class progressive number in index
*           STEP      -  RMS assignation steps:
*                           0 = not assigned
*                           1 = specific station class
*                           2 = mixed class
*           TDTLIM    -  Class limits (MJD, TDT)
*
* WARNING: if no valid class for the observation is found, the
*          corresponding output rms is set to a negative value
*
      SUBROUTINE crmscn(obscod,ads,mpctyp,tdt,rmsa,rmsd,avea,aved,
     +                  rmsaa,rmsda,idcl,step,tdtlim)
      IMPLICIT NONE

      INTEGER obscod,idcl(2),step(2)
      DOUBLE PRECISION rmsa,rmsd,avea,aved,rmsaa,rmsda,tdt,tdtlim(2,2)
      CHARACTER*(*) mpctyp,ads(2)

      INCLUDE 'parobc.h'
      INCLUDE 'parrms.h'

* NEEDED common blocks:
      INCLUDE 'rmscl.h'

      INTEGER ic,i,i1,i2,i3
      DOUBLE PRECISION rms1,ave1,rms1a

      IF(iiccrm.NE.36) STOP '**** crmscn: internal error (01) ****'
      IF(obscod.LT.obsc1 .OR. obscod.GT.obsc2)
     +    STOP '**** crmscn: input error (01) ****'

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
 1        CONTINUE
 2        CONTINUE
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
 5            CONTINUE
          END IF
 4        CONTINUE
      END IF
 3    CONTINUE

 10   CONTINUE
      IF(ic.EQ.1) THEN
          rmsa=rms1
          avea=ave1
          rmsaa=rms1a
      ELSE
          rmsd=rms1
          aved=ave1
          rmsda=rms1a
      END IF
 20   CONTINUE

      END
* Copyright (C) 2000 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: November 28, 2000
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         R D C O R M                           *
*  *                                                               *
*  *   Read the autocorrelation models for all the observatories   *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    FILE      -  Input file
*
* OUTPUT:   
*
      SUBROUTINE rdcorm(file)
      IMPLICIT NONE

      CHARACTER*(*) file

      INCLUDE 'parobc.h'
      INCLUDE 'parcor.h'
* Common blocks to be initialized:
* Autocorrelation models
      INCLUDE 'comcor.h'

      INCLUDE 'trig.h'

      INTEGER unit,io,line,ic,lf,ip,ivers
      CHARACTER rec*80

      INTEGER lench
      EXTERNAL lench

* Initialization
      DO 1 io=obsc1,obsc2
      pto2f(io,1,1)=0
      pto2f(io,2,1)=0
 1    CONTINUE
      nfunt=0
      npart=0
      line=0

* Reading model
      lf=lench(file)
      CALL filopl(unit,file)
      READ(unit,*) ivers
      IF(ivers.NE.1) THEN
          WRITE(*,201) file(1:lf)
          STOP '**** Abnormal end ****'
      END IF
 201  FORMAT('ERROR(rdcorm): unsupported version of file "',A,'"')
      READ(unit,*) aprmx
      READ(unit,*) maxdst
      minapw=1.D0/((aprmx*radsec)**2)

 2    CONTINUE
* Read a new observatory/coordinate
      READ(unit,100,END=10) rec
 100  FORMAT(A)
      line=line+1
* IC = coordinate (1=RA, 2=DEC)
      IF(rec(1:8).EQ.'TIME RA ') THEN
          ic=1
      ELSEIF(rec(1:8).EQ.'TIME DEC') THEN
          ic=2
      ELSE
          GOTO 20
      END IF
* IO = observatory code
      READ(rec(9:),*,ERR=20) io
      IF(io.LT.obsc1) STOP '**** rdcorm: obsc < obsc1 ****'
      IF(io.GT.obsc2) STOP '**** rdcorm: obsc > obsc2 ****'
      pto2f(io,ic,1)=nfunt+1
      nfo(io,ic)=0
      nparo(io,ic)=0

 3    CONTINUE
* Read the model for one coordinate (IC) of one observatory (IO)
      READ(unit,100,ERR=20) rec
      line=line+1
      IF(rec.EQ.'END') THEN
          IF(nfo(io,ic).LE.0) GOTO 20
          GOTO 2
      END IF
      nfunt=nfunt+1
      IF(nfunt.GT.nparx) STOP '**** rdcorm: nfunt > nparx ****'
      pto2f(io,ic,2)=nfunt
      nfo(io,ic)=nfo(io,ic)+1
* Determine function integer code (KP1) and number of parameters (NP1)
      CALL fcsfun(rec,kfun(nfunt),nparf(nfunt))
      IF(kfun(nfunt).LE.0) GOTO 20
      nparo(io,ic)=nparo(io,ic)+nparf(nfunt)
* Reading and storing parameter values and properties
      IF(npart+nparf(nfunt).GT.nparx)
     +    STOP '**** rdcorm: npart > nparx ****'
      ptf2p(1,nfunt)=npart+1
      DO 4 ip=1,nparf(nfunt)
      npart=npart+1
      READ(unit,*,ERR=20) par(npart)
 4    CONTINUE
      ptf2p(2,nfunt)=npart
      GOTO 3

* Regular end
 10   CONTINUE
      CALL filclo(unit,' ')
      iiccor=36
      RETURN

* Error termination
 20   CONTINUE
      WRITE(*,200) file(1:lf),line
 200  FORMAT('rdcorm: INPUT ERROR from file "',A,'" at record',I5)
      STOP '**** rdcorm: abnormal end ****'

      END
* Copyright (C) 2000 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: November 28, 2000
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         O B S C O R                           *
*  *                                                               *
*  * Correlation between two observations of the same observatory  *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    OBSCOD    -  Observatory code
*           T1        -  Time (MJD)                         \
*           ALPHA1    -  Right ascension (rad)               |
*           DEC1      -  Declination (rad)                   |  OBS 1
*           WRA1      -  A-priori weight for RA (rad**(-2))  |
*           WDEC1     -  A-priori weight for DEC (rad**(-2)) |
*           IDST1     -  Decision step (RA/DEC)             /
*           T2        -  Time (MJD)                         \
*           ALPHA2    -  Right ascension (rad)               |
*           DEC2      -  Declination (rad)                   |  OBS 2
*           WRA2      -  A-priori weight for RA (rad**(-2))  |
*           WDEC2     -  A-priori weight for DEC (rad**(-2)) |
*           IDST2     -  Decision step (RA/DEC)             /
*
* OUTPUT:   COVRA     -  Off-diagonal covariance in RA (rad**2)
*           COVDEC    -  Off-diagonal covariance in DEC (rad**2)
*
      SUBROUTINE obscor(obscod,t1,alpha1,dec1,wra1,wdec1,idst1,
     +                  t2,alpha2,dec2,wra2,wdec2,idst2,covra,covdec)
      IMPLICIT NONE

      INTEGER obscod,idst1(2),idst2(2)
      DOUBLE PRECISION t1,t2,alpha1,alpha2,dec1,dec2,wra1,wra2
      DOUBLE PRECISION wdec1,wdec2,covra,covdec

      INCLUDE 'parcor.h'
      INCLUDE 'parobc.h'
      INCLUDE 'comcor.h'

      INTEGER ic,ip1,if1,lf
      DOUBLE PRECISION dt
      CHARACTER*80 file
      LOGICAL first,fail,fail1,found,doit

      DATA first/.true./
      SAVE first,doit

      INTEGER lench
      EXTERNAL lench

* Initialization
      IF(first) THEN
          first=.false.
          fail=.false.
* Flag telling whether time correlations are to be computed at all
          doit=.false.
          CALL rdnlog('errmod.obscor.','use',doit,.false.,found,
     +                fail1,fail)
* Name of input file
          IF(doit) THEN
              file='obscor.mod'
              CALL rdncha('errmod.obscor.','file',file,.false.,found,
     +                    fail1,fail)
          END IF
          IF(fail) STOP '**** obscor: abnormal end ****'

* Input of correlation model
          IF(doit) THEN
              lf=lench(file)
              WRITE(*,210) file(1:lf)
              CALL rdcorm(file)
          ELSE
              WRITE(*,211)
          END IF
      END IF
 210  FORMAT('INFO(obscor): reading error correlation model ',
     +       'from file "',A,'"')
 211  FORMAT('INFO(obscor): NO error correlation model is used')

      covra=0
      covdec=0
      IF(.NOT.doit) RETURN

      IF(iiccor.NE.36) STOP '**** obscor: internal error (01) ****'
      IF(obscod.LT.obsc1.OR.obscod.GT.obsc2)
     +    STOP '**** obscor: internal error (02) ****'
      dt=ABS(t2-t1)

* Right ascension
      ic=1
      IF(idst1(ic).EQ.0 .OR. idst2(ic).EQ.0) GOTO 1
      IF(idst1(ic).GT.maxdst) GOTO 1
      IF(idst2(ic).GT.maxdst) GOTO 1
      IF(wra1.LT.minapw) GOTO 1
      IF(wra2.LT.minapw) GOTO 1
      if1=pto2f(obscod,ic,1)
      IF(if1.LE.0) GOTO 1
      ip1=ptf2p(1,if1)
      CALL fcorob(kfun(if1),nfo(obscod,ic),par(ip1),nparf(if1),dt,
     +            covra)
      covra=covra/SQRT(wra1*wra2)
 1    CONTINUE

* Declination
      ic=2
      IF(idst1(ic).EQ.0 .OR. idst2(ic).EQ.0) GOTO 2
      IF(idst1(ic).GT.maxdst) GOTO 2
      IF(idst2(ic).GT.maxdst) GOTO 2
      IF(wdec1.LT.minapw) GOTO 2
      IF(wdec2.LT.minapw) GOTO 2
      if1=pto2f(obscod,ic,1)
      IF(if1.LE.0) GOTO 2
      ip1=ptf2p(1,if1)
      CALL fcorob(kfun(if1),nfo(obscod,ic),par(ip1),nparf(if1),dt,
     +            covdec)
      covdec=covdec/SQRT(wdec1*wdec2)
 2    CONTINUE

      END
* Copyright (C) 2000 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: November 16, 2000
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         F C O R O B                           *
*  *                                                               *
*  *              Computation of correlation function              *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    KFUN      -  Function integer codes
*           NFUN      -  Number of elementary functions
*           PAR       -  Parameter values
*           NPARF     -  Number of parameters per function
*           DT        -  Time lag
*
* OUTPUT:   FUN       -  Function value
*
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

      INCLUDE 'fcfund.h'

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

 1    CONTINUE
      fun=fun+term

      END
* Copyright (C) 2000 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: November 3, 2000
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         F C S F U N                           *
*  *                                                               *
*  *              LS fit of covariance functions:                  *
*  *       function integer codes and number of parameters         *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    NAME      -  Function name
*
* OUTPUT:   KFUN      -  Function integer identificator (0=don't know)
*           NPAR      -  Number of parameters
*
      SUBROUTINE fcsfun(name,kfun,npar)
      IMPLICIT NONE

      INTEGER kfun,npar
      CHARACTER*(*) name

      kfun=0
      npar=0
* Multiplicative coefficient (constant)
      IF(name.EQ.'+COEF') THEN
          kfun=1
          npar=1
* Exponential function
      ELSEIF(name.EQ.'*EXP') THEN
          kfun=2
          npar=1
* Normal function
      ELSEIF(name.EQ.'*NORM') THEN
          kfun=3
          npar=1
* Parabola (to be multiplied by exponential or normal functions)
      ELSEIF(name.EQ.'*PARAB') THEN
          kfun=4
          npar=1
* ADD HERE NEW FUNCTIONS, using increasing integer identificator (KFUN)
* Remember to update accordingly also subroutines FCFUND and FCWPAR
      END IF

      END
* Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: October 31, 1997
* version 1.8.0, Steven Chesley, Dec. 14, 1998
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                          M P C I N                            *
*  *                                                               *
*  *   Input of astrometric+ radar  observations from  MPC format  *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    FILE      -  File name
*           NOBSX     -  Max dimension of vectors (iobs,t,alpha,delta,
*                        obscod,acct,acca,accd)
*
* OUTPUT:   IOBS      -  Observation type: 1=astrometry: alpha,delta 
*                                          2,3,4=radar: range,range-rate
*           TDT       -  Time (MJD, TDT)
*           TUTM      -  Time (MJD, UTM)
*           ALPHA     -  Right ascension (rad)-Range (AU) if radar
*           DELTA     -  Declination (rad)-Range-rate (AU/day) if radar
*           OBSCOD    -  Observatory code
*           ACCT      -  Accuracy of time (d)
*           ACCA      -  Accuracy of right ascension (rad) - of range (AU)
*           ACCD      -  Accuracy of declination (rad) - of range-rate (AU/day)
*           SMAG      -  Observed magnitude (string, incl. color)
*           NOBS      -  Number of observations read
c           OBS       -  successful input flag
*
      subroutine mpcin(obs,file,objid,iobs,tdt,tutm,alpha,delta,obscod,
     +                 acct,acca,accd,smag,nobs,nobsx)
      implicit none

      LOGICAL obs

      character*(*) file
      integer nobsx
      character*(*) objid(nobsx)
      integer iobs(nobsx)
      character*6 smag(nobsx)
      double precision tdt(nobsx),tutm(nobsx),alpha(nobsx),delta(nobsx)
      double precision acct(nobsx),acca(nobsx),accd(nobsx)
      integer nobs,obscod(nobsx)
c ==========================================
      logical error
      character*80 rec1,rec2
      character*1 obscha
      integer l,unit
      external lench
      integer lench
c use MPC radar obs if radmpc=.true.
      logical radmpc
      data radmpc/.false./

      call filopn(unit,file,'old')

      nobs=0
 1    continue
      read(unit,101,end=10) rec1
 101  format(a)
      l=lench(rec1)
      IF(l.eq.0)GOTO 1
      nobs=nobs+1
      if(nobs.gt.nobsx) stop ' **** mpcin: nobs > nobsx ****'
c ================= HANDLE RADAR ============================
      if(rec1(15:15).eq.'R'.or.rec1(15:15).eq.'r')then
c radmpc determines whether MPC radar observations are used or ignored.
         if(radmpc)then
            write(*,*)'Radar observation near record ',nobs
            read(unit,101,end=20) rec2
c  alpha is actually r, delta is rdot
            call mpcrad(rec1,rec2,iobs(nobs),
     +           objid(nobs),tdt(nobs),tutm(nobs),alpha(nobs),
     +           delta(nobs),obscod(nobs),acct(nobs),acca(nobs),
     +           accd(nobs),error)
            smag(nobs)='      '
c  generate range and range rate only if mixed observation
            if(mod(iobs(nobs),100).eq.1)then
               iobs(nobs)=iobs(nobs)+1
               iobs(nobs+1)=iobs(nobs)+1
               objid(nobs+1)=objid(nobs)
               tdt(nobs+1)=tdt(nobs)
               tutm(nobs+1)=tutm(nobs)
               alpha(nobs+1)=0d0
               delta(nobs+1)=delta(nobs)
               delta(nobs)=0d0
               obscod(nobs+1)=obscod(nobs)
               acct(nobs+1)=acct(nobs)
               acca(nobs+1)=-1d0
               accd(nobs+1)=accd(nobs)
               accd(nobs)=-1d0
               smag(nobs+1)='      '
               nobs=nobs+1
            endif
         else
            if(rec1(15:15).eq.'R')write(*,*) 'ignoring MPC radar obs'
            nobs=nobs-1
            goto 1
         endif
c ================= HANDLE SATELLITE ============================
      elseif(rec1(15:15).eq.'S'.or.rec1(15:15).eq.'s')then
         write(*,*)'Satellite observation at record ',nobs
         nobs=nobs-1
c         iobs(nobs)=3000+j
         goto 1
c ================= HANDLE OPTICAL ============================
      else
         call mpctr(rec1,objid(nobs),tdt(nobs),tutm(nobs),alpha(nobs),
     +        delta(nobs),obscod(nobs),acct(nobs),acca(nobs),
     +        accd(nobs),smag(nobs),error)
c assign iobs=1000+ichar(obscha), but replace blanks with P's
         obscha=rec1(15:15)
         if(obscha.eq.' ')obscha='P'
         if(.not.error)iobs(nobs)=1000+ichar(obscha)
c handle in a temporary way the roving observer problem
         IF(obscod(nobs).eq.247)THEN
            nobs=nobs-1
            GOTO 1
         ENDIF
      endif
      if(error) then
          l=lench(file)
          write(*,102) file(1:l),nobs
 102      format(' **** mpcin: input conversion error ****'/
     .       ' **** file = "',a,'", line',i5,' ****')
          nobs=0
          obs=.false.
          call filclo(unit,' ')
          RETURN
      end if
      goto 1
c regular ending
 10   continue
      call filclo(unit,' ')
      IF(nobs.eq.0)THEN
         obs=.false.
      ELSE
         obs=.true.
      ENDIF
      return
c erroneous radar ending
 20   write(*,*) 'mpcin: incomplete radar observation, end of file'
      stop
      end
* Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: March 12, 1997
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                          M P C T R                            *
*  *                                                               *
*  *   Transformation of an astrometric observation (MPC format)   *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    REC       -  MPC record
*
* OUTPUT:   OBJID     -  IAU Object designation
*           TDT       -  Time (MJD, TDT)
*           TUTM      -  Time (MJD, UTM)
*           ALPHA     -  Right ascension (rad)
*           DELTA     -  Declination (rad)
*           OBSCOD    -  Observatory code
*           ACCT      -  Accuracy of time (d)
*           ACCA      -  Accuracy of right ascension (rad)
*           ACCD      -  Accuracy of declination (rad)
*           SMAG      -  Apparent magnitude and color (string)
*           ERROR     -  Conversion error
*
      SUBROUTINE mpctr(rec,objid,tdt,tutm,alpha,delta,obscod,acct,acca,
     +                 accd,smag,error)
      IMPLICIT NONE
      CHARACTER*80 rec
      CHARACTER*(*) objid
      CHARACTER*6 smag
      DOUBLE PRECISION tdt,tutm,alpha,delta,acct,acca,accd
      INTEGER obscod
      LOGICAL error

      INCLUDE 'trig.h'

      INTEGER year,month,day,ll,pp,ndd,mjd,mjdt,errcod
      DOUBLE PRECISION sec,sect
      CHARACTER*12 mpcname
      CHARACTER*9 chdate
      CHARACTER*3 scale
      LOGICAL err1,err2,err3


      INTEGER lench
      DOUBLE PRECISION tjm1
      EXTERNAL lench,tjm1

      CHARACTER*3 obsstr

      error=.true.
      errcod=1
c Object Name
      mpcname=rec(1:12)
      call iaucod(mpcname,objid,err3)
      if (err3) errcod=17
c Date
      READ(rec,100,ERR=10) year,month,chdate
 100  FORMAT(15X,I4,1X,I2,1X,A9)
      ll=lench(chdate)
      IF(ll.LE.0) GOTO 10
* Position of the decimal point
      pp=INDEX(chdate,'.')
      IF(pp.EQ.0) THEN
          errcod=2
          READ(chdate,*,ERR=10) day
          sec=0
          acct=1
      ELSE
          errcod=3
          READ(chdate(1:pp-1),*,ERR=10) day
          READ(chdate(pp:),*,ERR=10) sec
          sec=sec*86400
          ndd=ll-pp
          acct=10.0d0**(-ndd)
      END IF
      IF(year.LT.1972) THEN
          scale='UT1'
      ELSE
          scale='UTC'
      END IF
      mjd=nint(tjm1(day,month,year,0.d0))
      CALL cnvtim(mjd,sec,scale,mjdt,sect,'TDT')
      tutm=mjd+sec/86400.d0
      tdt=mjdt+sect/86400.d0

      CALL rdanga(rec(33:44),alpha,acca,err1)
      CALL rdanga(rec(45:56),delta,accd,err2)
      errcod=4
      IF(err1.OR.err2.or.err3) GOTO 10

      alpha=alpha*radh
      acca =acca *radh
      delta=delta*radeg
      accd =accd *radeg

      errcod=5
      READ(rec(78:80),*,ERR=10) obsstr
      CALL statcode(obsstr,obscod)
      READ(rec(66:71),166,ERR=10)smag
 166  FORMAT(a6)
      error=.false.

 10   CONTINUE
      IF(error) THEN
          WRITE(*,101) errcod
          tdt=0
          alpha=0
          delta=0
          obscod=-1
          acct=1.d99
          acca=1.d99
          accd=1.d99
      END IF
 101  FORMAT(' mpctr: error code',I3)

      END
* Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: March 12, 1997
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                          M P C R A D                          *
*  *                                                               *
*  *   Transformation of a radar observation (MPC format)          *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    REC1,REC2 -  2 MPC records with R,r in column 15
*
* OUTPUT:   IOBS      -  Obs. type: 2001=r,rdot 2002=r 2003 rdot
*                        if surface bounce add 100
*           OBJID     -  IAU Object designation
*           TDT       -  Time (MJD, TDT)
*           TUTM      -  Time (MJD, UTM)
*           R         -  Range (AU) order 0
*           V         -  Range rate (AU/day) order 0
*           OBSCOD    -  Observatory code
*           ACCT      -  Accuracy of time (d)
*           ACCR      -  Accuracy of range
*           ACCV      -  Accuracy of range rate
*           ERROR     -  Conversion error
*
      SUBROUTINE mpcrad(rec1,rec2,iobs,objid,tdt,tutm,r,v,obscod,
     +          acct,accr,accv,error)
      IMPLICIT NONE
      CHARACTER*80 rec1,rec2
      CHARACTER*(*) objid
      DOUBLE PRECISION tdt,tutm,r,v,acct,accr,accv
      INTEGER obscod,iobs
      LOGICAL error

      INCLUDE 'trig.h'
      INCLUDE 'vlight.h'

      CHARACTER*80 recsw
      INTEGER year,month,day,year2,month2,obscod2
      INTEGER ll,pp,ndd,mjd,mjdt,errcod
      DOUBLE PRECISION sec,sect
      CHARACTER*12 mpcname1,mpcname2
      CHARACTER*9 chdate,chdate2
      CHARACTER*3 scale

      DOUBLE PRECISION hz,dt,df
      LOGICAL rng,vel

      INTEGER lench
      DOUBLE PRECISION tjm1
      EXTERNAL lench,tjm1

      INTEGER isec,iisec
c suface bounce
      LOGICAl surf

      error=.true.
      errcod=1
c sort R and r record
      IF(rec1(15:15).eq.'R')THEN
         IF(rec2(15:15).eq.'r')THEN
c all right, proceed
         ELSE
c not a good pair of records
            errcod=50
            GOTO 10
         ENDIF
      ELSEIF(rec1(15:15).eq.'r')THEN
         IF(rec2(15:15).eq.'R')THEN
c swap
            recsw=rec1
            rec1=rec2
            rec2=recsw
         ELSE
c not a good pair of records
            errcod=50
            GOTO 10
         ENDIF
      ELSE
c not a good pair of records, maybe no radar at all
         errcod=50
         GOTO 10
      ENDIF
c Object Name
      mpcname1=rec1(1:12)
      mpcname2=rec2(1:12)
      if(mpcname1.ne.mpcname2)then
         write(*,*)'name error radar obs. ',mpcname1,' ',mpcname2
         errcod=44
         goto 10
      endif
      call iaucod(mpcname1,objid,error)
      if (error) then
         errcod=17
         goto 10
      endif
c ========== HANDLE DATE AND TIME =================
      READ(rec1,100,ERR=10) year,month,chdate
 100  FORMAT(15X,I4,1X,I2,1X,A9)
      READ(rec2,100,ERR=10) year2,month2,chdate2
      if(year.ne.year2.or.month.ne.month2.or.chdate.ne.chdate2)then
         write(*,*)'date error radar obs. ',year,' ',month,' ',chdate
         write(*,*)'                      ',year2,' ',month2,' ',chdate2
         errcod=45
         goto 10
      endif
      ll=lench(chdate)
      IF(ll.LE.0) GOTO 10
* Position of the decimal point
      pp=INDEX(chdate,'.')
      IF(pp.EQ.0) THEN
          WRITE(*,*)'mpcrad: radar data without correct time'
          STOP
c          errcod=2
c         READ(chdate,*,ERR=10) day
c         sec=0
c         acct=1
      ELSE
          errcod=3
          READ(chdate(1:pp-1),*,ERR=10) day
          READ(chdate(pp:),*,ERR=10) sec
c convert to seconds and round, radar observations are normal points 
c made at integer seconds.
          sec=sec*86400.d0
          isec=nint(sec)
c but, if we get either 59 sec or 1 sec, we rather believe that
c MPC made a mess by rounding time to 10^-5 days 
          iisec=mod(isec,60)
          IF(iisec.eq.1)THEN
            sec=isec-1
          ELSEIF(iisec.eq.59)THEN      
            sec=isec+1
          ELSE
            sec=isec
          ENDIF
c          WRITE(*,*)sec,isec,iisec
          ndd=ll-pp
          acct=10.0d0**(-ndd)
      END IF
      IF(year.LT.1972) THEN
          scale='UT1'
      ELSE
          scale='UTC'
      END IF
      mjd=nint(tjm1(day,month,year,0.d0))
      CALL cnvtim(mjd,sec,scale,mjdt,sect,'TDT')
      tutm=mjd+sec/86400.d0
      tdt=mjdt+sect/86400.d0
c ========== READ TIME DELAY (CONVERT TO ZERO ORDER RANGE) =================
c time delay (in seconds, then convert to days) 
      read(rec1(33:47),FMT='(BZ,F15.10)')dt
      read(rec2(34:47),FMT='(BZ,F14.10)')accr
      dt=dt/86400.d0
      accr=accr/86400.d0
c ========== READ DOPPLER SHIFT =================
c doppler shift (megahertz)
c set obs. type: 2=r,rdot 3=r only 4=rdot only
      read(rec1(48:62),FMT='(BZ,F15.4)')df
      read(rec2(48:62),FMT='(BZ,F15.4)')accv
c ========== READ TRANSMITTER FREQUENCY IN HZ=================
      read(rec1(63:68),FMT='(BZ,F6.1)')hz
      hz=hz*1.d6
      IF(lench(rec2(63:68)).ne.0) THEN
c        continuation of frequency is not handled for now.
         errcod=37
         goto 10
      endif
c ========== CONVERT TO ZERO ORDER RANGE & RANGE-RATE =================
c order zero range
      r=dt*vlight/2.d0      
c order zero range rate
      v=-vlight*df/(hz*2.d0)
c scaling to zeroth order
      accr = accr*vlight/2.d0
      accv = accv*vlight/hz
c deal with surface vs. mass center return:
      if(rec2(33:33).eq.'S')then
c rough fix; it should be deduced from H magnitude, e.g. stored in a file
         surf=.true.
c anyway the weighting needs to be changed in this case???
      elseif(rec2(33:33).eq.'C')then
         surf=.false.         
c observation is already reduced to the asteroid center of mass
c no correction is needed
      else
         errcod=73
         goto 10
      endif
c ========== GET OBSERVATION TYPE =================
      rng=.true.
      vel=.true.
      IF(lench(rec1(33:47)).eq.0) rng=.false.
      IF(lench(rec1(48:62)).eq.0) vel=.false.
c Assign obs type: 2001=r,rdot 2002=r only 2003=rdot only
      IF(rng.and.vel)THEN
         iobs=2001
      ELSEIF(rng)THEN
         iobs=2002
         v=0.d0
         accv=-1.d0
      ELSEIF(vel)THEN
         iobs=2003
         r=0.d0
         accr=-1.d0
      ELSE
         errcod=31
         goto 10
      ENDIF
      IF(surf)iobs=iobs+100
c ============= HANDLE OBSERVATORY CODES ====================
      errcod=5
c obscodes should be identical:
      IF(rec1(69:80).ne.rec2(69:80))THEN
         errcod=46
         GOTO 10
      ENDIF
c read observatory codes
      READ(rec1(69:71),*,ERR=10) obscod
      READ(rec1(78:80),*,ERR=10) obscod2
c divide by 10000 to get the transmitter, remainder is the receiver 
      obscod=obscod*10000+obscod2
      error=.false.

 10   CONTINUE
      IF(error) THEN
         WRITE(*,101) errcod,rec1
 101     FORMAT(' mpcrad: error code',I3,'rec:',A)
         tdt=0
         r=-1.
         v=1.d99
         obscod=-1
         acct=1.d99
         accr=1.d99
         accv=1.d99
      END IF

      END
* Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: October 31, 1997
* version 1.8.0, Steven Chesley, Dec. 14, 1998
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                          J P L I N                            *
*  *                                                               *
*  *         Input of radar  observations from JPL format          *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    FILE      -  File name
*           NOBSX     -  Max dimension of vectors (iobs,t,alpha,delta,
*                        obscod,acct,acca,accd)
*
* OUTPUT:   IOBS      -  Observation type: 2000 for radar
*                                          +100 if surface return
*                                          +2 for time delay or +3 for Doppler
*           TDT       -  Time (MJD, TDT)
*           TUTM      -  Time (MJD, UTM)
*           R         -  Zero order range measurement (AU/day)
*           V         -  Zero order range rate measurement
*           OBSCOD    -  Observatory code: trxsta*10000+recsta
*           ACCT      -  Accuracy of time (d) (always 10^-10)
*           ACCR      -  Accuracy of of range (AU)
*           ACCV      -  Accuracy of range-rate (AU/day)
*           SMAG      -  Dummy string with six blanks
*           RMSMAG    -  0.d0 dummy
*           NOBS      -  Number of observations read
c           OBS       -  successful input flag
*
      subroutine jplin(obs,file,objid,iobs,tdt,tutm,r,v,obscod,
     +                 acct,accr,accv,smag,rmsmag,nobs,nobsx)
      implicit none

      LOGICAL obs

      character*(*) file
      integer nobsx
      character*(*) objid(nobsx)
      integer iobs(nobsx)
      character*6 smag(nobsx)
      double precision tdt(nobsx),tutm(nobsx),r(nobsx),v(nobsx)
      double precision acct(nobsx),accr(nobsx),accv(nobsx),rmsmag(nobsx)
      integer nobs,obscod(nobsx)


      logical error
      character*100 rec
      integer i,l,unit
      external lench
      integer lench

      call filopn(unit,file,'old')

      nobs=0
      do 1 i=1,nobsx
         read(unit,101,end=10) rec
 101     format(a)
         l=lench(rec)
         IF(l.eq.0)GOTO 1
         nobs=nobs+1
         if(nobs.gt.nobsx) stop ' **** mpcin: nobs > nobsx ****'
c alpha is actually r, delta is rdot
         call jplrad(rec,iobs(nobs),
     +        objid(nobs),tdt(nobs),tutm(nobs),r(nobs),
     +        v(nobs),obscod(nobs),acct(nobs),accr(nobs),
     +        accv(nobs),error)
         smag(nobs)='      '
         rmsmag(nobs)=0.d0
         if(error) then
            l=lench(file)
            write(*,102) file(1:l),nobs
 102        format(' **** mpcin: input conversion error ****'/
     .           ' **** file = "',a,'", line',i5,' ****')
            nobs=0
            obs=.false.
            call filclo(unit,' ')
            RETURN
         end if
 1    enddo

c regular ending
 10   continue
      call filclo(unit,' ')
      IF(nobs.eq.0)THEN
         obs=.false.
         write(*,*) 'Warning: Found ',nobs,' obs in ',file
      ELSE
         obs=.true.
         write(*,*) 'Found ',nobs,' obs in ',file
      ENDIF
      return
      end
* Copyright (C) 1999 ORBFIT Consortium
* Version: Sept. 28, 1999
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                          J P L R A D                          *
*  *                                                               *
*  *   Transformation of a radar observation (JPL format)          *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    REC       - Record from JPL (HTML) radar file
*
* OUTPUT:   IOBS      -  Obs. type: 2002=r 2003 rdot
*                        if surface bounce add 100
*           OBJID     -  IAU Object designation
*           TDT       -  Time (MJD, TDT)
*           TUTM      -  Time (MJD, UTM)
*           R         -  Range (AU) order 0
*           V         -  Range rate (AU/day) order 0
*           OBSCOD    -  Observatory code
*           ACCT      -  Accuracy of time (d)
*           ACCR      -  Accuracy of range
*           ACCV      -  Accuracy of range rate
*           ERROR     -  Conversion error
*
      SUBROUTINE jplrad(rec,iobs,objid,tdt,tutm,r,v,obscod,
     +          acct,accr,accv,error)
      IMPLICIT NONE
      CHARACTER*(*) rec
      CHARACTER*(*) objid
      DOUBLE PRECISION tdt,tutm,r,v,acct,accr,accv
      INTEGER obscod,iobs
      LOGICAL error

      INCLUDE 'trig.h'
      INCLUDE 'vlight.h'
      INTEGER errcod
c name
      character*6 number
      character*17 nametmp
      integer lnum,lnam
c time
      integer year,month,day,hour,min,isec,mjd,mjdt
      double precision sec,sect
      character*3 scale
c measurements
      double precision obs,rms,freq
      character*2 unit 
      character*3 surf
      logical range
c obscode
      character*9 trxstr,recstr
      integer iotr,iore
c functions
      INTEGER station,ix
      DOUBLE PRECISION tjm1
      EXTERNAL tjm1

      error=.true.
      errcod=1
c ========== HANDLE OBJECT NAME =================
      READ(rec,101,ERR=10) number,nametmp
 101  FORMAT(a5,1x,a17)
      call rmsp(number,lnum)
      call rmsp(nametmp,lnam)
      if(lnam.gt.9)lnam=9
      if(lnum.eq.0)then
c remove identified object if present ('1991AQ=1994RD')
         ix=index(nametmp,'=')
         if(ix.gt.0)lnam=ix-1
         objid=nametmp(1:lnam)
      else
         objid=number(1:lnum)
      endif
c ========== HANDLE DATE AND TIME =================
      READ(rec,102,ERR=10) year,month,day,hour,min,isec
 102  FORMAT(24X,I4,5(1X,I2))
      IF(year.LT.1972) THEN
          scale='UT1'
      ELSE
          scale='UTC'
      END IF
      sec=(hour*60d0+min)*60d0+isec
      mjd=nint(tjm1(day,month,year,0.d0))
      CALL cnvtim(mjd,sec,scale,mjdt,sect,'TDT')
      tutm=mjd+sec/86400.d0
      tdt=mjdt+sect/86400.d0
      acct=10d-10
c ========== READ MEASUREMENT =================
      READ(rec,103,ERR=10) obs,rms,unit,surf,freq
 103  FORMAT(44x,f13.2,1x,f7.3,1x,a2,1x,a3,1x,f5.0)
      freq=freq*1d6
      if(unit.eq.'us')then
         iobs=2002
         range=.true.
      elseif(unit.eq.'Hz')then
         iobs=2003
         range=.false.
      else
         errcod=37
         goto 10         
      endif
      if(surf.eq.'COM')then
         continue
      elseif(surf.eq.'PP ')then
         iobs=iobs+100
      else
         errcod=47
         goto 10         
      endif
      if(range)then
         r=obs*1d-6/86400d0*vlight/2.d0
         accr=rms*1d-6/86400d0*vlight/2.d0
         v=0d0
         accv=-1d0
      else
         v=-obs*vlight/(freq*2.d0)
         accv=rms*vlight/(freq*2.d0)
         r=0d0
         accr=-1d0
      endif
c This is the old code. Why 
c order zero range
c      r=dt*vlight/2.d0      
c order zero range rate
c      v=-vlight*df/(hz*2.d0)
c scaling to zeroth order
c      accr = accr*vlight/2.d0
c      accv = accv*vlight/hz
c deal with surface vs. mass center return:
c ======================= HANDLE OBSERVATORY CODES ====================
      READ(rec,104,ERR=10) trxstr,recstr
 104  FORMAT(79x,a9,1x,a9)
      errcod=5
      iotr=station(trxstr)
      iore=station(recstr)
      obscod=iotr*10000+iore
c ====================== CLEAN UP =====================================
      error=.false.
 10   CONTINUE
      IF(error) THEN
         WRITE(*,105) errcod,rec
 105     FORMAT(' jplrad: error code',I3,'rec:',A)
         tdt=0
         r=-1.
         v=1.d99
         obscod=-1
         acct=1.d99
         accr=1.d99
         accv=1.d99
      ENDIF
      END
c ====================== station function =============================      
c compute observatory code from string
      integer function station(stastr)
      implicit none
      character*(*) stastr
c We have the following stations:
c     'Arecibo  ' = 251
c     'DSS 13   ' = 252
c     'DSS 14   ' = 253
c     'Haystack ' = 254
c     'Evpatoria' = 255 (not in MPC file...)
c     'Greenbank' = 256
      if(stastr.eq.'Arecibo  ')then
         station=251
      elseif(stastr.eq.'DSS 13   ')then
         station=252
      elseif(stastr.eq.'DSS 14   ')then
         station=253
      elseif(stastr.eq.'Haystack ')then
         station=254
      elseif(stastr.eq.'Evpatoria')then
         station=255
      elseif(stastr.eq.'Greenbank')then
         station=256
      else
         station=500
         write(*,*)'jplrad: unknown station: ',stastr
      endif
      return
      end
* Copyright (C) 1996 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: August 27, 1996
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         S E S S A G                           *
*  *                                                               *
*  *       Transform an angle into sessagesimal notation           *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    ANG       -  Angle
*
* OUTPUT:   SIGA      -  Sign
*           INTA      -  Integer part
*           MINA      -  Minutes
*           SECA      -  Seconds
*
      SUBROUTINE sessag(ang,siga,inta,mina,seca)
      IMPLICIT NONE

      DOUBLE PRECISION ang,seca
      INTEGER inta,mina,truncat
      CHARACTER*1 siga

      DOUBLE PRECISION anga,u


      IF(ang.GE.0.d0) THEN
          anga=ang
          siga='+'
      ELSE
          anga=ABS(ang)
          siga='-'
      END IF

      inta=truncat(anga,1d-10)
      u=(anga-inta)*60.d0
      u=ABS(u)
      mina=truncat(u,1d-8)
      seca=(u-mina)*60.d0

      END

c Copyright Orbfit Consortium 1999
c this routine determines an asteroid radius if needed for radar
      subroutine astrad(objid,iobs,m)
      implicit none
      integer m,iobs(m)
      character*9 objid(m)
c observation numbers: maximum, space left
c      INCLUDE 'parobx.h'
c      INTEGER nlef
c ====================================

      integer unit,ln,lnam,lnum,isav,i
      character*18 oid,number,name
      double precision hmag,diam,exponent
      logical needed
c radius of asteroid common
      include 'radius.h'

c First find out if we need the radius at all
      needed=.false.
      do i=1,m
         if(iobs(i)/100.eq.21)then
            needed=.true.
            isav=i
            goto 100
         endif
      enddo
 100  continue
      if(.not.needed)then
c bail out
c         write (*,*) 'Asteroid radius not required.'
         radius=-1d0
         return
      else
c get radius
         call filopl(unit,'astorb.rad')
         oid=objid(isav)
         call rmsp(oid,ln)
c read records in an infinite loop (it's F(UGLY)77)         
 1       continue
            read(unit,101,end=109) number,name,hmag,diam
            call rmsp(number,lnum)
            call rmsp(name,lnam)
            if(  oid(1:ln).eq.number(1:lnum) .or.
     +           oid(1:ln).eq.name(1:lnam)    )THEN
c              found object in file
               if(diam.gt.0d0)then
                  radius=diam/2d0/1.4998d8
                  write(*,102)'RADAR: Using radius = ',
     +                 radius*1.4998d8,' km from IRAS diameter.'
               else
                  exponent=0.5d0*(6.3d0-log10(0.2d0)-0.4d0*hmag)
                  radius=10**exponent/2d0/1.4998d8
                  write(*,102)'RADAR: Using radius = ',
     +                 radius*1.4998d8,' km from absolute magnitude.'
               endif
c exit infinite loop
               goto 109
            endif
         goto 1
      endif
      write(*,*)'*** astrad warning: ',name,' not found in astorb.dat.'
      radius=1d0/1.4998d8
      write(*,*)'RADAR: Using radius = ',radius*1.4998d8,' km.'
 109  call filclo(unit,' ')
      return
c ==============================================================
c                  num   nam    comp   Hmag  Gmag    col   diam
 101        format(a5,1X,A18,1X,15x,1X,f5.2,1X,5x,1X,4x,1X,f5.1)
 102        format(a,f6.2,a)
c ==============================================================
      end
c
* Copyright (C) 1997 by Mario Carpino
* Version: February 24, 1997
*
*  *****************************************************************
*  *                                                               *
*  *                         R D A N G A                           *
*  *                                                               *
*  *    Read an angle with its accuracy from a character string    *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    STRING    -  Character string
*
* OUTPUT:   ANGLE     -  Angle (no unit conversion is performed)
*           ACC       -  Angle accuracy
*           ERROR     -  Conversion error
*
      SUBROUTINE rdanga(string,angle,acc,error)
      IMPLICIT NONE
      CHARACTER*(*) string
      DOUBLE PRECISION angle,acc
      LOGICAL error

* Max string length
      INTEGER lx
      PARAMETER (lx=200)

      CHARACTER*(lx) c1,c,field
      INTEGER l,isig,nf,i,pp,iv,ll
      LOGICAL nospli
      DOUBLE PRECISION fact,rv


      INTEGER lench,nitchs
      EXTERNAL lench,nitchs

      error=.true.
      IF(lench(string).GT.lx) STOP '**** rdanga: LEN(string) > lx ****'
      c1=string
      CALL norstr(c1,l)

* The sign may be separated from the value by blank space
      isig=1
      IF(c1(1:1).EQ.'+') THEN
          c=c1(2:)
          CALL norstr(c,l)
      ELSEIF(c1(1:1).EQ.'-') THEN
          isig=-1
          c=c1(2:)
          CALL norstr(c,l)
      ELSE
          c=c1
      END IF

      nf=nitchs(c)
      IF(nf.LT.1.OR.nf.GT.3) RETURN
      angle=0
      fact=1
      DO 1 i=1,nf
          CALL stspli(c,' ',field,nospli)
          IF(nospli) RETURN
          pp=INDEX(field,'.')
          IF(pp.GT.0) THEN
              IF(i.NE.nf) RETURN
              READ(field,*,err=10,end=10) rv
          ELSE
              READ(field,*,err=10,end=10) iv
              rv=iv
          END IF
          angle=angle+fact*rv
          IF(i.EQ.nf) THEN
              CALL norstr(field,ll)
              pp=INDEX(field,'.')
              IF(pp.EQ.0) THEN
                  acc=1
              ELSE
                  ll=lench(field)
                  acc=10.0d0**(pp-ll)
              END IF
              acc=acc*fact
          ELSE
              fact=fact/60.d0
          END IF
 1    CONTINUE
      angle=angle*isig
      error=.false.

 10   CONTINUE
      END

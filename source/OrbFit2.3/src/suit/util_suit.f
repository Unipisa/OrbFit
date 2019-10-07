* Copyright (C) 2000 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: November 30, 2000
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         A N G V C F                           *
*  *                                                               *
*  *           Conversion factor for angular velocities            *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    UNIT      -  Description of unit
*
* OUTPUT:   CF        -  Conversion factor from rad/d to "unit"
*           ERROR     -  Error flag
*
      SUBROUTINE angvcf(unit,cf,error)
      IMPLICIT NONE

      CHARACTER*(*) unit
      DOUBLE PRECISION cf
      LOGICAL error

      INCLUDE 'trig.h'

      DOUBLE PRECISION ct
      CHARACTER*20 unang,untim
      LOGICAL nospli

      error=.true.

* Input string must be of the form 'unang/untim',
* where: "unang" is the angular unit (rad/deg/arcmin/arcsec)
*        "untim" is the time unit (d/h/min/s);
* for instance: 'deg/d' or 'arcsec/min'
      untim=unit
      CALL stspli(untim,'/',unang,nospli)
      IF(nospli) RETURN
      CALL locase(unang)
      CALL locase(untim)

* List of supported angular units
      IF(unang.EQ.'rad') THEN
          cf=1
      ELSEIF(unang.EQ.'deg') THEN
          cf=degrad
      ELSEIF(unang.EQ.'arcmin') THEN
          cf=degrad*60
      ELSEIF(unang.EQ.'''') THEN
          cf=degrad*60
      ELSEIF(unang.EQ.'arcsec') THEN
          cf=degrad*3600
      ELSEIF(unang.EQ.'"') THEN
          cf=degrad*3600
      ELSE
          RETURN
      END IF

* List of supported time units
      IF(untim.EQ.'d') THEN
          ct=1
      ELSEIF(untim.EQ.'day') THEN
          ct=1
      ELSEIF(untim.EQ.'h') THEN
          ct=24
      ELSEIF(untim.EQ.'hour') THEN
          ct=24
      ELSEIF(untim.EQ.'min') THEN
          ct=24*60
      ELSEIF(untim.EQ.'s') THEN
          ct=24*3600
      ELSE
          RETURN
      END IF
      cf=cf/ct
      error=.false.

      END
c =====================================================================
c TEE (as in unix shell)
c =====================================================================
      subroutine tee(iun,string)
      implicit none
      integer iun
      character*(*) string
      integer le
      le=index(string,'=')
      write(*,*)string(1:le-1)
      write(iun,*)string(1:le-1)
      return
      end
c ================================================================
c MENU
c ================================================================
      SUBROUTINE menu(ifl,menunam,nopt,s,
     +             s1,s2,s3,s4,s5,s6,s7,s8,s9,s10)
      IMPLICIT NONE
      INTEGER ifl,nopt,ll,iunit,lench
      CHARACTER*(*) s,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10
      CHARACTER*20 menunam
      CHARACTER*120 helpfi,ddocd1
      INCLUDE 'doclib.h'
c
      IF(nopt.lt.2.or.nopt.gt.10)then
         WRITE(*,*) ' this menu can only handle betw. 2 and 10 options'
         ifl=0
         return
      ENDIF
 3    continue
      ll=index(s,'=')
      WRITE(*,*) s(1:ll-1)
      ll=index(s1,'=')
      WRITE(*,*)' 1 = ', s1(1:ll-1)
      ll=index(s2,'=')
      WRITE(*,*)' 2 = ', s2(1:ll-1)
      IF(nopt.lt.3) goto 2
      ll=index(s3,'=')
      WRITE(*,*)' 3 = ', s3(1:ll-1)
      IF(nopt.lt.4) goto 2
      ll=index(s4,'=')
      WRITE(*,*)' 4 = ', s4(1:ll-1)
      IF(nopt.lt.5) goto 2
      ll=index(s5,'=')
      WRITE(*,*)' 5 = ', s5(1:ll-1)
      IF(nopt.lt.6) goto 2
      ll=index(s6,'=')
      WRITE(*,*)' 6 = ', s6(1:ll-1)
      IF(nopt.lt.7) goto 2
      ll=index(s7,'=')
      WRITE(*,*)' 7 = ', s7(1:ll-1)
      IF(nopt.lt.8) goto 2
      ll=index(s8,'=')
      WRITE(*,*)' 8 = ', s8(1:ll-1)
      IF(nopt.lt.9) goto 2
      ll=index(s9,'=')
      WRITE(*,*)' 9 = ', s9(1:ll-1)
      IF(nopt.lt.10) goto 2
      ll=index(s10,'=')
      WRITE(*,*)'10 = ', s10(1:ll-1)
c
c room to increase
c
 2    WRITE(*,103)
 103  format(' 0 = exit; -1=help')
      WRITE(*,*)' selection?  '
      read(*,*,err=3)ifl
c wrong flag and exit
 4    IF(ifl.lt.-1.or.ifl.gt.nopt)THEN
            WRITE(*,*)ifl,' option not understood'
            goto 3
      ELSEIF(ifl.eq.-1)THEN
         ddocd1=ddocd
         ll=lench(ddocd1)
         helpfi=ddocd1(1:ll)//'/'//menunam
         CALL rmsp(helpfi,ll)
         helpfi=helpfi(1:ll)//'.help'
         CALL filopn(iunit,helpfi,'OLD')
         CALL filcat(iunit)
         CALL filclo(iunit,' ')
         WRITE(*,*)' selection?  '
         read(*,*,err=3)ifl
         GOTO 4
      ENDIF
      RETURN
      END
c ===========================================
      SUBROUTINE filcat(iunit)
      IMPLICIT NONE
      INTEGER iunit,i,imax,ll,lcom,lench
      PARAMETER (imax=100)
      CHARACTER*100 record
c =========================================
      DO i=1,imax
        READ(iunit,100,END=2)record
 100    FORMAT(a)
        ll=lench(record)
c        WRITE(*,*)ll
c comments begin with %, for TeX compatibility
        lcom=index(record,'%')
        IF(lcom.gt.1)THEN
           WRITE(*,100)record(1:lcom-1)
        ELSEIF(lcom.eq.0)THEN
           IF(ll.gt.0)THEN
              WRITE(*,100)record(1:ll)
           ELSE
              WRITE(*,*)
           ENDIF
        ENDIF
      ENDDO
 2    RETURN
      END
c ===============================================
c  NIGHTS
c function computing number of nights of observations
c Copyright A. Milani, OrbFit consortium, 21/9/2000
c ===============================================
      INTEGER FUNCTION nights(m,iobs,aln,den,tut,idsta,sel,
     +     rmsa,rmsd) 
      IMPLICIT NONE
c =========OBSERVATIONS =========================
      INCLUDE 'parobx.h'
c number of observations
      INTEGER m
c observations: alpha, delta, time (UT), station code 
      double precision aln(m),den(m),tut(m)
      integer idsta(m),iobs(m)
c RMS of observation error
      double precision rmsa(m),rmsd(m)
c selection flags; number of observations not discarded
      INTEGER sel(m)
c ===========END INTERFACE=======================
      double precision t1
      integer i
c ===============================================
      t1=tut(1)
      nights=1
      DO i=2,m
c This is the most trivial algorithm: if there is a 16 hours interval
c between two observations, they belong to different nights. 
c But a night containing only observations being discarded does not count
        IF(tut(i)-t1.gt.0.66d0.and.sel(i).gt.0)THEN
           nights=nights+1
           t1=tut(i)
        ENDIF
      ENDDO
      RETURN
      END
c ==================================================================
c this should be improved to cover the following "strange" cases:
c 1) an asteroid which is observed around the clock; the above algorithm
c    would rate it as a 3-nighter after 48 hours!
c 2) an isolate observation, maybe even with degraded accuracy, 
c    might not qualify as 'one night'
c 3) a single observation with radar (maybe, in the future, from a space misison)
c    could be considered equivalent to many, many nights
c All the data are available in the function call to add these improvements, 
c e.g. the RMS (rmsa, rmsd), the obs. type (iobs, 200x for radar)
c and the outlier rejection flag (sel)
c Note: the computation is currently done using UT; a leap second does not matter






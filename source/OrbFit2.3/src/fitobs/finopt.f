c ===================================================================
c FINOPT
c ===================================================================
c input options, for the propagator and the specific main program
c input: progna = program name (6 characters)
c        run    = run identifier (80 characters, up to 76 non-blank) 
      SUBROUTINE finopt(progna,run,astna0,astnap,iun20,iun8)
      implicit none
      character*6 progna
      character*80 run
      integer iun20,iun8
      CHARACTER*(*) astna0,astnap
c ==========END INTERFACE============================================
      integer le,iunit
      character*90 file
      logical found
      LOGICAL ireq,fail1,fail
      CHARACTER*7 prognp
      INCLUDE 'proout.h'
c =============================
c read option for propagator
      CALL libini
      CALL namini
c default options 
      CALL filopl(iunit,'fitobs.def')
      CALL rdnam(iunit)
      CALL filclo(iunit,' ')
c =============================
c particular options for this run
      CALL rmsp(run,le)
      file=run(1:le)//'.fop'
      INQUIRE(FILE=file,EXIST=found)
      IF(found) THEN
        CALL filopn(iunit,file,'OLD')
        CALL rdnam(iunit)
        CALL filclo(iunit,' ')
      ELSE
        write(*,*)'Option file not found: ',file
        write(*,*)'*** Using default options. ***'
c        stop 
      ENDIF
c =============================
c check for non-existing options
      call rmsp(progna,le)
      file=progna(1:le)//'.key'
      CALL rdklst(file)
      CALL chkkey
c ==============================
c initialisations for Gauss method
      CALL iodini
c      IF(fail) STOP 'finopt: error on gaussian initialization'
c =====================================================================
c Output files: for control and results, for covariance
      CALL rmsp(run,le)
      file=run(1:le)//'.fou'
      call filopn(iun20,file,'UNKNOWN')
      file=run(1:le)//'.fga'
      call filopn(iun8,file,'UNKNOWN')
c =====================================================================
c Output files: for errors, close approaches, for propagator parameters
      CALL rmsp(run,le)
      file=run(1:le)//'.err'
      call filopn(ierrou,file,'UNKNOWN')
      numerr=0
      file=run(1:le)//'.clo'
      call filopn(iuncla,file,'UNKNOWN')
      numcla=0
      file=run(1:le)//'.pro'
      call filopn(ipirip,file,'UNKNOWN')
c =============================
c read option for physical model and integration method
      CALL rmodel(run)
c asteroid name(s)
c =============SOME ASTEROID NAME NEEDED===================
      fail=.false.
**      ireq=.true.
      ireq=.false.
      prognp=progna//'.'
      CALL rmsp(prognp,le)
      CALL rdncha(prognp,'astna0',astna0,ireq,
     +            found,fail1,fail)
      IF(.not.found.or.fail)THEN
**         WRITE(*,*)' finopt: first arc ast. name required'
**         STOP
         CALL rmsp(run,le)
         WRITE(*,*)' finopt: first arc ast. name not found, ',
     +        'using ',run(1:le)
         astna0=run(1:le)
      ENDIF
      ireq=.false.
      CALL rdncha(prognp,'astnap',astnap,ireq,
     +            found,fail1,fail)
      IF(.not.found.or.fail)THEN
         WRITE(*,*)' asteroid ',astna0
         astnap=' '
      ELSE
         WRITE(*,*)' two asteroids  ', astna0,astnap
      ENDIF
      RETURN
      END


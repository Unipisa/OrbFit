c ===================================================================
c OPTPRO
c ===================================================================
c input options, for the propagator and the specific main program
c input: progna = program name (6 characters)
c        run    = run identifier (80 characters, up to 76 non-blank) 
      SUBROUTINE optpro(progna,eledir)
      implicit none
      character*6 progna
      character*(*) eledir
c ==========END INTERFACE============================================
      CHARACTER*80 run
      integer le,iunit
      character*12 file
      logical found
      LOGICAL ireq,fail1,fail
      CHARACTER*7 prognp
c =============================
c read option for propagator
      CALL libini
      CALL namini
c default options are only for propag
      CALL filopl(iunit,'propag.def')
      CALL rdnam(iunit)
      CALL filclo(iunit,' ')
c =============================
c particular options for this run
      call rmsp(progna,le)
      file=progna(1:le)//'.nop'
      INQUIRE(FILE=file,EXIST=found)
      IF(found) THEN
        CALL filopn(iunit,file,'OLD')
        CALL rdnam(iunit)
        CALL filclo(iunit,' ')
      ELSE
        write(*,*)'**** file not found: ',file
        write(*,*)'******* ',progna,' abnormal end ****'
        stop 
      ENDIF
c =============================
c check for non-existing options
      call rmsp(progna,le)
      file=progna(1:le)//'.key'
      CALL rdklst(file)
      CALL chkkey
c ==============================
c read option for physical model and integration method
      CALL rmodel(progna)
c ============= WHERE TO FIND/PUT ASTEROID ELEMENTS===================
      fail=.false.
      ireq=.true.
      prognp=progna//'.'
      CALL rdncha(prognp,'eledir',eledir,ireq,
     +            found,fail1,fail)
      IF(fail)THEN
         WRITE(*,*)' optpro: asteroid elements directory required'
         STOP
      ENDIF
c ============= CONTROL OF SEVERAL OPTIONAL SERVICES ===================

c ================================================
      RETURN
      END


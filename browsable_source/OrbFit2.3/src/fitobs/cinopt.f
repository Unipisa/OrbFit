c ===================================================================
c CINOPT
c ===================================================================
c input options, for the catalog propagator main program
c input: progna = program name (6 characters)
c        run    = run identifier (80 characters, up to 76 non-blank) 
      SUBROUTINE cinopt(progna,run,iun20,tref,catnam0,catnam1,covpro)
      implicit none
      character*6 progna
      character*80 run
      integer iun20
      double precision tref
      CHARACTER*80 catnam0,catnam1
      logical covpro
c ==========END INTERFACE============================================
      integer le,iunit,l0,l1
      character*20 file
      LOGICAL ireq,fail1,fail,found,ex
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
      CALL rmsp(run,le)
      file=run(1:le)//'.mop'
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
c initialisations for Gauss method
c      CALL iodini
c =====================================================================
c Output files: for control and results, for covariance
      CALL rmsp(run,le)
      file=run(1:le)//'.mou'
      call filopn(iun20,file,'UNKNOWN')
c =============================
c read option for physical model and integration method
      CALL rmodel(run)
c ================reference time============================
      ireq=.true.
      fail=.false.
      prognp=progna//'.'
      CALL rmsp(prognp,le)
      CALL rdnrea(prognp,'tref',tref,ireq,
     +            found,fail1,fail)
      WRITE(*,*)  'reference time ', tref
      WRITE(iun20,*) 'reference time ', tref
c ====================================================
c read the name of the elements file, inquire
      ireq=.true. 
      CALL rdncha(prognp,'catnam0',catnam0,ireq,found,
     +     fail1,fail)
      CALL rmsp(catnam0,l0)
      INQUIRE(file=catnam0(1:l0),exist=ex)
      IF(.not.found.or.fail.or..not.ex)THEN
         write(*,*) found,fail,ex
         write(*,*) 'catalog not found:',catnam0(1:l0)
         STOP
      ENDIF
      CALL rdncha(prognp,'catnam1',catnam1,ireq,found,
     +     fail1,fail)
      CALL rmsp(catnam1,l1)
      IF(.not.found.or.fail)THEN
         write(*,*) 'catalog name not found:',catnam1(1:l1)
         STOP
      ENDIF
c ================generate catalog with covariances=============
      ireq=.true.
      CALL rdnlog(prognp,'covpro',covpro,ireq,
     +            found,fail1,fail)
      IF(fail)THEN
         WRITE(*,*)' cinopt: generate covariance logical required'
         STOP
      ENDIF
      WRITE(*,*)  ' generate catalog with covariances is ',covpro
      WRITE(iun20,*)' generate catalog with covariances is ',covpro
c ====================================================
      RETURN
      END

c ================================================================== 
c FINOBS
c ================================================================== 
c Arc  observations input
      subroutine finobs(progna,iar,astna0,objid,obs0,m,iobs,tau,aln,den,
     +    tut,idsta,sel,rmsa,rmsd,rmsmag,smag,nlef,rwofi0,iun20)
      implicit none
c =======INPUT==================================
c ========= file names, i/o control ============
c program name
      character*(*) progna
c asteroid name (9 characters)
      character*(*) astna0
c messages unit
      integer iun20
c flag for arcs 1-2
      integer iar
c ===== observational data ===========================
c observation numbers: maximum, space left
      include 'parobx.h'
      integer nlef
c =======OUTPUT==================================
c residuals and weights file name
      character*60 rwofi0
      character*(*) objid(nlef)
c successful input flag
      logical obs0
c ===== observational data ===========================
c observation number
      integer m
c observations: alpha, delta, time (ET and UT), station code, type
      double precision aln(nlef),den(nlef),tau(nlef),tut(nlef)
      integer idsta(nlef),iobs(nlef)
c app. magnitude
      character*6 smag(nlef) 
c selection flag 0=discard 1=select 2=prelim
      integer sel(nlef)
c RMS of observation error, of magnitude
      double precision rmsa(nlef),rmsd(nlef),rmsmag(nlef)
c change flag
      logical change
c ============END INTERFACE==========================
      CHARACTER*18 nam0
      INTEGER lnam
c file names, , for observations
      character*60 obsdir0
c input flags 
      logical fail, fail1,found,ireq
c string length (after blank removal)
      integer le
c loop indexes
      integer ll,lench
c program name with dot
      character*7 prognp
c directory char
      INCLUDE 'sysdep.h'
c precedence to .obs file with respect to .rwo file: false for use in fitobs
      LOGICAL precob
      precob=.false.
c      precob=.true.
c ====================================================
c exclude non existing arcs
      nam0=astna0
      CALL rmsp(nam0,lnam)
      IF(lnam.eq.0)THEN
c this arc does not exist
         obs0=.false.
         m=0
         RETURN
      ENDIF
c find directory of observations
      prognp=progna//'.'
      ireq=.false.
      fail=.false.
      if(iar.eq.1)then
         CALL rdncha(prognp,'obsdir0',obsdir0,ireq,found,
     +            fail1,fail)
      elseif(iar.eq.2)then
         CALL rdncha(prognp,'obsdirp',obsdir0,ireq,found,
     +            fail1,fail)
      endif
      IF(.not.found.or.fail) then
         write(*,*)' obs. directory option not found, using ''obsdata'''
         obsdir0='obsdata'
      ENDIF
      CALL rmsp(obsdir0,le)
c find observations file name
      fail=.false.
      IF(iar.eq.1)THEN
         CALL rdncha(prognp,'nam0',nam0,ireq,found,
     +            fail1,fail)
      ELSEIF(iar.eq.2 )THEN
            CALL rdncha(prognp,'namp',nam0,ireq,found,
     +            fail1,fail)
      ENDIF
      IF(.not.found.or.fail)THEN
         nam0=astna0
         CALL rmsp(nam0,lnam)
         WRITE(*,*)' obs. filename option not found, using name ',
     +        nam0(1:lnam),';'
      ENDIF
c input data 
      CALL inobs(obsdir0,nam0,precob,objid,obs0,m,iobs,tau,
     +    aln,den,tut,idsta,sel,rmsa,rmsd,rmsmag,smag,nlef,iun20,change)
c compute name of .rwo file
      ll=lench(obsdir0)
      rwofi0=obsdir0(1:ll)//dircha
      ll=lench(rwofi0)
      rwofi0=rwofi0(1:ll)//nam0//'.rwo'
      CALL rmsp(rwofi0,ll)
      RETURN
      END






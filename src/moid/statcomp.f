c     ********************************************************************
c     *******************  PROGRAM   S T A T C O M P  ********************
c     ********************************************************************
c     ************* written by GIOVANNI F. GRONCHI (2001) ****************
c     ********************************************************************
c     *************** E-MAIL gronchi@mail.dm.unipi.it ********************
c     ********** Department of Mathematics, UNIVERSITY of PISA ***********
c     ********************************************************************
c     ====================================================================
      PROGRAM statcomp2
      IMPLICIT NONE
c     ====================================================================
c     asteroid counter
      INTEGER astnum,nast,nastmax,num_obj,loop
      PARAMETER (nastmax=2000)
c     asteroid name
      CHARACTER*9  name
c     orbital elements
      DOUBLE PRECISION elkep(6),elem0(6)
c     for subroutine rdelem
      DOUBLE PRECISION cov0(6,6),norm0(6,6), mass
      DOUBLE PRECISION hmag,gmag,enne
c     planets data
      INCLUDE 'pldata.h'
      CHARACTER*3 eltype
      CHARACTER*80 comele
c     successful input flags
      LOGICAL defcov,defelem,neodys_elem
c     current time (JD)
      DOUBLE PRECISION t1
c     elements
      DOUBLE PRECISION elem1(6)
c     available data
      LOGICAL ok
c     file names 
      CHARACTER*60 file,elefil,eledir,elemfiles(4)
      INTEGER le,lnam,num_files
      CHARACTER*6 progna
c     file identifiers
      INTEGER iunlog,iunit,iunhand,iun1,iun2,iun3,iun4
c     loop index
      INTEGER i
c     Earth orbital elements
      DOUBLE PRECISION eqp(6),ekpl(6)
      DOUBLE PRECISION apl,epl,Ipl,omegapl,Ompl,lpl
      INCLUDE 'trig.h'
      INCLUDE 'proout.h'
      INCLUDE 'neoopt.h'
      CHARACTER*60 filcat
c     input options
      progna='statco'      
      CALL optpro(progna,eledir)
c     =======================================================

c     check files
      CALL filopn(iun4,'CHECK.hi-st','unknown')
      CALL filopn(iun3,'CHECK.morse_weier','unknown')
      CALL filopn(iun2,'CHECK.warning','unknown')
      CALL filopn(iun1,'CHECK.solvsys','unknown')
      CALL filopn(iunlog,'statcomp.log','unknown')
      CALL filopn(iunhand,'statcomp.names.fail','unknown')

      astnum = 0
c     =========== ENTER MAIN LOOP ===========================
      DO 1 loop=1,nastmax
         astnum = astnum + 1
c     read next name
         write(*,*) 'Name (Ctrl-D to quit)?'
 5       READ(*,100,end=999) name
c     5       READ(iunlst,100,end=999) name
 100     FORMAT(a9)
         CALL rmsp(name,lnam)
         IF(name(1:1).eq.'!')THEN
            WRITE(iunhand,*) name,'\t\tCommented near line ',nast
            GOTO 5
         ENDIF
         defelem=.false.
         defcov=.false.
c     error file to be opened
         CALL filnam('./err',name,'err',file,le)
         CALL filopn(ierrou,file(1:le),'unknown')
         numerr=0
c     put flag in each log file
         WRITE(iunlog,*)'***************************************'
         WRITE(iunlog,*)' Processing ',name,' number ',astnum
         WRITE(iunlog,*)'***************************************'
         WRITE(*,*)'\n ***************************************'
         WRITE(*,*)'  Processing ',name,' number ',astnum
         WRITE(*,*)'***************************************\n'
c     input file
         CALL filnam(eledir,name,'eq1',elefil,le)         
         INQUIRE(file=elefil(1:le),exist=neodys_elem)
         IF(neodys_elem)THEN
            elemfiles(1) = elefil(1:le)
            num_files = 1
         ELSE
            WRITE(ierrou,*)' elements not found for ',name
            numerr=numerr+1
            GOTO 2
         ENDIF
         num_obj=1
         CALL rdelem(iunlog,name,num_obj,elemfiles,num_files,defelem,
     +        defcov,eltype,t1,elem1,cov0,norm0,mass,hmag,gmag,comele)
         IF(.not.defelem)THEN
            WRITE(iunlog,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            WRITE(iunlog,*)'!!WARNING! WARNING! WARNING! WARNING!!'
            WRITE(iunlog,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            WRITE(iunlog,*)'Asteroid ',name,' number ',astnum, 
     +           ' not found in asteroid catalog.'
            WRITE(iunhand,*)name,'\t\tnot in catalogs.'
            GOTO 2
         ENDIF
         
         CALL coocha(elem1,eltype,bigg,elkep,'KEP',enne)
c     output file 
         CALL filnam('.',name,'statpts',elefil,le)         
         CALL filopn(iunit,elefil(1:le),'unknown')
c     elements of the Earth
         CALL earth(t1,eqp)
         CALL coocha(eqp,'EQU',bigg,ekpl,'KEP',enne)

c     ******************************************************
         CALL writestat(name,lnam,iunit,iun1,iun2,iun3,iun4,
     *    ekpl,elkep)    
c     ******************************************************

c     close output file
         CALL filclo(iunit,' ')
c     clean up err and clo files
 2       CONTINUE
c     close error file
         IF(numerr.gt.0)THEN
            CALL filclo(ierrou,' ')
         ELSE
            CALL filclo(ierrou,'DELETE')
         ENDIF
c     end of main loop
 1    ENDDO
      STOP ' Warning, too many objects. Reached end of loop.'
c     EOF exit point
 999  CONTINUE
      CALL filclo(iunlog,' ')
      CALL filclo(iunhand,' ')
      CALL filclo(iun4,' ')
      CALL filclo(iun3,' ')
      CALL filclo(iun2,' ')
      CALL filclo(iun1,' ')
      
      STOP
      END

c ===========MODULE FILE_OPER=====================
c FILE OPENING/CLOSE AND UNIT ASSIGNMENT:
c CONTAINS
c SUBROUTINES
c
c filopn	file opening with unit assignment
c filclo	file closing with unit release
c filass	unit assignment without file opening
c filopl	unit assignment and file opening (with library search)
c filopf	unit assignment and file opening (with library search and
c		no abort if the file does not exist)
c dlifex	delete a file if esists
c libini	inizialization of default library directory
c filnam        composition of file name from dir, name suffix
c splinam       split asteroid name (identification, multiple solution)
c fidinam       find path for file
c dircom2       used by fidinam
c
c
* Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: December 5, 1997
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         F I L O P N                           *
*  *                                                               *
*  *               Unit allocation and file opening                *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    NAME      -  File name to be opened
*           STATUS    -  Open status
*
* OUTPUT:   IUN       -  Allocated unit
*
      SUBROUTINE filopn(iun,name,status)
      IMPLICIT NONE

      INCLUDE 'comfil.h'

      INTEGER iun,ll,ls,i
      CHARACTER*(*) name,status
      LOGICAL opnd

      INTEGER lench
      EXTERNAL lench

      INQUIRE(FILE=name,OPENED=opnd)
      IF(opnd) THEN
          ll=lench(name)
          WRITE(*,102) name(1:ll)
          STOP '**** filopn: abnormal end ****'
      END IF
 102  FORMAT(' **** filopn: internal error (01) ****'/
     +       ' **** FILE: ',A,' ****')

      IF(iicfil.NE.36) THEN
          DO 1 i=iunf1,iunf2
 1        allunt(i)=.false.
          iicfil=36
      END IF

      DO 2 iun=iunf1,iunf2
      IF(allunt(iun)) GOTO 2
      OPEN(iun,FILE=name,STATUS=status,ERR=3)
      filnam(iun)=name
      allunt(iun)=.true.
      RETURN
 2    CONTINUE

      STOP '**** filopn: all units are already allocated ****'

 3    CONTINUE
      ll=lench(name)
      ls=lench(status)
      WRITE(*,101) name(1:ll),status(1:ls)
 101  FORMAT(' **** filopn: cannot OPEN file "',A,'" (status=',A,
     +       ') ****')
      STOP '**** filopn: abnormal end ****'
      END
* Copyright (C) 1996-1998 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: January 12, 1998
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         F I L C L O                           *
*  *                                                               *
*  *                File closing and unit release                  *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    IUN       -  Unit to be closed
*           STATUS    -  Close status (if blank, it is not used)
*
      SUBROUTINE filclo(iun,status)
      IMPLICIT NONE

      INCLUDE 'comfil.h'

      INTEGER iun,ls
      CHARACTER*(*) status

      INTEGER lench
      EXTERNAL lench

      IF(iicfil.NE.36) STOP '**** filclo: internal error (01) ****'
      IF(iun.LT.iunf1.OR.iun.GT.iunf2)
     +          STOP '**** filclo: internal error (02) ****'
      IF(.NOT.allunt(iun)) THEN
          WRITE(*,200) iun
          STOP
      END IF
 200  FORMAT(' **** filclo: unit',i4,' is not opened ****')

      ls=lench(status)
      IF(ls.LE.0) THEN
          CLOSE(iun)
      ELSE
          CLOSE(iun,STATUS=status)
      END IF

      allunt(iun)=.false.
      filnam(iun)=' '

      END
* Copyright (C) 1996 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: September 19, 1996
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         F I L A S S                           *
*  *                                                               *
*  *               Unit allocation for file opening                *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    NAME      -  File name to be opened
*
* OUTPUT:   IUN       -  Allocated unit
*
      SUBROUTINE filass(iun,name)
      IMPLICIT NONE

      INCLUDE 'comfil.h'

      INTEGER iun,i
      CHARACTER*(*) name

      IF(iicfil.NE.36) THEN
          DO 1 i=iunf1,iunf2
 1        allunt(i)=.false.
          iicfil=36
      END IF

      DO 2 iun=iunf1,iunf2
      IF(allunt(iun)) GOTO 2
      filnam(iun)=name
      allunt(iun)=.true.
      RETURN
 2    CONTINUE

      STOP '**** filass: all units are already allocated ****'

      END
* Copyright (C) 1996 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: September 19, 1996
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         F I L O P L                           *
*  *                                                               *
*  *        Unit allocation and file opening (STATUS='old')        *
*  *         (searching also in default library directory)         *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    NAME      -  File name to be opened
*
* OUTPUT:   IUN       -  Allocated unit
*
      SUBROUTINE filopl(iun,name)
      IMPLICIT NONE

      INCLUDE 'comfil.h'

* NEEDED common blocks:
      INCLUDE 'comlib.h'

      INTEGER iun,ll,i
      CHARACTER*(*) name
      CHARACTER*120 nam1
      LOGICAL found

      INTEGER lench
      EXTERNAL lench

      IF(iiclib.NE.36) STOP '**** filopl: internal error (01) ****'

      IF(iicfil.NE.36) THEN
          DO 1 i=iunf1,iunf2
 1        allunt(i)=.false.
          iicfil=36
      END IF

      nam1=name
      INQUIRE(FILE=nam1,EXIST=found)
      IF(.NOT.found) THEN
          nam1=libdir(1:lenld)//name
          INQUIRE(FILE=nam1,EXIST=found)
          IF(.NOT.found) THEN
              ll=lench(name)
              WRITE(*,102) name(1:ll)
              STOP '**** filopl: abnormal end ****'
          END IF
      END IF
 102  FORMAT(' **** filopl: cannot find file "',a,'" ****')

      DO 2 iun=iunf1,iunf2
      IF(allunt(iun)) GOTO 2
      OPEN(iun,FILE=nam1,STATUS='old',ERR=3)
      filnam(iun)=nam1
      allunt(iun)=.true.
      RETURN
 2    CONTINUE

      STOP '**** filopl: all units are already allocated ****'

 3    CONTINUE
      ll=lench(nam1)
      WRITE(*,101) nam1(1:ll)
 101  FORMAT(' **** filopl: cannot OPEN file "',a,'" (STATUS=old) ****')
      STOP '**** filopl: abnormal end ****'
      END
* Copyright (C) 1996-1998 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: January 12, 1998
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         L I B I N I                           *
*  *                                                               *
*  *         Inizialization of default library directory           *
*  *                                                               *
*  *****************************************************************
*
      SUBROUTINE libini
      IMPLICIT NONE

      INCLUDE 'sysdep.h'
      INCLUDE 'parlib.h'

* Common blocks to be initialized:
      INCLUDE 'comlib.h'

      LOGICAL found
      INTEGER unit

      INTEGER lench
      EXTERNAL lench

* The file 'libdir.dat' can be used to modify the path of the
* library directory (with respect to the built-in value contained
* in parlib.h) without need of compiling again the software:
* it is searched only in the current (working) directory
      INQUIRE(FILE='libdir.dat',EXIST=found)

      IF(found) THEN
          CALL filopn(unit,'libdir.dat','old')
          READ(unit,100,END=10,ERR=10) libdir
          CALL filclo(unit,' ')
      ELSE
          libdir=dlibd
      END IF
 100  FORMAT(A)
      lenld=lench(libdir)
      IF(libdir(lenld:lenld).NE.dircha) THEN
          lenld=lenld+1
          libdir(lenld:lenld)=dircha
      END IF

      iiclib=36
      RETURN

 10   CONTINUE
      STOP '**** libini: error reading file "libdir.dat" ****'

      END
* Copyright (C) 1996 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: September 19, 1996
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         F I L O P F                           *
*  *                                                               *
*  *        Unit allocation and file opening (STATUS='old')        *
*  *         (searching also in default library directory)         *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    NAME      -  File name to be opened
*
* OUTPUT:   IUN       -  Allocated unit
*           FOUND     -  Was the file found?
*
      SUBROUTINE filopf(iun,name,found)
      IMPLICIT NONE

      INCLUDE 'comfil.h'

* NEEDED common blocks:
      INCLUDE 'comlib.h'

      INTEGER iun,ll,i
      CHARACTER*(*) name
      CHARACTER*120 tname
      LOGICAL found

      INTEGER lench
      EXTERNAL lench

      IF(iiclib.NE.36) STOP '**** filopf: internal error (01) ****'

      IF(iicfil.NE.36) THEN
          DO 1 i=iunf1,iunf2
 1        allunt(i)=.false.
          iicfil=36
      END IF

      tname=name
      INQUIRE(FILE=tname,EXIST=found)
      IF(.NOT.found) THEN
          tname=libdir(1:lenld)//name
          INQUIRE(FILE=tname,EXIST=found)
          IF(.NOT.found) THEN
              iun=0
              RETURN
          END IF
      END IF

      DO 2 iun=iunf1,iunf2
      IF(allunt(iun)) GOTO 2
      OPEN(iun,FILE=tname,STATUS='old',ERR=3)
      filnam(iun)=tname
      allunt(iun)=.true.
      RETURN
 2    CONTINUE

      STOP '**** filopf: all units are already allocated ****'

 3    CONTINUE
      ll=lench(tname)
      WRITE(*,101) tname(1:ll)
 101  FORMAT(' **** filopf: cannot OPEN file "',a,'" (STATUS=old) ****')
      STOP '**** filopf: abnormal end ****'
      END
* Copyright (C) 1996 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: August 27, 1996
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         D L I F E X                           *
*  *                                                               *
*  *                Delete a file (if it exists)                   *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    FILE      -  Name of the file to be deleted
*
      SUBROUTINE dlifex(file)
      IMPLICIT NONE

      CHARACTER*(*) file

      INTEGER unit
      LOGICAL found

      INQUIRE(FILE=file,EXIST=found)
      IF(found) THEN
          CALL filopn(unit,file,'old')
          CALL filclo(unit,'delete')
      END IF

      END
c =======================================
c FILNAM
c
c computes  file name in assigned directory
c ========================================
      SUBROUTINE filnam(eledir,astnam,suffix,file,le)
      IMPLICIT NONE

      CHARACTER*(*) eledir
      CHARACTER*(*) astnam
      CHARACTER*(*) suffix
      CHARACTER*(*) file
      INTEGER le
c==== end interface
      INTEGER ldir,l1,l2,l3
      PARAMETER (ldir=60)
      CHARACTER*200 tmpstr
      INTEGER lench
      INCLUDE 'sysdep.h'

      l1=lench(eledir)
      l2=lench(astnam)
      l3=lench(suffix)
      IF(l1+l2+l3+2.gt.200)THEN
         WRITE(*,*) 'filnam: Filename too long:',
     +        eledir(1:l1),dircha,astnam(1:l2),'.',suffix(1:l3)
         STOP
      ENDIF
      tmpstr=eledir(1:l1)//dircha//astnam(1:l2)//'.'//suffix(1:l3)
      CALL rmsp(tmpstr,le)
      file=tmpstr(1:le)
      IF(le.gt.ldir-4)THEN
         WRITE(*,*)'filnam: possible file name truncation',astnam,file
      ENDIF
      RETURN
      END

c ==============================================
c handling of name with multiple solutions
      SUBROUTINE splinam(name0,name1,m)
      IMPLICIT NONE
      CHARACTER*19 name0
      CHARACTER*9 name1
      INTEGER m
      INTEGER i
c find separator
      i=index(name0,'_')
      IF(i.le.0)THEN
c        WRITE(*,*)' name splitting problem for ', name0
         name1=name0(1:9)
         m=1
         RETURN
      ENDIF
      name1=' '
      name1=name0(1:i-1)
      READ(name0(i+1:i+4),101)m
c     WRITE(*,*)name1,m
 101  FORMAT(i4)
      RETURN
      END
c =======================================
c FIDINAM
c
c computes  file name in assigned directory
c ========================================
      SUBROUTINE fidinam(eledir,astnam,suffix,file,le)
      IMPLICIT NONE
      INTEGER ldir
      PARAMETER (ldir=60)
      CHARACTER*60 eledir
      CHARACTER*(*) astnam
      CHARACTER*9 nam0,namp
      CHARACTER*3 suffix
      CHARACTER*60 file
      INTEGER le,ld,ll,ld1
      INCLUDE 'sysdep.h'
      INTEGER j
      CHARACTER*10 dircom2,direct
      CALL rmsp(eledir,le)
c convert name in case it is of the form nam0=namp
      call rmsp(astnam,ld1)
      ll=index(astnam,'=')-1
      IF(ll.lt.0)THEN
         nam0=astnam(1:9)
      ELSE
         nam0=astnam(1:ll)
         namp=astnam(ll+2:ld1)
      ENDIF
      CALL rmsp(nam0,ld1)
      direct=dircom2(nam0,ld)
      file=eledir(1:le)//dircha//direct(1:ld)//dircha//astnam
      CALL rmsp(file,le)
      file=file(1:le)//'.'//suffix
      IF(le.gt.ldir-4)THEN
         WRITE(*,*)'fidinam: possible file name truncation ',astnam
         WRITE(*,*)file
      ENDIF
      CALL rmsp(file,le)
      DO j=le+1,60
        file(j:j)=' '
      ENDDO
      RETURN
      END

c =============================
c DIRCOM
c 
c computes subdirectory name
      CHARACTER*10 FUNCTION dircom2(name0,le)
      IMPLICIT NONE
      CHARACTER*9 name0
      INTEGER le
      INCLUDE 'proout.h'
      CHARACTER*9 name
      CHARACTER*1 nm
      INTEGER number,ll,ld,nyea,nodir
c make sure it is left aligned
      name=name0
      CALL rmsp(name,ll)
c find if numbered
      READ(name0,FMT='(I9)',ERR=10)number
c numbered: directories of 1000 each
      nodir=number/1000
      dircom2=' '
      IF(nodir.gt.9)THEN
         WRITE(dircom2,'(I2)')nodir
         le=2
      ELSE
         WRITE(dircom2,'(I1)')nodir
         le=1
      ENDIF
      RETURN
c unnumbered
 10   CONTINUE
c select surveys
      ld=index(name,'-')
      IF(ld.gt.0)THEN      
c surveys
         IF(index(name,'P-L').ne.0)THEN
c Palomar-Leiden
            dircom2='P-L'
            le=3
         ELSEIF(index(name,'T-').ne.0)THEN
c trojan surveys
            le=3
            IF(index(name,'T-1').ne.0)THEN
               dircom2='T-1'
            ELSEIF(index(name,'T-2').ne.0)THEN
               dircom2='T-2'
            ELSEIF(index(name,'T-3').ne.0)THEN
               dircom2='T-3'
            ELSE
               WRITE(*,*)'dircom2: name parsing error',name0
               WRITE(ierrou,*)'dircom2: name parsing error',name0
               numerr=numerr+1
               dircom2='unknown'
               le=7
            ENDIF
         ELSE
            WRITE(*,*)'dircom2: name parsing error',name0
            WRITE(ierrou,*)'dircom2: name parsing error',name0
            numerr=numerr+1
            dircom2='unknown'
            le=7
         ENDIF
      ELSE
c by year
         READ(name(1:4),FMT='(I4)',ERR=20) nyea
         IF(nyea.lt.1920)THEN
            dircom2='old'
            le=3
         ELSEIF(nyea.lt.1930)THEN
            dircom2='1920'
            le=4
         ELSEIF(nyea.lt.1940)THEN
            dircom2='1930'
            le=4
         ELSEIF(nyea.lt.1950)THEN
            dircom2='1940'
            le=4
         ELSEIF(nyea.lt.1960)THEN
            dircom2='1950'
            le=4
         ELSEIF(nyea.lt.1970)THEN
            dircom2='1960'
            le=4
         ELSEIF(nyea.lt.1980)THEN
            dircom2='1970'
            le=4
         ELSEIF(nyea.ge.1998)THEN
            nm=name(5:5)
            WRITE(dircom2,199)nyea,nm
 199        FORMAT(I4,A1)
            le=5
         ELSEIF(nyea.lt.0.or.nyea.gt.2100)THEN
            WRITE(*,*)'dircom2: name parsing error',name0
            WRITE(ierrou,*)'dircom2: name parsing error',name0
            numerr=numerr+1
            dircom2='unknown'
            le=7
         ELSE
            WRITE(dircom2,FMT='(I4)') nyea
            le=4
         ENDIF
      ENDIF
      RETURN
c error parsing year
 20   CONTINUE
      WRITE(*,*)'dircom2: name parsing error',name0
      WRITE(ierrou,*)'dircom2: name parsing error',name0
      numerr=numerr+1
      dircom2='unknown'
      le=7
      RETURN
      END

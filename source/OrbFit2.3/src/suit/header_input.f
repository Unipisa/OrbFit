c ===================MODULE header_input============================
c HEADER NAMELIST: Carpino 1996-1998
c CONTAINS:
c SUBROUTINES
c rdfnam		read and store header namelist
c rdfcha		read a character value
c rdfint		read an integer value
c rdflog		read a logical value
c rdfrea		read a real value
c rdftim		read a date/time value
c rdfref		read a reference system description
c *splkvc		split a record into keyword+value+comment
c *chkfln		check field length for input records
c *getrsc		get a data record skipping comment lines
c
c
c
* Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: March 8, 1997
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         R D F N A M                           *
*  *                                                               *
*  *            Read simplified file-header namelist               *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    UNIT      -  Input unit
*           FILE      -  Input filename (for error messages)
*
* OUTPUT:   NR        -  Records read so far
*
      SUBROUTINE rdfnam(unit,file,nr)
      IMPLICIT NONE

      INCLUDE 'parnam.h'
      INCLUDE 'parch.h'
      INCLUDE 'parcmc.h'

* Common blocks to be initialized:
      INCLUDE 'comfnm.h'

      INTEGER unit,nr
      INTEGER lf,lk,i
      CHARACTER*(lchx) rec,key1,val1,comm
      CHARACTER*(*) file
      LOGICAL skip,end

      INTEGER lench
      EXTERNAL lench

      iicfnm=0

      nfne=0
      kuorlf=0
      nmif=file
      hnfuni=unit
      lf=lench(file)
      nr=0

* Read a record from the namelist
 2    READ(unit,100,end=11) rec
 100  FORMAT(a)
      nr=nr+1
      CALL splkvc(rec,key1,val1,comm,skip,end)

* Detection of the end of the namelist
      IF(end) THEN
          iicfnm=36
          RETURN
      END IF

* Skip comment lines
      IF(skip) GOTO 2

* Look if the key is already present in the namelist
      lk=lench(key1)
      DO 4 i=1,nfne
      IF(key1(1:lk).EQ.keysf(i)) THEN
          WRITE(*,210) key1(1:lk),file(1:lf),krcfnm(i),nr
          STOP '**** rdfnam: abnormal end ****'
      END IF
 4    CONTINUE
 210  FORMAT(' rdfnam: duplicate definition of keyword "',a,
     +       '" in file "',a,'":'/
     +       '         first  occurence at line',i5/
     +       '         second occurence at line',i5)

      nfne=nfne+1
      CALL chkpdf(nfne,nfnex,'nfnex')
      i=nfne
      keysf(i)=key1
      valsf(i)=val1
      krcfnm(i)=nr
      kuorf(i)=0
      GOTO 2

 11   CONTINUE

* Abort when no "END_OF_HEADER" record is found
*     WRITE(*,200)'end',file(1:lf),nr
*     STOP '**** rdfnam: abnormal end ****'
*200  FORMAT(' rdfnam: cannot find namelist ',a,':'/
*    +       '         unexpected END of file "',a,'" after record',i6)

      iicfnm=36

      END
* Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: March 8, 1997
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         R D F C H A                           *
*  *                                                               *
*  *                  Read a character string                      *
*  *           from simplified file-header namelist                *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    UNIT      -  Input unit
*           KEY       -  Keyword
*           REQRD     -  If true, when the keyword is not found in the
*                        namelist, execution stops with an error message
*
* OUTPUT:   C         -  Character string
*           FOUND     -  True when the keyword has been found in the
*                        namelist
*           KR        -  Record number in the input file
*
      SUBROUTINE rdfcha(unit,key,reqrd,c,found,kr)
      IMPLICIT NONE

      INCLUDE 'parnam.h'

* NEEDED common blocks:
      INCLUDE 'comfnm.h'

      CHARACTER*(*) key,c
      LOGICAL reqrd,found
      INTEGER kr,unit

      INTEGER i,lk,lz,lf
      CHARACTER vartyp*30,rest*(lcvx)
      LOGICAL error

      INTEGER lench
      EXTERNAL lench

      IF(iicfnm.NE.36) STOP '**** rdfcha: internal error (01) ****'
      IF(unit.NE.hnfuni) STOP '**** rdfcha: internal error (02) ****'

      vartyp='CHARACTER'

      kr=0
      DO 1 i=1,nfne
      IF(key.EQ.keysf(i)) THEN
          found=.true.
          CALL strcnt(valsf(i),c,rest,error)
          IF(error) GOTO 10
          kr=krcfnm(i)
          kuorlf=kuorlf+1
          kuorf(i)=kuorlf
          RETURN
      END IF
 1    CONTINUE

      found=.false.
      IF(reqrd) THEN
          lk=lench(key)
          lz=lench(vartyp)
          lf=lench(nmif)
          WRITE(*,100) key(1:lk),vartyp(1:lz),nmif(1:lf)
          STOP '**** rdfcha: abnormal end ****'
      END IF
 100  FORMAT(' ERROR: missing definition of keyword "',a,'" (type: ',
     +       a,')'/8x,'in file "',a,'"')

      RETURN

* Error in reading namelist value
 10   CONTINUE
      lk=lench(key)
      lz=lench(vartyp)
      lf=lench(nmif)
      WRITE(*,101) key(1:lk),vartyp(1:lz),krcfnm(i),nmif(1:lf)
 101  FORMAT(' ERROR in reading value for keyword "',a,'" (type: ',
     +       a,')'/7x,'at record',i4,' in file "',a,'"')
      STOP '**** rdfcha: abnormal end ****'
      END
* Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: March 8, 1997
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         R D F I N T                           *
*  *                                                               *
*  *                  Read an integer quantity                     *
*  *           from simplified file-header namelist                *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    UNIT      -  Input unit
*           KEY       -  Keyword
*           REQRD     -  If true, when the keyword is not found in the
*                        namelist, execution stops with an error message
*
* OUTPUT:   K         -  Integer value
*           FOUND     -  True when the keyword has been found in the
*                        namelist
*           KR        -  Record number in the input file
*
      SUBROUTINE rdfint(unit,key,reqrd,k,found,kr)
      IMPLICIT NONE

      INCLUDE 'parnam.h'

* NEEDED common blocks:
      INCLUDE 'comfnm.h'

      CHARACTER*(*) key
      CHARACTER vartyp*30
      LOGICAL reqrd,found
      INTEGER k,kr,i,lk,lz,lf,unit


      INTEGER lench
      EXTERNAL lench

      IF(iicfnm.NE.36) STOP '**** rdfint: internal error (01) ****'
      IF(unit.NE.hnfuni) STOP '**** rdfint: internal error (02) ****'

      vartyp='INTEGER'

      kr=0
      DO 1 i=1,nfne
      IF(key.EQ.keysf(i)) THEN
          found=.true.
          READ(valsf(i),*,ERR=10) k
          kr=krcfnm(i)
          kuorlf=kuorlf+1
          kuorf(i)=kuorlf
          RETURN
      END IF
 1    CONTINUE

      found=.false.
      IF(reqrd) THEN
          lk=lench(key)
          lz=lench(vartyp)
          lf=lench(nmif)
          WRITE(*,100) key(1:lk),vartyp(1:lz),nmif(1:lf)
          STOP '**** rdfint: abnormal end ****'
      END IF
 100  FORMAT(' ERROR: missing definition of keyword "',a,'" (type: ',
     +       a,')'/8x,'in file "',a,'"')

      RETURN

* Error in reading namelist value
 10   CONTINUE
      lk=lench(key)
      lz=lench(vartyp)
      lf=lench(nmif)
      WRITE(*,101) key(1:lk),vartyp(1:lz),krcfnm(i),nmif(1:lf)
 101  FORMAT(' ERROR in reading value for keyword "',a,'" (type: ',
     +       a,')'/7x,'at record',i4,' in file "',a,'"')
      STOP '**** rdfint: abnormal end ****'
      END
* Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: March 8, 1997
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         R D F L O G                           *
*  *                                                               *
*  *                  Read a logical quantity                      *
*  *           from simplified file-header namelist                *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    UNIT      -  Input unit
*           KEY       -  Keyword
*           REQRD     -  If true, when the keyword is not found in the
*                        namelist, execution stops with an error message
*
* OUTPUT:   FLAG      -  Logical value
*           FOUND     -  True when the keyword has been found in the
*                        namelist
*           KR        -  Record number in the input file
*
      SUBROUTINE rdflog(unit,key,reqrd,flag,found,kr)
      IMPLICIT NONE

      INCLUDE 'parnam.h'

* NEEDED common blocks:
      INCLUDE 'comfnm.h'

      CHARACTER*(*) key
      CHARACTER vartyp*30
      LOGICAL reqrd,found,flag
      INTEGER kr,i,lk,lz,lf,unit


      INTEGER lench
      EXTERNAL lench

      IF(iicfnm.NE.36) STOP '**** rdflog: internal error (01) ****'
      IF(unit.NE.hnfuni) STOP '**** rdflog: internal error (02) ****'

      vartyp='LOGICAL'

      kr=0
      DO 1 i=1,nfne
      IF(key.EQ.keysf(i)) THEN
          found=.true.
          READ(valsf(i),*,ERR=10) flag
          kr=krcfnm(i)
          kuorlf=kuorlf+1
          kuorf(i)=kuorlf
          RETURN
      END IF
 1    CONTINUE

      found=.false.
      IF(reqrd) THEN
          lk=lench(key)
          lz=lench(vartyp)
          lf=lench(nmif)
          WRITE(*,100) key(1:lk),vartyp(1:lz),nmif(1:lf)
          STOP '**** rdflog: abnormal end ****'
      END IF
 100  FORMAT(' ERROR: missing definition of keyword "',a,'" (type: ',
     +       a,')'/8x,'in file "',a,'"')

      RETURN

* Error in reading namelist value
 10   CONTINUE
      lk=lench(key)
      lz=lench(vartyp)
      lf=lench(nmif)
      WRITE(*,101) key(1:lk),vartyp(1:lz),krcfnm(i),nmif(1:lf)
 101  FORMAT(' ERROR in reading value for keyword "',a,'" (type: ',
     +       a,')'/7x,'at record',i4,' in file "',a,'"')
      STOP '**** rdflog: abnormal end ****'
      END
* Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: March 8, 1997
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         R D F R E A                           *
*  *                                                               *
*  *                    Read a real quantity                       *
*  *           from simplified file-header namelist                *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    UNIT      -  Input unit
*           KEY       -  Keyword
*           REQRD     -  If true, when the keyword is not found in the
*                        namelist, execution stops with an error message
*
* OUTPUT:   V         -  Real value
*           FOUND     -  True when the keyword has been found in the
*                        namelist
*           KR        -  Record number in the input file
*
      SUBROUTINE rdfrea(unit,key,reqrd,v,found,kr)
      IMPLICIT NONE

      INCLUDE 'parnam.h'

* NEEDED common blocks:
      INCLUDE 'comfnm.h'

      CHARACTER*(*) key
      CHARACTER vartyp*30
      LOGICAL reqrd,found
      DOUBLE PRECISION v
      INTEGER kr,i,lk,lz,lf,unit


      INTEGER lench
      EXTERNAL lench

      IF(iicfnm.NE.36) STOP '**** rdfrea: internal error (01) ****'
      IF(unit.NE.hnfuni) STOP '**** rdfrea: internal error (02) ****'

      vartyp='REAL'

      kr=0
      DO 1 i=1,nfne
      IF(key.EQ.keysf(i)) THEN
          found=.true.
          READ(valsf(i),*,ERR=10) v
          kr=krcfnm(i)
          kuorlf=kuorlf+1
          kuorf(i)=kuorlf
          RETURN
      END IF
 1    CONTINUE

      found=.false.
      IF(reqrd) THEN
          lk=lench(key)
          lz=lench(vartyp)
          lf=lench(nmif)
          WRITE(*,100) key(1:lk),vartyp(1:lz),nmif(1:lf)
          STOP '**** rdfrea: abnormal end ****'
      END IF
 100  FORMAT(' ERROR: missing definition of keyword "',a,'" (type: ',
     +       a,')'/8x,'in file "',a,'"')

      RETURN

* Error in reading namelist value
 10   CONTINUE
      lk=lench(key)
      lz=lench(vartyp)
      lf=lench(nmif)
      WRITE(*,101) key(1:lk),vartyp(1:lz),krcfnm(i),nmif(1:lf)
 101  FORMAT(' ERROR in reading value for keyword "',a,'" (type: ',
     +       a,')'/7x,'at record',i4,' in file "',a,'"')
      STOP '**** rdfrea: abnormal end ****'
      END
* Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: March 8, 1997
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         R D F T I M                           *
*  *                                                               *
*  *                Read a time quantity (MJD+SEC)                 *
*  *           from simplified file-header namelist                *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    UNIT      -  Input unit
*           KEY       -  Keyword
*           REQRD     -  If true, when the keyword is not found in the
*                        namelist, execution stops with an error message
*
* OUTPUT:   TIMSTR    -  Time as a string
*           MJD       -  Modified Julian Date (integer part)
*           SEC       -  Seconds within the day
*           SCALE     -  Time scale
*           FOUND     -  True when the keyword has been found in the
*                        namelist
*           KR        -  Record number in the input file
*
      SUBROUTINE rdftim(unit,key,reqrd,timstr,mjd,sec,scale,found,kr)
      IMPLICIT NONE

      INCLUDE 'parnam.h'

* NEEDED common blocks:
      INCLUDE 'comfnm.h'

      CHARACTER*(*) key,timstr,scale
      CHARACTER*30 vartyp
      LOGICAL reqrd,found,error
      DOUBLE PRECISION sec
      INTEGER mjd,kr,i,lk,lz,lf,unit

      INTEGER lench
      EXTERNAL lench

      IF(iicfnm.NE.36) STOP '**** rdftim: internal error (01) ****'
      IF(unit.NE.hnfuni) STOP '**** rdftim: internal error (02) ****'

      vartyp='TIME'

      kr=0
      DO 1 i=1,nfne
      IF(key.EQ.keysf(i)) THEN
          found=.true.
          CALL ch2tim(valsf(i),mjd,sec,scale,error)
          IF(error) GOTO 10
          timstr=valsf(i)
          kr=krcfnm(i)
          kuorlf=kuorlf+1
          kuorf(i)=kuorlf
          RETURN
      END IF
 1    CONTINUE

      found=.false.
      IF(reqrd) THEN
          lk=lench(key)
          lz=lench(vartyp)
          lf=lench(nmif)
          WRITE(*,100) key(1:lk),vartyp(1:lz),nmif(1:lf)
          STOP '**** rdftim: abnormal end ****'
      END IF
 100  FORMAT(' ERROR: missing definition of keyword "',a,'" (type: ',
     +       a,')'/8x,'in file "',a,'"')

      RETURN

* Error in reading namelist value
 10   CONTINUE
      lk=lench(key)
      lz=lench(vartyp)
      lf=lench(nmif)
      WRITE(*,101) key(1:lk),vartyp(1:lz),krcfnm(i),nmif(1:lf)
 101  FORMAT(' ERROR in reading value for keyword "',a,'" (type: ',
     +       a,')'/7x,'at record',i4,' in file "',a,'"')
      STOP '**** rdftim: abnormal end ****'
      END
* Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: March 8, 1997
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         R D F R E F                           *
*  *                                                               *
*  *            Read a reference system description                *
*  *           from simplified file-header namelist                *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    UNIT      -  Input unit
*           KEY       -  Keyword
*           REQRD     -  If true, when the keyword is not found in the
*                        namelist, execution stops with an error message
*
* OUTPUT:   RSYS      -  Reference system type (EQUM/EQUT/ECLM)
*           EPOCH     -  Epoch specification (J2000/OFDATE)
*           FOUND     -  True when the keyword has been found in the
*                        namelist
*           KR        -  Record number in the input file
*
      SUBROUTINE rdfref(unit,key,reqrd,rsys,epoch,found,kr)
      IMPLICIT NONE

      INCLUDE 'parnam.h'

* NEEDED common blocks:
      INCLUDE 'comfnm.h'

      CHARACTER*(*) key,rsys,epoch
      LOGICAL reqrd,found
      INTEGER kr,unit

      INTEGER i,lk,lz,lf
      CHARACTER vartyp*30
      LOGICAL error

      INTEGER lench
      EXTERNAL lench

      IF(iicfnm.NE.36) STOP '**** rdfref: internal error (01) ****'
      IF(unit.NE.hnfuni) STOP '**** rdfref: internal error (02) ****'

      vartyp='REFERENCE SYSTEM'

      kr=0
      DO 1 i=1,nfne
      IF(key.EQ.keysf(i)) THEN
          found=.true.
          CALL ch2ref(valsf(i),rsys,epoch,error)
          IF(error) GOTO 10
          kr=krcfnm(i)
          kuorlf=kuorlf+1
          kuorf(i)=kuorlf
          RETURN
      END IF
 1    CONTINUE

      found=.false.
      IF(reqrd) THEN
          lk=lench(key)
          lz=lench(vartyp)
          lf=lench(nmif)
          WRITE(*,100) key(1:lk),vartyp(1:lz),nmif(1:lf)
          STOP '**** rdfref: abnormal end ****'
      END IF
 100  FORMAT(' ERROR: missing definition of keyword "',a,'" (type: ',
     +       a,')'/8x,'in file "',a,'"')

      rsys=' '
      epoch=' '

      RETURN

* Error in reading namelist value
 10   CONTINUE
      lk=lench(key)
      lz=lench(vartyp)
      lf=lench(nmif)
      WRITE(*,101) key(1:lk),vartyp(1:lz),krcfnm(i),nmif(1:lf)
 101  FORMAT(' ERROR in reading value for keyword "',a,'" (type: ',
     +       a,')'/7x,'at record',i4,' in file "',a,'"')
      STOP '**** rdfref: abnormal end ****'
      END
* Copyright (C) 1996 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: December 31, 1996
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         S P L K V C                           *
*  *                                                               *
*  *     Split a namelist record into keyword+value+comment        *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    REC       -  Record
*
* OUTPUT:   KEY       -  Keyword field
*           VAL       -  Value field
*           COMM      -  Comment field
*           SKIP      -  "Only comment" flag
*           END       -  "End of namelist" flag
*
      SUBROUTINE splkvc(rec,key,val,comm,skip,end)
      IMPLICIT NONE

      INCLUDE 'parcmc.h'

      CHARACTER*(*) rec,key,val,comm
      LOGICAL skip,end

      CHARACTER*200 rec1,tmp
      INTEGER lr,lk,lv,ipc,ipu

      INTEGER lench
      EXTERNAL lench

      key=' '
      val=' '
      comm=' '

* Detection of the end of the namelist
      rec1=rec
      CALL rmsp(rec1,lr)
      IF(rec1(1:lr).EQ.'END_OF_HEADER') THEN
          end=.true.
          skip=.false.
          RETURN
      ELSE
          end=.false.
      END IF

* Compute length excluding comments
      lr=lench(rec)
      IF(lr.LT.1) THEN
          skip=.true.
          RETURN
      END IF
      ipc=INDEX(rec(1:lr),comcha)
      IF(ipc.EQ.1) THEN
          skip=.true.
          comm=rec(2:)
          RETURN
      END IF

* Separate comment from keyword+value
      IF(ipc.EQ.0) THEN
          rec1=rec(1:lr)
      ELSE
          rec1=rec(1:ipc-1)
          comm=rec(ipc+1:)
          lr=lench(rec1)
          IF(lr.LT.1) THEN
              skip=.true.
              RETURN
          END IF
      END IF
      skip=.false.

* Keyword field
      ipu=INDEX(rec1(1:lr),'=')
      IF(ipu.EQ.0) THEN
          CALL rmsp(rec1,lk)
          IF(lk.GT.LEN(key)) STOP '**** splkvc: lk > LEN(key) ****'
          key=rec1(1:lk)
          RETURN
      END IF
      tmp=rec1(1:ipu-1)
      CALL rmsp(tmp,lk)
      IF(lk.GT.LEN(key)) STOP '**** splkvc: lk > LEN(key) ****'
      key=tmp(1:lk)

* Value field
      tmp=rec1(ipu+1:)
      CALL norstr(tmp,lv)
      IF(lv.GT.LEN(val)) STOP '**** splkvc: lv > LEN(val) ****'
      val=tmp(1:lv)

      END
* Copyright (C) 1996-1998 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: January 13, 1998
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         C H K F L N                           *
*  *                                                               *
*  *                     Check field length                        *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    NA        -  Actual value required for the parameter
*           NX        -  Value declared for the parameter
*           NAME      -  Parameter name
*           KR        -  Record number in input file
*           FILE      -  Input file name
*
      SUBROUTINE chkfln(na,nx,name,kr,file)
      IMPLICIT NONE

      INTEGER na,nx,kr,ll,lf
      CHARACTER*(*) name,file

      INTEGER lench
      EXTERNAL lench

      IF(na.LE.nx) RETURN

      ll=lench(name)
      lf=lench(file)
      WRITE(*,100) name(1:ll),nx,kr,file(1:lf),na
 100  FORMAT(' **** Insufficient PARAMETER definition ****'/
     + ' Sorry, the present version of the program does not',
     + ' allow'/
     + ' a length of the "',a,'" field longer than',i4,' characters'/
     + ' Please correct line',i4,' in file ',a/
     + ' (length of the field =',i4,' characters)')
      STOP '**** chkfln: abnormal end ****'
      END

* Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: March 8, 1997
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         G E T R S C                           *
*  *                                                               *
*  *         Get a data record skipping comment lines              *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    UNIT      -  Input unit
*           NR        -  Number of records read so far
*
* OUTPUT:   REC       -  Record
*           NR        -  Number of records read so far (updated)
*           END       -  End of file reached
*
      SUBROUTINE getrsc(unit,rec,nr,end)
      IMPLICIT NONE

      INCLUDE 'parcmc.h'

      INTEGER unit,nr
      CHARACTER*(*) rec
      LOGICAL end

      end=.false.

 1    CONTINUE
      READ(unit,100,END=2) rec
 100  FORMAT(A)
      nr=nr+1
      IF(rec(1:1).EQ.comcha) GOTO 1
      RETURN

 2    CONTINUE
      end=.true.
      rec=' '

      END

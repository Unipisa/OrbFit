* Copyright (C) 1997-1999 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: February 16, 1999
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         R D O P T O                           *
*  *                                                               *
*  *                   Read options for ORBFIT                     *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    RUN       -  Run name
*           NIFX      -  Physical dimension of arrays ELFT,ELF1
*           CNIFX     -  Real name of NIFX parameter
*
* OUTPUT:   OP        -  Selection of execution steps:
*                          1 = print residual/weight file (obsolete)
*                          2 = initial orbit determination
*                          3 = differential correction
*                          4 = identifications
*                          5 = ephemerides
*                          6 = magnitude determination
*           NAME      -  Object names
*           NAMEO     -  Names to be used for looking in orb. el. files
*           NAMOF     -  Names of observation files (w/o extension)
*           DIR       -  Directories where observations are
*           NOBJ      -  Number of objects
*           ELFT      -  Input orbital element files (for all objects)
*           NELFT     -  Number of input element files (for all objects)
*           ELF1      -  Input orbital element files (for each object)
*           NELF1     -  Number of input element files (for each object)
*           OEFILE    -  Output orbital element file (possibly blank)
*           OEPTIM    -  Epoch of output elements (MJD, TDT)
*           OEPSET    -  Flag stating that an output epoch is requested
*           OETYPE    -  Type of output elements (CAR/EQU/KEP/EQP)
*
      SUBROUTINE rdopto(run,op,name,nameo,namof,dir,nobj,elft,nelft,
     +                  elf1,nelf1,oefile,oeptim,oepset,oetype,
     +                  nifx,cnifx)
      IMPLICIT NONE

      INCLUDE 'sysdep.h'

      INTEGER op(10),nobj,nelft,nelf1(2),nifx
      CHARACTER*(*) run,dir(2),elf1(nifx,2)
      CHARACTER*(80)name(2),nameo(2),namof(2)
      CHARACTER*(*) elft(nifx),oefile,oetype,cnifx
      DOUBLE PRECISION oeptim
      LOGICAL oepset

      INTEGER i,lr,mjd,mjde,ln,ld
      DOUBLE PRECISION sec,sece
      LOGICAL found,fail1,fail
      CHARACTER cc*100,scale*10

      INTEGER lench
      EXTERNAL lench

* Required execution steps (by default, all are selected but ephemeris)
      DO 1 i=1,10
      op(i)=1
 1    CONTINUE
      op(5)=0
      fail=.false.
***   CALL rdnint('operations.','print_res',op(1),.false.,
***  +            found,fail1,fail)
      CALL rdnint('operations.','init_orbdet',op(2),.false.,
     +            found,fail1,fail)
      CALL rdnint('operations.','diffcor',op(3),.false.,
     +            found,fail1,fail)
      CALL rdnint('operations.','ident',op(4),.false.,
     +            found,fail1,fail)
      CALL rdnint('operations.','ephem',op(5),.false.,
     +            found,fail1,fail)
      CALL rdnint('operations.','magfit',op(6),.false.,
     +            found,fail1,fail)

* Input orbital element file (common to all objects)
      nelft=0
      CALL rdmcha('input_files.','incond',elft,nelft,nifx,
     +            cnifx,.false.,found,fail1,fail)

* Output orbital element file
      lr=lench(run)
      oefile=run(1:lr)//'.oel'
      CALL rdncha('output_files.','elem',oefile,.false.,
     +            found,fail1,fail)

      nobj=1
* First object
      name(1)=run
      CALL rdncha('object1.','name',name(1),.false.,
     +            found,fail1,fail)
      IF(.NOT.fail1) THEN
          dir(1)=' '
          CALL rdncha('object1.','obs_dir',dir(1),.false.,
     +                found,fail1,fail)
      END IF
      nelf1(1)=0
      CALL rdmcha('object1.','inc_files',elf1(1,1),nelf1(1),nifx,
     +            cnifx,.false.,found,fail1,fail)
      IF(.NOT.found) CALL rdodin(name(1),nelf1(1),elf1(1,1),nifx)
      CALL rdncha('object1.','inc_name',nameo(1),.false.,
     +            found,fail1,fail)
      IF(.NOT.found) nameo(1)=name(1)
      CALL rmsp(nameo(1),ln)
      CALL rdncha('object1.','obs_fname',namof(1),.false.,
     +            found,fail1,fail)
      IF(.NOT.found) THEN
          namof(1)=name(1)
          CALL rmsp(namof(1),ln)
      END IF

* Second object
      name(2)=' '
      nameo(2)=' '
      namof(2)=' '
      CALL rdncha('object2.','name',name(2),.false.,
     +            found,fail1,fail)
      IF(found) THEN
          nobj=2
          dir(2)=' '
          CALL rdncha('object2.','obs_dir',dir(2),.false.,
     +                found,fail1,fail)
          nelf1(2)=0
          CALL rdmcha('object2.','inc_files',elf1(1,2),nelf1(2),nifx,
     +                cnifx,.false.,found,fail1,fail)
          IF(.NOT.found) CALL rdodin(name(2),nelf1(2),elf1(1,2),nifx)
          CALL rdncha('object2.','inc_name',nameo(2),.false.,
     +                found,fail1,fail)
          IF(.NOT.found) nameo(2)=name(2)
          CALL rmsp(nameo(2),ln)
          CALL rdncha('object2.','obs_fname',namof(2),.false.,
     +                found,fail1,fail)
          IF(.NOT.found) THEN
              namof(2)=name(2)
              CALL rmsp(namof(2),ln)
          END IF
      END IF

* Epoch and type of output elements
      oeptim=0.d0
      CALL rdntim('output.','epoch',cc,mjd,sec,scale,.false.,
     +            oepset,fail1,fail)
      IF(oepset) THEN
          CALL cnvtim(mjd,sec,scale,mjde,sece,'TDT')
          oeptim=mjde+sece/86400.d0
      END IF

      oetype='EQU'
      CALL rdncha('output.','elements',oetype,.false.,
     +            found,fail1,fail)

      IF(fail) STOP '**** rdopto: abnormal end ****'

      DO 2 i=1,nobj
      ld=lench(dir(i))
      IF(ld.LE.0) GOTO 2
      IF(dir(i)(ld:ld).NE.dircha) dir(i)(ld+1:ld+1)=dircha
 2    CONTINUE

      END
*
* DEFAULT NAMES FOR INPUT ORBITAL ELEMENT FILES
*
      SUBROUTINE rdodin(name,nelf,elf,nifx)
      IMPLICIT NONE

      INTEGER nelf,nifx
      CHARACTER*(*) name,elf(nifx)

      INTEGER ln
      CHARACTER*80 file
      LOGICAL found

      INTEGER lench
      EXTERNAL lench

      IF(nelf.NE.0) RETURN

      ln=lench(name)
      file=name(1:ln)//'.ele'
      INQUIRE(FILE=file,EXIST=found)
      IF(found .AND. nelf.LT.nifx) THEN
          nelf=nelf+1
          elf(nelf)=file
      END IF

      file='astorb.dat'
      file=name(1:ln)//'.ele'
      INQUIRE(FILE=file,EXIST=found)
      IF(found .AND. nelf.LT.nifx) THEN
          nelf=nelf+1
          elf(nelf)=file
      END IF

      END

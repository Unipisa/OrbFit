PROGRAM prop_tax
 USE propel_mod
 USE name_rules, ONLY: name_len
 IMPLICIT NONE
 CHARACTER*6 progna ! option control: progname
 CHARACTER*80 run ! option control: current run 
 INTEGER ler ! character counters
 INTEGER in_pro ! input units
 INTEGER iun_log,nearout ! output units
 INTEGER err_line, npr
 INTEGER n1, j, i, i1, i2
 DOUBLE PRECISION d_max ! control
 DOUBLE PRECISION pro(3), hmag
 CHARACTER*(name_len) name
! **********BEGIN EXEXCUTION******************************** 
! progna='protax'
 WRITE(*,*)' run name?'
 READ(*,*) run 
 CALL rmsp(run,ler)
 WRITE(*,*)' control d_max?'
 READ(*,*) d_max 
 WRITE(*,*)' control ', d_max 
! open log file
 CALL filopn(iun_log,run(1:ler)//'.prta','unknown')
 WRITE(iun_log,*)' control ', d_max 
! open input file 
 CALL filopn(in_pro,'numball.pro','old')
! read proper elements
 npr=0
 CALL  input_propels(in_pro,npr,err_line)
 CALL filclo(in_pro,' ')
 WRITE(iun_log,*)' input ',npr,' prop. els of numbered'
! open input file 
 CALL filopn(in_pro,'multiall.pro','old')
! read proper elements
 CALL  input_propels(in_pro,npr,err_line)
 WRITE(iun_log,*)' input ',npr,' prop. els of numbered and multiopp'
 CALL filclo(in_pro,' ')
! error return?
 IF(err_line.gt.0)THEN
    STOP
 ENDIF
! open output file 
 CALL filopn(nearout,run(1:ler)//'.near','unknown')
 WRITE(nearout,200)
200 FORMAT('  name   ',1X,'  H  ',1X,'  name   ',1X,'  H  ',1X,'  sigma   ',1X,&
&      '   da/a   ',1X,'    de    ',1X,'   dsinI   ')
!  FORMAT(A9,1X,F5.2,2X,A9,1X,F5.2,1X,F10.7,3(1X,F10.7))
! sort by sinI
 CALL heapsort(propel(1:npro)%pr_el(3),npro,isrti)
 WRITE(iun_log,*)' end sorting by sinI'
! main loop searching for neighbours
 n1=0
! select nearby cases
 DO j=1,npro-1
    pro=propel(j)%pr_el
    name=propel(j)%name
    hmag=propel(j)%hmag
    CALL taxostep(j,pro,name,hmag,d_max,nearout,n1)
    IF(mod(j,1000).eq.0) WRITE(iun_log,*)' tried ', j, ' objects, couples ',n1 
 ENDDO
 WRITE(iun_log,*)' end selection of close couples, n1=',n1
! sort by sigma
 CALL heapsort(clos(1:n1)%sigma,n1,isrts)
! write .near file
 DO j=1,n1
    i=isrts(j)
    i1=clos(i)%addr(1)
    i2=clos(i)%addr(2)
    WRITE(nearout,100)clos(i)%names(1),propel(i1)%hmag,clos(i)%names(2),propel(i2)%hmag,&
&           clos(i)%sigma,clos(i)%delta_el(1),clos(i)%delta_el(2),clos(i)%delta_el(3)
100    FORMAT(A9,1X,F5.2,2X,A9,1X,F5.2,1X,F10.7,3(1X,F10.7))
 ENDDO
 CALL filclo(nearout,' ')
! summary
 WRITE(*,*)' found ', n1,' couples closer than ', d_max
 WRITE(iun_log,*)' found ', n1,' couples closer than ', d_max
END PROGRAM prop_tax

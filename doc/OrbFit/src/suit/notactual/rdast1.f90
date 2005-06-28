! Copyright (C) 1997-1999 by Mario Carpino (carpino@brera.mi.astro.it)  
! Version: February 12, 1999                                            
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                         R D A S T 1                           *    
!  *                                                               *    
!  *          Read orbital elements for a list of objects          *    
!  *       from a file written in Bowell's astorb.dat format       *    
!  *          (version reading OLD, pre-1999 format)               *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! INPUT:    UNIT      -  Input FORTRAN unit (must be already opened)    
!           FILNAM    -  Input file name (for error messages)           
!           OBJNAM    -  Object names                                   
!           NOBJ      -  Number of objects                              
!                                                                       
! OUTPUT:   DEFORB    -  Tells whether orbital elements are defined     
!           DEFCN     -  Tells whether covariance/normal matrices       
!                            are defined                                
!           ELTYPE    -  Type of orbital elements (EQU/KEP/CAR)         
!           TELEM     -  Epoch of orbital elements (MJD, TDT)           
!           ELEM      -  Orbital elements (ECLM J2000)                  
!           COVE      -  Covariance matrix of orbital elements          
!           NORE      -  Normal matrix of orbital elements              
!           MASS      -  Mass (solar masses)                            
!           H         -  H absolute magnitude (if <-100, missing)       
!           G         -  G slope parameter                              
!           COMELE    -  Comment on orbital elements                    
!                                                                       
! WARNING: the routine assumes that objects having DEFORB=.true.        
!          have already orbital elements defined (possibly from another 
!          input file) and does not overwrite them                      
!                                                                       
! OBJECT NAME TRANSLATION: all names of objects appearing in the        
! file are modified (removing all blanks) before comparison with        
! the name requested by the calling module                              
!                                                                       
SUBROUTINE rdast1(unit,filnam,objnam,nobj,deforb,defcn,           &
     &                  eltype,telem,elem,cove,nore,mass,h,g,comele)
  USE fund_const    
  IMPLICIT NONE                                                                       
  INTEGER unit,nobj 
  DOUBLE PRECISION telem(nobj),elem(6,nobj),cove(6,6,nobj) 
  DOUBLE PRECISION nore(6,6,nobj),mass(nobj),h(nobj),g(nobj) 
  CHARACTER*(*) filnam,objnam(nobj),eltype(nobj),comele(nobj) 
  LOGICAL deforb(nobj),defcn(nobj)                                                                       
  INTEGER ln,nr,lf,nrem,k,flags(6),year,month,day,lc 
  DOUBLE PRECISION el1(6) 
  CHARACTER n1*4,n2*18,name*18,hc*5,gc*5,krc*10                                                   
  INTEGER lench 
  DOUBLE PRECISION tjm1 
  EXTERNAL lench,tjm1                                                                 
! Number of remaining object (orbit not yet found)                      
  nrem=0 
  DO 10 k=1,nobj 
     IF(deforb(k)) GOTO 10 
     nrem=nrem+1 
10 END DO
  IF(nrem.LE.0) RETURN 
  lf=lench(filnam)                                                                       
  nr=0 
1 CONTINUE 
  READ(unit,100,END=2) n1,n2 
  nr=nr+1 
  IF(n1.EQ.'    ') THEN 
     name=n2 
  ELSE 
     name=n1 
  END IF
  CALL rmsp(name,ln) 
  IF(ln.LE.0) THEN 
     WRITE(*,200) filnam(1:lf),nr 
200  FORMAT('ERROR in reading file "',A,'": no object name at record',I6)
     GOTO 1 
  END IF
  DO 3 k=1,nobj 
     IF(deforb(k)) GOTO 3 
     IF(name.EQ.objnam(k)) THEN 
        BACKSPACE(unit) 
        READ(unit,100) n1,n2,hc,gc,flags,year,month,day,el1
100     FORMAT(A4,1X,A18,17X,A5,1X,A5,17X,6I4,12X,I4,2I2,1X,   &
     &       3(F10.6,1X),F9.6,1X,F10.8,1X,F12.8)       
        deforb(k)=.true. 
        defcn(k)=.false. 
        eltype(k)='KEP' 
        telem(k)=tjm1(day,month,year,0.D0) 
        elem(1,k)=el1(6) 
        elem(2,k)=el1(5) 
        elem(3,k)=el1(4)*radeg 
        elem(4,k)=el1(3)*radeg 
        elem(5,k)=el1(2)*radeg 
        elem(6,k)=el1(1)*radeg 
        mass(k)=0.d0 
        IF(hc.EQ.'     ') THEN 
           h(k)=-1.D9 
        ELSE 
           READ(hc,101) h(k) 
        END IF
        IF(gc.EQ.'     ') THEN 
           g(k)=0.15D0 
        ELSE 
           READ(gc,101) g(k) 
101        FORMAT(F5.2) 
        END IF
        WRITE(krc,107) nr 
107     FORMAT(I6)
        CALL rmsp(krc,lc) 
        comele(k)='read from file "'//filnam(1:lf)//                  &
     &              '" at record '//krc(1:lc) 
        nrem=nrem-1 
        IF(nrem.LE.0) RETURN 
     END IF
3 END DO                                                                      
  GOTO 1 
2 CONTINUE
END SUBROUTINE rdast1

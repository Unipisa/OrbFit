C                                                                               
C      ASC2EPH creates a direct access JPL Planetary Ephemeris file from        
C      one or more ascii text files.                                            
C                                                                               
C      This program, 'asc2eph', requires (via standard input) an ascii          
C      header file ('header.XXX'), followed by one or more ascii ephemeris      
C      data files ('ascSYYYY.XXX').  All files must have the same ephemeris     
C      number, XXX.  Further, the data files must be consecutive in time        
C      with no gaps between them.                                               
C                                                                               
C      By default, the output ephemeris will span the same interval as the input
C      text file(s).  If you are interested in only a portion of data, set the  
C      below T1 and T2 to the begin and end times of the span you desire.  T1   
C      and T2 must be specified in  Julian Ephemeris Days (ET).                 
C                                                                               
C      A sample sequence of files might be:                                     
C                                                                               
C        header.405  asc+1920.405 asc+1940.405 asc+1960.405 asc+1980.405        
C                                                                               
C      (The data files for DE200 and DE405 contain 20 years each; for DE406, 50 
C                                                                               
C                                                                               
C      This program is written in standard Fortran-77.                          
C                                                                               
C ******************************************************************************
C                                                                               
C                                    *** NOTE ***                               
C                                                                               
C      However, the units in which the length of a direct access record is speci
C      are PROCESSOR DEPENDENT.  The parameter NRECL, the number of units per wo
C      controls the length of a record in the direct access ephemeris.          
C                                                                               
C *****  The user must choose one of the following statements  *****            
C                                                                               
       PARAMETER ( NRECL = 4 )                                                  
*       PARAMETER ( NRECL = 1 )                                                  
                                                                                
C ******************************************************************************
                                                                                
       CHARACTER*6  CNAM (400)                                                  
       CHARACTER*6  TTL  (14,3)                                                 
       CHARACTER*12 HEADER                                                      
                                                                                
       DOUBLE PRECISION  AU                                                     
       DOUBLE PRECISION  CVAL (400)                                             
       DOUBLE PRECISION  DB(3000), DB2Z
       DOUBLE PRECISION  EMRAT                                                  
       DOUBLE PRECISION  SS   (3)                                               
       DOUBLE PRECISION  T1, T2                                                 
                                                                                
       INTEGER  I, IPT (3,12), IRECSZ, IN                                       
       INTEGER  J                                                               
       INTEGER  K,KSIZE                                                         
       INTEGER  LPT (3)                                                         
       INTEGER  N, NCON, NROUT, NCOEFF, NRW, NUMDE                              
       INTEGER  OUT                                                             
                                                                                
       LOGICAL  FIRST                                                           
       DATA     DB2Z  /0.d0/
       DATA     FIRST / .TRUE. /                                                
                                                                                
C ***********************************************************************       
C     By default, the output ephemeris will span the same interval as the       
C     input ascii data file(s).  The user may reset these to other JED's.       
C                                                                               
       T1 = 0.D0                                                                
       T2 = 9999999.d0                                                          
                                                                                
                                                                                
      if(nrecl .eq. 0) then                                                     
      write(*,*)'*** ERROR: User did not set NRECL ***'                         
      stop                                                                      
      endif                                                                     
                                                                                
                                                                                
C      Write a fingerprint to the screen.                                       
                                                                                
       WRITE(*,*) ' JPL ASCII-TO-DIRECT-I/O program. ' //                       
     .            ' Last modified 31-Mar-1993.'                                 
                                                                                
                                                                                
C                                                                               
C      Read the size and number of main ephemeris records.                      
C                                                                               
       READ  (*,'(6X,I6)')  KSIZE                                               
       WRITE (*,100)              KSIZE                                         
100    FORMAT(/'KSIZE =',I6)                                                    
                                                                                
       IRECSZ = NRECL * KSIZE                                                   
                                                                                
                                                                                
C                                                                               
C      Now for the alphameric heading records (GROUP 1010)                      
C                                                                               
       CALL  NXTGRP ( HEADER )                                                  
                                                                                
       IF (HEADER .NE. 'GROUP   1010')  THEN                                    
          CALL  ERRPRT ( 1010, 'NOT HEADER' )                                   
       ENDIF                                                                    
                                                                                
       READ (*,'(14A6)')    TTL                                                 
       WRITE(*,'(/(14A6))') TTL                                                 
                                                                                
                                                                                
C                                                                               
C      Read start, end and record span  (GROUP 1030)                            
C                                                                               
       CALL  NXTGRP ( HEADER )                                                  
                                                                                
       IF ( HEADER .NE. 'GROUP   1030' ) THEN                                   
          CALL ERRPRT ( 1030, 'NOT HEADER' )                                    
       ENDIF                                                                    
                                                                                
       READ (*,'(3D12.0)')  SS                                                  
                                                                                
                                                                                
C                                                                               
C      Read number of constants and names of constants (GROUP 1040/4).          
C                                                                               
       CALL  NXTGRP ( HEADER )                                                  
                                                                                
       IF ( HEADER .NE. 'GROUP   1040 ') THEN                                   
          CALL ERRPRT ( 1040, 'NOT HEADER' )                                    
       ENDIF                                                                    
                                                                                
       READ (*,'(I6)')    N                                                     
       READ (*,'(10A8)') (CNAM(I),I=1,N)                                        
                                                                                
       NCON = N                                                                 
                                                                                
                                                                                
C                                                                               
C      Read number of values and values (GROUP 1041/4)                          
C                                                                               
       CALL  NXTGRP ( HEADER )                                                  
                                                                                
       IF ( HEADER .NE. 'GROUP   1041' ) THEN                                   
          CALL ERRPRT ( 1041, 'NOT HEADER' )                                    
       ENDIF                                                                    
                                                                                
       READ (*,'(I6)')       N                                                  
       READ (*,'(3D26.18)')  (CVAL(I),I=1,N)                                    
                                                                                
       DO  I = 1, N                                                             
           IF ( CNAM(I) .EQ. 'AU    ' )  AU    = CVAL(I)                        
           IF ( CNAM(I) .EQ. 'EMRAT ' )  EMRAT = CVAL(I)                        
           IF ( CNAM(I) .EQ. 'DENUM ' )  NUMDE = CVAL(I)                        
       END DO                                                                   
                                                                                
       WRITE (*,'(500(/2(A8,D24.16)))') (CNAM(I),CVAL(I),I=1,N)                 
                                                                                
                                                                                
C                                                                               
C      Read pointers needed by INTERP (GROUP 1050)                              
C                                                                               
       CALL  NXTGRP ( HEADER )                                                  
                                                                                
       IF ( HEADER .NE. 'GROUP   1050' ) THEN                                   
          CALL ERRPRT ( 1050, 'NOT HEADER' )                                    
       ENDIF                                                                    
                                                                                
       READ (*,'(13I6)')   ((IPT(I,J),J=1,12),LPT(I),I=1,3)                     
       WRITE(*,'(/(3I5))') IPT, LPT                                             
                                                                                
                                                                                
C                                                                               
C      Open direct-access output file ('JPLEPH')                                
C                                                                               
       OPEN ( UNIT   = 12,                                                      
     .        FILE   = 'JPLEPH',                                                
     .        ACCESS = 'DIRECT',                                                
     .        FORM   = 'UNFORMATTED',                                           
     .        RECL   = IRECSZ,                                                  
     .        STATUS = 'NEW' )                                                  
                                                                                
                                                                                
C                                                                               
C     Read and write the ephemeris data records (GROUP 1070).                   
C                                                                               
      CALL  NXTGRP ( HEADER )                                                   
                                                                                
      IF ( HEADER .NE. 'GROUP   1070' ) CALL ERRPRT(1070,'NOT HEADER')          
                                                                                
      NROUT  = 0                                                                
      IN     = 0                                                                
      OUT    = 0                                                                
                                                                                
      READ (*,'(2I6,3000(/3D26.18))',IOSTAT =IN)                                
     .      NRW, NCOEFF, (DB(K),K=1,NCOEFF)                                     
                                                                                
                                                                                
      DO WHILE (       ( IN    .EQ. 0 )                                         
     .           .AND. ( DB(2) .LT. T2) )                                       
                                                                                
          IF ( 2*NCOEFF .NE. KSIZE ) THEN                                       
             CALL ERRPRT(NCOEFF,' 2*NCOEFF not equal to KSIZE')                 
          ENDIF                                                                 
                                                                                
C                                                                               
C         Skip this data block if the end of the interval is less               
C         than the specified start time or if the it does not begin             
C         where the previous block ended.                                       
C                                                                               
          IF  ( (DB(2) .GE. T1) .AND. (DB(1) .GE. DB2Z) ) THEN                  
                                                                                
             IF ( FIRST ) THEN                                                  
C                                                                               
C               Don't worry about the intervals overlapping                     
C               or abutting if this is the first applicable                     
C               interval.                                                       
C                                                                               
                DB2Z  = DB(1)                                                   
                FIRST = .FALSE.                                                 
             ENDIF                                                              
                                                                                
             IF (DB(1) .NE. DB2Z ) THEN                                         
C                                                                               
C               Beginning of current interval is past the end                   
C               of the previous one.                                            
C                                                                               
                CALL ERRPRT (NRW, 'Records do not overlap or abut')             
             ENDIF                                                              
                                                                                
             DB2Z  = DB(2)                                                      
             NROUT = NROUT + 1                                                  
                                                                                
             WRITE (12,REC=NROUT+2,IOSTAT=OUT) (DB(K),K=1,NCOEFF)               
                                                                                
             IF ( OUT .NE. 0 ) THEN                                             
                CALL ERRPRT (NROUT,                                             
     .                     'th record not written because of error')            
             ENDIF                                                              
                                                                                
C                                                                               
C            Save this block's starting date, its interval span, and its end    
C            date.                                                              
C                                                                               
             IF (NROUT .EQ. 1) THEN                                             
                SS(1) = DB(1)                                                   
                SS(3) = DB(2) - DB(1)                                           
             ENDIF                                                              
             SS(2) = DB(2)                                                      
                                                                                
C                                                                               
C            Update the user as to our progress every 10th block.               
C                                                                               
             IF ( MOD(NROUT,10) .EQ. 1 ) THEN                                   
                IF ( DB(1) .GE. T1 ) THEN                                       
                   WRITE (*,271) NROUT, DB(2)                                   
                ELSE                                                            
                   WRITE (*,*)                                                  
     .           ' Searching for first requested record...'                     
                ENDIF                                                           
             ENDIF                                                              
271          FORMAT (I6,                                                        
     .              ' EPHEMERIS RECORDS WRITTEN.  LAST JED = ',                 
     .              F12.2)                                                      
                                                                                
          ENDIF                                                                 
                                                                                
          READ (*,'(2I6,3000(/3D26.18))',IOSTAT =IN)                            
     .          NRW, NCOEFF, (DB(K),K=1,NCOEFF)                                 
                                                                                
       END DO                                                                   
                                                                                
       WRITE (*,275) NROUT, DB(2)                                               
275    FORMAT(/I6, ' EPHEMERIS RECORDS WRITTEN.  LAST JED = ',F12.2)            
                                                                                
                                                                                
C                                                                               
C      Write header records onto output file.                                   
C                                                                               
       NROUT = 1                                                                
       WRITE(12,REC=1,IOSTAT=OUT) TTL,CNAM,SS,NCON,AU,                          
     .                            EMRAT,IPT,NUMDE,LPT                           
       IF ( OUT .NE. 0 ) THEN                                                   
          CALL ERRPRT ( NROUT, 'th record not written because of error')        
       ENDIF                                                                    
                                                                                
       NROUT = 2                                                                
       WRITE(12,REC=2,IOSTAT=OUT)CVAL                                           
       IF ( OUT .NE. 0 ) THEN                                                   
          CALL ERRPRT ( NROUT, 'th record not written because of error')        
       ENDIF                                                                    
                                                                                
                                                                                
C                                                                               
C      We're through.  Wrap it up.                                              
C                                                                               
       CLOSE (12)                                                               
       STOP ' OK'                                                               
                                                                                
                                                                                
       END                                                                      
                                                                                
C      $Header: asc2eph.oth,v 1.1 93/08/23 14:53:39 faith Exp $                 
C                                                                               
       SUBROUTINE  ERRPRT (I, MSG)                                              
                                                                                
       CHARACTER*(*)  MSG                                                       
       INTEGER        I                                                         
                                                                                
       WRITE (*,200)  I, MSG                                                    
200    FORMAT('ERROR #',I8,2X,A50)                                              
                                                                                
       STOP ' ERROR '                                                           
       END                                                                      
                                                                                
                                                                                
      SUBROUTINE  NXTGRP ( HEADER )                                             
                                                                                
      CHARACTER*(*)  HEADER                                                     
      CHARACTER*12   BLANK                                                      
                                                                                
C                                                                               
C     Start with nothing.                                                       
C                                                                               
      HEADER = ' '                                                              
                                                                                
C                                                                               
C     The next non-blank line we encounter is a header record.                  
C     The group header and data are seperated by a blank line.                  
C                                                                               
      DO WHILE ( HEADER .EQ. ' ' )                                              
          READ (*,'(A)') HEADER                                                 
      ENDDO                                                                     
                                                                                
C     Found the header.  Read the blank line so we can get at the data.         
      READ (*, '(A)') BLANK                                                     
                                                                                
      RETURN                                                                    
      END                                                                       
                                                                                

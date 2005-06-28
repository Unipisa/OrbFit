! ===================================================================   
! TRIVOPT minimum option routine, without propagation!!!
! ===================================================================   
! input options, for the propagator and the specific main program       
! input: progna = program name (6 characters)                           
!        run    = run identifier (80 characters, up to 76 non-blank)    
      SUBROUTINE trivopt(progna,run,iun20) 
      implicit none 
      character*6 progna 
      character*80 run 
      integer iun20 
! ==========END INTERFACE============================================   
      integer le,iunit 
      character*100 file 
      logical found 
! =============================                                         
! read option for propagator                                            
      CALL libini 
      CALL namini 
! default options are only for propag                                   
      CALL filopl(iunit,'propag.def') 
      CALL rdnam(iunit) 
      CALL filclo(iunit,' ') 
! =============================                                         
! particular options for this run                                       
      file=run//'.mop' 
      CALL rmsp(file,le) 
      INQUIRE(FILE=file(1:le),EXIST=found) 
      IF(found) THEN 
        CALL filopn(iunit,file(1:le),'OLD') 
        CALL rdnam(iunit) 
        CALL filclo(iunit,' ') 
      ELSE 
        write(*,*)'**** file not found: ',file 
        write(*,*)'******* ',progna,' only default ****' 
      ENDIF 
! =============================                                         
! check for non-existing options                                        
      call rmsp(progna,le) 
      file=progna(1:le)//'.key' 
      CALL rdklst(file) 
      CALL chkkey                                                      
! ===================================================================== 
! Output files: for control and results, for covariance                 
      CALL rmsp(run,le) 
      file=run(1:le)//'.mou' 
      call filopn(iun20,file,'UNKNOWN') 
! =============================                                         
! read option for physical model and integration method                 
!     CALL rmodel                                                  
! ====================================================                  
      END SUBROUTINE trivopt                                          

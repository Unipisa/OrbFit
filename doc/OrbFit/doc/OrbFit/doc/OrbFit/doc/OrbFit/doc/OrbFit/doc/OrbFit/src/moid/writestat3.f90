! ************************************************************      
! ** C O M P U T I N G   S T A T I O N A R Y    P O I N T S **      
! ***  and   W R I T I N G   on   O U T P U T    F I L E S ***      
! ************************************************************      
! ********** written by GIOVANNI F. GRONCHI (2001) ***********      
! ******* Department of Mathematics, UNIVERSITY of PISA ******      
! last modified June 2004 (GFG)
! ============================================================      
  SUBROUTINE writestat3(name,lnam,iunit,iun1,iun2,iun3,iun4,iun5, &
       & ekpl,elkep)                                                   
    USE output_control
    IMPLICIT NONE 
! asteroid name and name lenght                                     
    CHARACTER*9,INTENT(IN) ::  name 
    INTEGER,INTENT(IN) :: lnam 
! file identifiers                                                  
    INTEGER,INTENT(IN) :: iunit,iun1,iun2,iun3,iun4,iun5 
! elements of the Earth and of the asteroid                         
    DOUBLE PRECISION,INTENT(IN) :: ekpl(6),elkep(6) 
! ------------- end interface --------------------------------
! eccentric anomalies                                               
    DOUBLE PRECISION u(20),upl(20) 
! SQUARED DISTANCE function                                         
    DOUBLE PRECISION D2(20) 
!
    INTEGER nstat ! number of stationary points (maximum = 20)   
    INTEGER nummin,nummax ! number of relative minima/maxima found
    INTEGER :: answer(20) ! type of singular point : 
!                             answer(j) =  1    MAXIMUM
!                             answer(j) = -1    MINIMUM                
!                             answer(j) =  0    SADDLE                 
!                             answer(j) =  2    CANNOT DECIDE          
    CHARACTER*7 :: singtype 

    LOGICAL :: morse ! check with Morse theory 
!                        morse = true      OK   
!                        morse = false     ERROR  
    LOGICAL :: weier ! check with Weierstrass theory:
!                        weier = true      OK 
!                        weier = false     ERROR
    LOGICAL :: warnflag(3) ! program warning flags:
!      warnflag(1) = true    OK
!      warnflag(1) = false   leading coefficient of resultant is very small
!      warnflag(2) = true    OK                                          
!      warnflag(2) = false   higher degree terms in resultant are not small
!      warnflag(3) = true    OK                                           
!      warnflag(3) = false   low precision for resultant coefficients     
    LOGICAL :: sflag(6) ! solving system messages:
!         sflag(1) = true    OK
!         sflag(1) = false   there are two good solutions             
!         sflag(2) = true    OK                                               
!         sflag(2) = false   neither of the two evaluation is close to 0
!         sflag(3) = true    OK                                               
!         sflag(3) = false   the s-component of the solution is complex
!         sflag(4) = true    OK    
!         sflag(4) = false   leading coeff. of the 2nd degree pol. is small
!         sflag(5) = true    OK    
!         sflag(5) = false   1st and 2nd coeffs except of 2nd degree poly
!                            are small                       
!         sflag(6) = true    OK    
!         sflag(6) = false   2nd degree poly have all coefficients small
    LOGICAL :: hzflag ! hzflag = .true.  OK!
                      ! hzflag = .false. abs(root)>10^5
    LOGICAL :: hwflag ! hwflag = .true.  OK!
                      ! hwflag = .false. abs(root)>10^5
    LOGICAL :: multfl ! multfl = .true.  OK!
                      ! multfl = .false. 0 has multiplicity > 4

! loop index                                                        
    INTEGER :: j 
! ==================================================================
                                                                        
!   flags initialization                                              
    weier = .true.
!    weier = .false. ! dummy
    morse = .true.
    sflag(1) = .true. 
    sflag(2) = .true. 
    sflag(3) = .true.
    sflag(4) = .true.
    sflag(5) = .true.
    sflag(6) = .true.
    warnflag(1) = .true.
    warnflag(2) = .true.
    warnflag(3) = .true.
    hzflag = .true.
    hwflag = .true.
    multfl = .true.
                                                                        
    CALL compute_statpts(ekpl,elkep,u,upl,D2,nstat,nummin,nummax,answer, &
    & morse,weier,warnflag,sflag)                                    
                                                                        
!     ******************************************************************
!     ******************** WRITING OUTPUT ON SCREEN ********************
!     ******************************************************************
!      WRITE(*,*)'######################################################
!     *########'                                                        
!      WRITE(*,*)'######### S T A T I O N A R Y   P O I N T S ##########
!     *########'                                                        
!      WRITE(*,*)'######################################################
!     *########'                                                        
!      WRITE(*,*)'       u              upl             DIST          TY
!     *E'                                                               
!      WRITE(*,*)'======================================================
!     *========'                                                        

!   loop on number of stationary points                               
    DO j = 1,nstat 
!         write(*,*)                                                    
       IF (answer(j).eq.1) THEN 
          singtype='MAXIMUM' 
       ELSEIF (answer(j).eq.-1) THEN 
          singtype='MINIMUM' 
       ELSEIF (answer(j).eq.0) THEN 
          singtype='SADDLE' 
       ELSEIF (answer(j).eq.2) THEN 
          singtype='ERROR' 
       ENDIF
                                                                        
!   writing on screen                                                 
!       WRITE(*,105)u(j),upl(j),dsqrt(D2(j)),singtype                 
!   writing on file iunit.statpts                                     
       WRITE(iunit,105)u(j),upl(j),dsqrt(D2(j)),singtype 
105    FORMAT(2x,f13.8,3x,f13.8,3x,f13.8,5x,a7) 
                                                                        
!   ====================================================              
!   writing on stat_data file for matlab scripts                      
       WRITE(iun5,106)nstat,nummin,answer(j),u(j),upl(j) 
106    FORMAT(1x,i3,2x,i3,2x,i3,2x,f11.6,2x,f11.6) 
!   ====================================================              
                                                                        
!   END MAIN DO LOOP                                                  
    ENDDO
                                                                        
!   ******************************************************************
!       *******  W R I T I N G  on  C H E C K   F I L E S  *******    
!   ******************************************************************

    IF(verb_moid.ge.20) THEN                                                                        
!   ============== SOLVING SYSTEM WARNINGS ==============             
       IF(.not.sflag(1)) THEN 
!          WRITE(iun1,*)'++++++++++++++++++++++++++++++++++++++++++++' 
          WRITE(iun1,*)'writestat3: THERE ARE TWO EQUALLY GOOD SOLUTIONS' 
!          WRITE(iun1,*)'NEGATIVE CHECK FOR ',name(1:lnam) 
!          WRITE(iun1,*)'++++++++++++++++++++++++++++++++++++++++++++' 
       ENDIF
       IF(.not.sflag(2)) THEN 
!       WRITE(iun1,*)'++++++++++++++++++++++++++++++++++++++++++++' 
!       WRITE(iun1,*)'        SOLVING SYSTEM WARNING!        ' 
         WRITE(iun1,*)'writestat3: NEITHER OF THE EVALUATIONS IS CLOSE TO ZERO' 
!       WRITE(iun1,*)' SELECTING THE ONE WITH THE SMALLEST VALUE ' 
!       WRITE(iun1,*)'NEGATIVE CHECK FOR ',name(1:lnam) 
!       WRITE(iun1,*)'++++++++++++++++++++++++++++++++++++++++++++' 
       ENDIF
                                                                        
!   ============= ADDITIONAL WARNING FLAGS ==============             
                                                                        
       IF(.not.warnflag(1)) THEN 
!          WRITE(iun2,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' 
          WRITE(iun2,*)'writestat3: LEADING COEFFICIENT = 0 !' 
!          WRITE(iun2,*)'NEGATIVE CHECK FOR ',name(1:lnam) 
!          WRITE(iun2,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' 
       ENDIF
       IF(.not.warnflag(2)) THEN 
!          WRITE(iun2,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!' 
!          WRITE(iun2,*)'!!!!!!!! WARNING !!!!!!!!!!' 
          WRITE(iun2,*)'writestat3: HIGHER DEGREE TERMS NOT SMALL !' 
!          WRITE(iun2,*)'NEGATIVE CHECK FOR ',name(1:lnam) 
!          WRITE(iun2,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!' 
       ENDIF
       IF(.not.warnflag(3)) THEN 
!          WRITE(iun2,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' 
!          WRITE(iun2,*)'FFT WARNING: DIFFERENCE IN THE EVALUATIONS' 
!          WRITE(iun2,*)'OF THE RESULTANT IS NOT SMALL!!!:' 
          WRITE(iun2,*)'writestat3: DIFFERENCE IN THE RESULTANT &
               & EVALUATIONS NOT SMALL' 
!          WRITE(iun2,*)'NEGATIVE CHECK FOR ',name(1:lnam) 
!          WRITE(iun2,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' 
       ENDIF
    
!   ============= CHECK WITH MORSE THEORY ===============             
       IF (.not.morse) THEN 
!   writing on screen                                                 
          WRITE(*,*)'writestat3: ERROR! MORSE THEORY VIOLATION!' 
!          WRITE(*,*)'FAILED COMPUTATION FOR ',name(1:lnam) 
!   writing on CHECK.morse_weier file                                 
          WRITE(iun3,*)'writestat3: ERROR! MORSE THEORY VIOLATION!' 
!          WRITE(iun3,*)'FAILED COMPUTATION FOR ',name(1:lnam) 
          WRITE(iun3,*)'NUMSTAT=',nstat 
          WRITE(iun3,*)'NUMMIN=',nummin 
          WRITE(iun3,*)'NUMMAX=',nummax 
       ENDIF
!   =====================================================             
                                                                        
!   ============ CHECK WITH WEIERSTRASS THEOREM =========             
       IF(.not.weier)THEN 
!   writing on screen                                                 
          WRITE(*,*)'writestat3: ERROR! WEIERSTRASS THEOREM VIOLATION!' 
!          WRITE(*,*)'FAILED COMPUTATION FOR ',name(1:lnam) 
!   writing on CHECK.morse_weier file                                 
          WRITE(iun3,*)'writestat3: ERROR! WEIERSTRASS THEOREM VIOLATION!' 
!          WRITE(iun3,*)'FAILED COMPUTATION FOR ',name(1:lnam) 
       ENDIF
!   ===============================================================   
                                                                        
!   ===============================================================   
!   IF PREVIOUS CHECK ARE SUCCESSFULL, THEN CHECK FOR HIGH NUMBER     
!   OF MINIMA (>2) AND STATIONARY POINTS (>=8)                        
!   ===============================================================   
       IF((weier).and.(morse))THEN 
!   write on file CHECK.hi-st                                         
          IF (nummin.gt.2) THEN 
             write(iun4,*)'THERE ARE',nummin,' MINIMA for ',name(1:lnam) 
          ENDIF
          IF ((nummin+nummax)*2.ge.8) THEN 
             write(iun4,*)'THERE ARE',(nummin+nummax)*2,' STATPTS for ', &
                  & name(1:lnam)
             write(iun4,*)'====================================' 
          ENDIF
       ENDIF

    ENDIF
!   ================================================================  
                                                                        
!    WRITE(*,*)                                                       
!    WRITE(*,*)'======================================================
!   *========'                                                        
!    WRITE(*,*)'######################################################
!   *########'                                                        

  END SUBROUTINE writestat3

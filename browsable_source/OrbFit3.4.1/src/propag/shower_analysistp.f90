!===================FILE shower_analysis.f=======================          
! Created by Giacomo Tommei, 25 November 2002   
! TP version 3.2, Milani and Tommei, December 2004                         
!===============================================================
!                                                                       
! SUBROUTINES:                                                          
!              showret3
!                    CONTAINS: 
!                             sort_time2                   
!                             sort_index2 
!================================================================              
!   
! =============================================================== 
! SHOWRET                                                               
! decomposition of list of close approaches (to the same planet)        
! into showers and returns                                              
! =============================================================== 
                       
SUBROUTINE showret3tp(iunlog,no,vas_trace,dt,tgap,isho,nsho,iret,nret)
  USE tp_trace  
  IMPLICIT NONE 
! ===INPUT=========================================
  INTEGER, INTENT(IN) :: iunlog           ! log output unit
  INTEGER, INTENT(IN) :: no               ! total number of 
                                              ! close approaches
  TYPE(tp_point), DIMENSION(no), INTENT(INOUT) :: vas_trace
                                              ! close approach record
  DOUBLE PRECISION, INTENT(IN) :: dt      ! max time span 
                                              ! defining a shower
  DOUBLE PRECISION, INTENT(IN) :: tgap    ! desired time gap 
!==================================================================
!=================== OUTPUT========================================
  INTEGER, PARAMETER :: nshx=10000    ! maximum number of showers
  INTEGER, PARAMETER :: nretx=50000   ! maximum number of returns
  INTEGER :: nshoret                       
  INTEGER, INTENT(OUT) :: nsho        ! showers number
  INTEGER, INTENT(OUT) :: isho(nshx)  ! showers index
  INTEGER, INTENT(OUT) :: iret(nretx) ! returns (trails) index
  INTEGER, INTENT(OUT) :: nret        ! returns (trails) number
!==================================================================
!==================================================================
  INTEGER, PARAMETER :: nox=200000              ! max close app number, total
  DOUBLE PRECISION :: tclo(nox)      ! close approach times 
  INTEGER :: lsho                   ! length of shower
  DOUBLE PRECISION :: ts1,ts2       ! time span of showers
  INTEGER :: j,js                   ! loop indexes 
!==================================================================                        
  tclo(1:no)=vas_trace(1:no)%tcla 
! first step: sort by time, which is column 1                           
  CALL sort_time2(tclo,vas_trace,no) 
! second step: cut by time span                                         
  nsho=1 
  isho(1)=1 
  DO j=2,no 
     IF(tclo(j).gt.tclo(isho(nsho))+dt.or.tclo(j).gt.tclo(j-1)+tgap)THEN
! new shower                                                            
        IF(nsho+1.gt.nshx)THEN 
           WRITE(*,*)' showret: increase nshx, was ',nshx 
           STOP ' **** SHOWRET FAILED ' 
        ENDIF
        nsho=nsho+1 
        isho(nsho)=j 
! check for time gap                                                    
        IF(tclo(j).gt.tclo(j-1)+tgap)THEN 
           WRITE(iunlog,*)'showret: shower gap ',                   &
     &              tclo(j)-tclo(j-1),j,nsho                            
           WRITE(iunlog,*)vas_trace(j)%rindex,tclo(j),vas_trace(j-1)%rindex,tclo(j-1) 
        ENDIF
     ENDIF
  ENDDO
  isho(nsho+1)=no+1 
! loop on showers                                                       
  nret=0 
  iret(1)=1 
  DO 1 js=1,nsho 
! record time span of shower                                            
     ts1=tclo(isho(js)) 
     ts2=tclo(isho(js+1)-1) 
     nshoret=1 
     lsho=isho(js+1)-isho(js) 
! third step: sort by multiple sol. index inside the shower             
!        WRITE(*,*)js,lsho                                              
     CALL sort_index2(tclo,vas_trace,no,isho(js),lsho) 
! first return begins at beginning of shower                            
     nret=nret+1 
     iret(nret)=isho(js) 
! cut by non consecutive indexes                                        
     DO j=isho(js)+1,isho(js+1)-1 
        IF(vas_trace(j)%rindex.eq.vas_trace(j-1)%rindex)THEN 
           WRITE(iunlog,*)'duplicate point, shower ', js,' return ',nret 
           WRITE(iunlog,*) vas_trace(j)%rindex,tclo(j),tclo(j-1) 
        ELSEIF(vas_trace(j)%rindex.gt.vas_trace(j-1)%rindex+1)THEN 
! new return                                                            
           nret=nret+1 
           nshoret=nshoret+1 
           iret(nret)=j 
        ENDIF
     ENDDO
     WRITE(iunlog,199)js,nshoret,ts1,ts2 
199  FORMAT(' Shower ',i4,' split in ',i4,' returns',               &
     &   ' time from ',f11.5 , ' to ', f11.5)                           
! end loop on showers                                                   
1 ENDDO
  iret(nret+1)=no+1 
     
CONTAINS
                                           
! ===================================================================== 
! SUBROUTINE sort_time2                                                             
! sorts a list of close approaches by time
! =====================================================================
  SUBROUTINE sort_time2(t,a,no) 
!================INPUT==================================================        
    INTEGER, INTENT(IN) :: no                        ! actual no of records  
    DOUBLE PRECISION, INTENT(INOUT) :: t(no)         ! sort key
    TYPE(tp_point), DIMENSION(no),INTENT(INOUT) :: a !records        
!=============END INTERFACE============================================
    DOUBLE PRECISION  :: ts(nox)         ! workspace: sorted time
    TYPE(tp_point), DIMENSION(nox) :: as ! workspace: sorted records
!    INTEGER :: is(SIZE(t)) 
    INTEGER :: ipo(nox)                  ! pointers 
    LOGICAL :: sorted                    ! done sorting?
    INTEGER :: j,k                       ! loop indexes 
!=======================================================================
! sort by time, get pointers                                            
    CALL heapsort (t,no,ipo)                                
! reorder                                                               
    DO j=1,no 
       ts(j)=t(ipo(j)) 
       as(j)=a(ipo(j))
    ENDDO
! output in input arrays                                                
    DO j=1,no 
       t(j)=ts(j) 
       a(j)=as(j) 
    ENDDO
 
  END SUBROUTINE sort_time2
                                           

! =====================================================================
! SUBROUTINE sort_index2
! sorts a list of close approaches by index of multiple solution        
! but only starting from record number isho
! =====================================================================
  SUBROUTINE sort_index2(t,a,no,isho,lsho) 
!===============INPUT============================             
    INTEGER, INTENT(IN) :: no                ! total no of records
    INTEGER, INTENT(IN) :: isho              ! starting row
    INTEGER, INTENT(IN) :: lsho              ! number of rows
    DOUBLE PRECISION,INTENT(INOUT) :: t(no)  ! time
    TYPE(tp_point), INTENT(INOUT) :: a(no)   !records
    DOUBLE PRECISION ts(nox)                 ! workspace: sorted time
    TYPE(tp_point), DIMENSION(nox) :: as     ! workspace: sorted records
!===============END INTERFACE============================  
    INTEGER i(nox)                   ! sort key 
    INTEGER ipo(nox)                ! pointers 
    LOGICAL sorted                  ! done sorting?
    INTEGER j,k,jj                  ! loop indexes
!=======================================================
! sort by index, get pointers 
    DO j=1,no
       i(j)=a(j)%rindex
    END DO
    CALL heapsorti(i(isho),lsho,ipo)               
! reorder                                                               
    DO j=1,lsho 
       jj=j+isho-1 
       ts(jj)=t(ipo(j)+isho-1) 
       as(jj)=a(ipo(j)+isho-1) 
    ENDDO
! output in input arrays                                                
    DO j=1,lsho 
       jj=j+isho-1 
       t(jj)=ts(jj)  
       a(jj)=as(jj) 
    ENDDO
  END SUBROUTINE sort_index2

END SUBROUTINE showret3tp
!==============================================================


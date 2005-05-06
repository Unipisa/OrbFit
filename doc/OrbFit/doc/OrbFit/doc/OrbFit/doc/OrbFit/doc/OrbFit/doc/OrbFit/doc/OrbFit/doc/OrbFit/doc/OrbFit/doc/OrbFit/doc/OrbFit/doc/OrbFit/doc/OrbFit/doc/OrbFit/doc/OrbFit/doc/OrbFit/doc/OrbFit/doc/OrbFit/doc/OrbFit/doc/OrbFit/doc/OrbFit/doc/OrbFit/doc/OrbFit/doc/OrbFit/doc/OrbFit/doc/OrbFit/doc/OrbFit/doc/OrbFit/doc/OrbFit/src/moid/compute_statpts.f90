! ************************************************************      
! ** C O M P U T I N G   S T A T I O N A R Y    P O I N T S **      
! ************************************************************      
! ********** written by GIOVANNI F. GRONCHI (2001) ***********      
! ******* Department of Mathematics, UNIVERSITY of PISA ******      
! last modified August 2003 (GFG)
! NOTE! to be used by writestat3.f90 (neodys)
! ============================================================      
  SUBROUTINE compute_statpts(ekpl,elkep,u,upl,D2,nstat,nummin,nummax,answer, &
       & morse,weier,warnflag,sflag)
    USE fund_const                                    
     IMPLICIT NONE 
! elements of the Earth and of the asteroid                         
    DOUBLE PRECISION,INTENT(IN) :: ekpl(6),elkep(6) 
! eccentric anomalies                                               
    DOUBLE PRECISION,INTENT(OUT) :: u(20),upl(20) 
! SQUARED DISTANCE function                                         
    DOUBLE PRECISION,INTENT(OUT) :: D2(20) 
! number of stationary points (maximum = 20)                        
    INTEGER,INTENT(OUT) :: nstat 
! number of relative minima/maxima found                            
    INTEGER,INTENT(OUT) :: nummin,nummax 
! type of singular point : answer(j) =  1    MAXIMUM                
!                          answer(j) = -1    MINIMUM                
!                          answer(j) =  0    SADDLE                 
!                          answer(j) =  2    CANNOT DECIDE          
    INTEGER,INTENT(OUT) :: answer(20) 
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
                      ! (large root of the resultant poly)
    LOGICAL :: hwflag ! hwflag = .true.  OK!
                      ! hwflag = .false. abs(root)>10^5
                      ! (large root of 2nd deg poly)
    LOGICAL :: multfl ! multfl = .true.  OK!
                      ! multfl = .false. 0 has multiplicity > 4
! --------------------------------- end interface ------------------
! auxiliary variables                                               
    DOUBLE PRECISION :: uu,uupl,chu,chupl 
    DOUBLE PRECISION :: DD2 
! mutual reference system                                           
    DOUBLE PRECISION :: mutI,mutom,mutompl 
! mutual orbital elements                                           
    DOUBLE PRECISION :: a,e,i,om,apl,epl,ompl 
    DOUBLE PRECISION :: beta,betapl 
    COMMON/bbpl/beta,betapl  
! ordered 3-uples                                                   
    DOUBLE PRECISION :: ordu(20),ordupl(20),ordD2(20) 
    INTEGER :: ordans(20) 
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
                                                                        
!     mutual reference system                                           
         CALL mutualrefcha(ekpl,elkep,mutI,mutom,mutompl)

!     rename variables in mutual reference frame                        
      a    = elkep(1) 
      e    = elkep(2) 
      apl  = ekpl(1) 
      epl  = ekpl(2) 
      i    = mutI 
      om   = mutom 
      ompl = mutompl 
                                                                        
      beta=dsqrt(1-e**2) 
      betapl=dsqrt(1-epl**2) 
                                                                        
!     ==================================================================
!     COMPUTING FUNCTIONS of the CONSTANT ORBITAL ELEMENTS              
!     ==================================================================
      CALL aical(a,e,i,om,apl,epl,ompl) 
                                                                        
!     ==================================================================
!     COMPUTING STATIONARY POINTS FOR D2 AND DECIDING WHICH ARE         
!     RELATIVE MAXIMA/MINIMA BETWEEN THEM                               
!     ==================================================================
      CALL comp_heart_rot(u,upl,nstat,nummin,nummax, &
           & answer,warnflag,sflag,morse,weier,hzflag,hwflag,multfl)
!      CALL comp_heart(u,upl,nstat,nummin,nummax,answer,                 &
!     &     warnflag,sflag)                                              
                                                                        
!     loop on number of stationary points                               
      DO j = 1,nstat 
                                                                        
!     conversion in radians for D2eval                                  
         uu=u(j)*radeg 
         uupl=upl(j)*radeg 
                                                                        
!     ==================================================================
!     COMPUTE SQUARED DISTANCE in the points u,upl                      
!     ==================================================================
         CALL D2eval(uu,uupl,DD2) 
!     ==================================================================
                                                                        
         D2(j)=DD2 
                                                                        
!     angles between 0 and 360 degrees                                  
         CALL choosedeg(u(j),chu) 
         CALL choosedeg(upl(j),chupl) 
                                                                        
         u(j) = chu 
         upl(j) = chupl 
                                                                        
!     END MAIN DO LOOP                                                  
      ENDDO 
                                                                        
!     =====================================================             
!     sorting 3-uples (u,upl,D2) according to value of D2               
      CALL sortd2(nstat,u,upl,D2,answer) 
!     =====================================================             
                                                                        
!     ******************************************************************
!           ******************** C H E C K S ********************       
!     ******************************************************************
!     ============= CHECK WITH MORSE THEORY ===============             
      IF (2*(nummax+nummin).ne.nstat) THEN 
         morse = .false. 
      ENDIF 
!      WRITE(*,*)'MORSE=',morse                                         
!     =====================================================             
!     ============ CHECK WITH WEIERSTRASS THEOREM =========             
      IF((nummax.lt.1).or.(nummin.lt.1))THEN 
         weier = .false. 
      ENDIF 
!      WRITE(*,*)'WEIER=',morse                                         
!     ===============================================================   
                                                                        
    END SUBROUTINE compute_statpts

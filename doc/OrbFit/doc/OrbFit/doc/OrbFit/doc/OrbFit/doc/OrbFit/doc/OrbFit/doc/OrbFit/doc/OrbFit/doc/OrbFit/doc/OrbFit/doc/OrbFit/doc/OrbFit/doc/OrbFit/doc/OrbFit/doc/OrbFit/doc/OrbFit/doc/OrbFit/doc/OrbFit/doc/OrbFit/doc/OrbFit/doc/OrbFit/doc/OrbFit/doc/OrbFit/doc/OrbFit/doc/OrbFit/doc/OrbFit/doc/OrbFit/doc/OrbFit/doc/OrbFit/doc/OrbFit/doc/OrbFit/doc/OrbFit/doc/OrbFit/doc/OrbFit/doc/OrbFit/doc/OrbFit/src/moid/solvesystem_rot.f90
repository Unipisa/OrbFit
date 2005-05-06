! ******************************************************************
! GIVEN A VALUE OF t, THIS ROUTINE COMPUTES THE CORRESPONDING VALUE 
! OF s THAT SATISFIES THE SYSTEM
! ******************************************************************
! ************ written by GIOVANNI F. GRONCHI (2001) ***************
! ********** Department of Mathematics, UNIVERSITY of PISA *********
! last modified June 2004 (GFG)
! ==================================================================
  SUBROUTINE solvesystem_rot(util,upltil,nroots,zzero,wzero,sflag,hwflag) 
    USE output_control
    IMPLICIT NONE 
    DOUBLE PRECISION,INTENT(IN) :: util,upltil
    INTEGER,INTENT(INOUT) :: nroots 
    DOUBLE PRECISION, INTENT(INOUT), DIMENSION(nroots) :: zzero,wzero 
    LOGICAL,INTENT(INOUT) :: sflag(6) ! solving system messages:
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
    LOGICAL,INTENT(INOUT) :: hwflag ! hwflag = .true.  OK!
                                    ! hwflag = .false. abs(root)>10^5
!   --------- end interface ------------------------------------------
    DOUBLE PRECISION :: sutil,cutil,supltil,cupltil
    DOUBLE PRECISION :: sx,cx,sy,cy
    DOUBLE PRECISION :: stmp(2) 
    DOUBLE PRECISION :: su,cu,u,supl(2),cupl(2),upl(2) 
    DOUBLE PRECISION :: suplchk(20),cuplchk(20) 
    DOUBLE PRECISION :: suplchktmp(20,20),cuplchktmp(20,20) 
    DOUBLE PRECISION :: suchk(20),cuchk(20) 
    DOUBLE PRECISION :: uu(20),v 
    DOUBLE PRECISION :: deg2w,deg1w,deg0w
    DOUBLE PRECISION :: deg2wtmp,deg0wtmp 
!                                                                       
    DOUBLE PRECISION :: A1,A2,A3,A4,A5,A6,A7,A8,A9,A10 
    DOUBLE PRECISION :: A11,A12,A13,A14,A15 
    COMMON/Aj1to15/ A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15 
!                                                                       
    DOUBLE PRECISION :: z,w,evalst(2) 
!   choice flag ( choice = 0 ------> init value                       
!               ( choice = 1 ------> stmp(1)                          
!               ( choice = 2 ------> stmp(2)                          
    INTEGER :: choice 
!   smallness parameter                                               
    DOUBLE PRECISION,PARAMETER :: eps=1.d-5  
!   auxiliary
    DOUBLE PRECISION, DIMENSION(20) :: zzerotmp(20),wzerotmp(20) 
    INTEGER :: discheck(20),nrealroots
!   loop indexes                                                      
    INTEGER :: i,j,k 
!   ==================================================================
                                                                        
    choice = 0 
    
    sutil = sin(util)
    cutil = cos(util)
    supltil = sin(upltil)
    cupltil = cos(upltil)
    
    DO i = 1,nroots 
       
       z=zzero(i) 
                                                                        
!   ==============================================                    
!   deg2w(z)*(w**2) + deg1w(z)*(w) + deg0w(z) = 0                     
       deg2wtmp = -A8*cutil*cupltil + 4*A1*z**3 + 4*A3*z + &
            & A11*cutil - A12*sutil - 8*A3*cutil**2*z - &
            & 2*A12*z**3*cutil - A9*sutil*supltil*z**4 + &
            & 2*A10*z**3*cutil*cupltil - &
            & A10*sutil*cupltil*z**4 + &
            & 2*A8*sutil*z**3*cupltil + &
            & 2*A9*z**3*cutil*supltil + &
            & 2*A7*sutil*z**3*supltil + &
            & A7*cutil*supltil*z**4 + &
            & A8*cutil*cupltil*z**4 + &
            & 12*A3*cutil*sutil*z**2 - &
            & 2*A3*cutil*sutil*z**4 - &
            & 12*A1*cutil*sutil*z**2 + &
            & 2*A1*cutil*sutil*z**4

       deg2w = deg2wtmp + 2*A9*z*cutil*supltil + &
            & 2*A10*z*cutil*cupltil + 2*A8*sutil*z*cupltil - &
            & 2*A3*cutil*sutil + 2*A1*cutil*sutil + &
            & A9*sutil*supltil - A7*cutil*supltil + &
            & 8*A1*cutil**2*z - 2*A12*z*cutil - &
            & 2*A11*sutil*z - 8*A1*cutil**2*z**3 + &
            & A12*sutil*z**4 + 8*A3*cutil**2*z**3 - &
            & A11*cutil*z**4 - 2*A11*sutil*z**3 + &
            & A10*sutil*cupltil - 4*A1*z- 4*A3*z**3 + &
            & 2*A7*sutil*z*supltil
       
       
       deg1w = -2*(1 + z**2)*(z**2*A7*cutil*cupltil + &
            & z**2*A10*sutil*supltil - z**2*A8*cutil*supltil-&
            & z**2*A9*sutil*cupltil + 2*A7*sutil*z*cupltil - &
            & 2*A8*sutil*z*supltil + 2*A9*z*cutil*cupltil - &
            & 2*A10*z*cutil*supltil + A8*cutil*supltil + &
            & A9*sutil*cupltil - A10*sutil*supltil - &
            & A7*cutil*cupltil)
       
       deg0wtmp = -2*A3*cutil*sutil - 2*A12*z**3*cutil - &
            & A7*cutil*supltil*z**4 - 2*A7*sutil*z**3*supltil-&
            & 2*A11*sutil*z**3 - 2*A9*z**3*cutil*supltil - &
            & 2*A8*sutil*z**3*cupltil + &
            & A10*sutil*cupltil*z**4 - &
            & 2*A10*z**3*cutil*cupltil - A11*cutil*z**4 + &
            & 8*A3*cutil**2*z**3 + A9*sutil*supltil*z**4 + &
            & A12*sutil*z**4 - 8*A1*cutil**2*z**3 - &
            & 2*A11*sutil*z - 2*A12*z*cutil + &
            & A8*cutil*cupltil - 8*A3*cutil**2*z + &
            & 8*A1*cutil**2*z
       deg0w = deg0wtmp - 2*A7*sutil*z*supltil + &
            & A7*cutil*supltil - A9*sutil*supltil - &
            & A10*sutil*cupltil + 2*A1*cutil*sutil -&
            & A8*cutil*cupltil*z**4 - 2*A8*sutil*z*cupltil - &
            & 2*A10*z*cutil*cupltil - 2*A9*z*cutil*supltil + &
            & 2*A1*cutil*sutil*z**4 - &
            & 12*A1*cutil*sutil*z**2 - &
            & 2*A3*cutil*sutil*z**4 + &
            & 12*A3*cutil*sutil*z**2 - A12*sutil+4*A1*z**3 - &
            & 4*A1*z - 4*A3*z**3 + 4*A3*z + A11*cutil
                                                                        

!       write(*,*)'deg2w,deg1w,deg0w',deg2w,deg1w,deg0w

!     CHECK SMALLNESS of the COEFFICIENTS
       IF(abs(deg2w).lt.1.d-15)THEN 
          sflag(4) = .false.
          if(verb_moid.ge.20) then
             WRITE(*,*)'solvesystem_rot: leading coefficient of 2nd degree &
                  & polynomial has magnitude less than 1.d-15',deg2w 
          endif
          IF(abs(deg1w).lt.1.d-15)THEN 
             sflag(5) = .false.
             if(verb_moid.ge.20) then
                WRITE(*,*)'solvesystem_rot: also 2nd coefficient of 2nd &
                     & degree polynomial has magnitude less than 1.d-15',deg1w 
             endif
             IF(abs(deg0w).lt.1.d-15)THEN 
                sflag(6) = .false.
                if(verb_moid.ge.20) then
                   WRITE(*,*)'solvesystem_rot: almost completely degenerate &
                        & 2nd degree polynomial',deg0w,deg1w,deg2w
                endif
             ENDIF
          ENDIF
          GOTO 124 
       ENDIF

!     NORMALIZING COEFFICIENTS (the order is IMPORTANT!)                
       deg0w = deg0w/deg2w 
       deg1w = deg1w/deg2w 
       deg2w = deg2w/deg2w 

!     check the discriminant
       IF(deg1w**2 - 4.d0*deg0w*deg2w.lt.0.d0) THEN
          discheck(i) = 1
          sflag(3) = .false.
          write(*,*)'solvesystem_rot: negative discriminant',&
               & deg1w**2 - 4.d0*deg0w*deg2w
          GOTO 133
       ELSE
          discheck(i) = 0
       ENDIF

!     solutions                                                         
       stmp(1) = (-deg1w + dsqrt(deg1w**2 - 4.d0*deg0w*deg2w))/       &
            &        (2.d0*deg2w)                                              
       stmp(2) = (-deg1w - dsqrt(deg1w**2 - 4.d0*deg0w*deg2w))/       &
            &        (2.d0*deg2w)                                              
                                                                        
       IF((abs(stmp(1)).gt.1.d5).or.(abs(stmp(2)).gt.1.d5)) THEN
          hwflag = .false.
          if(verb_moid.ge.20) then
             write(*,*)'solvesystem_rot: large roots of 2nd deg poly'
          endif
          GOTO 124
       ENDIF

!     ======================================================            
!     deg4w(z)*(w**4) + deg3w(z)*(w**3) + deg2w(z)*(w**2) +             
!     deg1w(z)*(w) + deg0w(z) = 0                                       
                                                                        
!     conversion into eccentric anomalies                               
       sx = 2.d0*zzero(i)/(1.d0+(zzero(i))**2) 
       cx = (1.d0-(zzero(i))**2)/(1.d0+(zzero(i))**2) 
       su = sx*cutil + cx*sutil
       cu = cx*cutil - sx*sutil
       u = datan2(su,cu) 
       uu(i) = u 
!         suchk(i) = su 
!         cuchk(i) = cu 
       suchk(i) = sx 
       cuchk(i) = cx 
!                                                                       
       DO j = 1,2 
          sy =  2.d0*stmp(j)/(1.d0+(stmp(j))**2) 
          cy = (1.d0-(stmp(j))**2)/(1.d0+(stmp(j))**2) 
          supl(j) = sy*cupltil + cy*supltil
          cupl(j) = cy*cupltil - sy*supltil
          upl(j) = datan2(supl(j),cupl(j)) 
          v = upl(j) 
          
          suplchktmp(i,j) =  sy 
          cuplchktmp(i,j) =  cy 
          
          evalst(j) = 2.d0*A4*sin(v)*cos(v)-2.d0*A6*cos(v)*sin(v)+    &
               & A7*sin(u)*cos(v) - A8*sin(u)*sin(v) + A9*cos(u)*cos(v) &
               & - A10*cos(u)*sin(v) + A13*cos(v) - A14*sin(v)          
       ENDDO
                                                                        

!     SELECTING the CORRESPONDING SOLUTION                              
!     ************ checking smallness of the evaluations *************  
       IF((abs(evalst(1)).lt.eps).or.(abs(evalst(2)).lt.eps))THEN 

!     choosing the one that gives the minimum value                     
          IF (abs(evalst(1)).le.abs(evalst(2))) THEN 
             choice = 1 
             wzero(i) = stmp(1) 
!               suplchk(i) =  supl(1) 
!               cuplchk(i) =  cupl(1) 

             suplchk(i) = suplchktmp(i,1) 
             cuplchk(i) = cuplchktmp(i,1)  

          ELSEIF (abs(evalst(1)).gt.abs(evalst(2))) THEN 
             choice = 2 
             wzero(i) = stmp(2) 
!               suplchk(i) =  supl(2) 
!               cuplchk(i) =  cupl(2) 

             suplchk(i) = suplchktmp(i,2) 
             cuplchk(i) = cuplchktmp(i,2)  

          ENDIF

!     if they are BOTH CLOSE TO ZERO                                    
          IF((abs(evalst(1)).lt.eps).and.(abs(evalst(2)).lt.eps))THEN
!     **************************                                         
             sflag(1) = .false.
!     **************************                                           
             if(verb_moid.ge.20) then
!     writing on screen                                                 
!               WRITE(*,*)'++++++++++++++++++++++++++++++++++++++++++++'
               WRITE(*,*)'solvesystem_rot: WARNING! TWO GOOD SOLUTIONS'        
               WRITE(*,*)'EVALST(1)=',evalst(1),'  EVALST(2)=',evalst(2)
!               WRITE(*,*)'++++++++++++++++++++++++++++++++++++++++++++'
            endif
         ENDIF
                                                                        
!     ================== DOUBLE CORRESPONDENCE CHECK ===================
!     if sflag(1) has become 0, then this check has to be done for      
!     all the following values of i                                     
                                                                        
!     -------------------------------------------------------------     
!     CASE 1 (i=1)                                                      
          IF((.not.sflag(1)).and.(i.eq.1)) THEN 
!     definitively choosing previously chosen value of s                
!     -------------------------------------------------------------     
!     CASE 2 (i>1)                                                      
          ELSEIF((.not.sflag(1)).and.(i.gt.1))THEN 
!     at first choosing previously chosen value of s                    
             DO k = 1,i-1 
                
                IF( (abs(suchk(i)-suchk(k)).lt.1.d-1).and.     &
                     &  (abs(cuchk(i)-cuchk(k)).lt.1.d-1).and.         &
                     &  (abs(suplchk(i)-suplchk(k)).lt.1.d-1).and.     &
                     &  (abs(cuplchk(i)-cuplchk(k)).lt.1.d-1) )THEN      
!     write(*,*)'OK, MODIFYING wzero to',stmp(2)                        
                                                                        
                   IF(choice.eq.1) THEN 
                      wzero(i) = stmp(2) 
                      suplchk(i) = supl(2) 
                      cuplchk(i) = cupl(2) 
                   ELSEIF(choice.eq.2) THEN 
                      wzero(i) = stmp(1) 
                      suplchk(i) = supl(1) 
                      cuplchk(i) = cupl(1) 
                   ENDIF
                   
                ELSE 
!     leaving previous value of s                                       
                   
                ENDIF
             ENDDO
          ENDIF
          
!     neither of the evaluation is close to zero                        
       ELSE 
!     ***********************                                           
          sflag(2) = .false. 
!     ***********************                                           
          if(verb_moid.ge.20) then
!          WRITE(*,*)'++++++++++++++++++++++++++++++++++++++++++++'   
!          WRITE(*,*)'        SOLVING SYSTEM WARNING!        '        
             WRITE(*,*)'solvesystem_rot: NEITHER OF THE EVALUATIONS &
                  & IS CLOSE TO ZERO'   
!          WRITE(*,*)' SELECTING THE ONE WITH THE SMALLEST VALUE '    
!          WRITE(*,*)'++++++++++++++++++++++++++++++++++++++++++++'   
             WRITE(*,*)'EVALST(1)=',evalst(1),'EVALST(2)=',evalst(2)    
!            WRITE(*,*)'THE SMALLEST IS',min(abs(evalst(1)),            
!     *           abs(evalst(2)))                                       
          endif
!     choosing the one that gives the minimum value                     
          IF (abs(evalst(1)).le.abs(evalst(2))) THEN 
             wzero(i) = stmp(1) 
          ELSEIF (abs(evalst(1)).gt.abs(evalst(2))) THEN 
             wzero(i) = stmp(2) 
          ENDIF
                                                                        
       ENDIF

133    CONTINUE
       
    ENDDO

!     selecting roots with both real components
    nrealroots = 0
    DO i=1,nroots
       IF(discheck(i).eq.0) THEN
          nrealroots=nrealroots+1
          zzerotmp(nrealroots)=zzero(i)
          wzerotmp(nrealroots)=wzero(i)
       ENDIF
    ENDDO

! renumbering zzero,wzero,nroots
    DO i = 1,nrealroots
       zzero(i) = zzerotmp(i)
       wzero(i) = wzerotmp(i)
    ENDDO
    nroots = nrealroots
                                                                        
124 CONTINUE ! to skip computations in case of degeneration of 2nd deg poly
             ! or in case of large roots
    RETURN 
  END SUBROUTINE solvesystem_rot

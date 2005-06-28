  SUBROUTINE solvesystem_ta_shift(Vpltil,vcomtil,nroots,wzero,zzero,sflag,hwflag) 
    USE output_control
    IMPLICIT NONE 
    REAL(KIND=8),INTENT(IN) :: Vpltil,vcomtil
    INTEGER,INTENT(INOUT) :: nroots 
    REAL(KIND=8),INTENT(INOUT), DIMENSION(nroots) :: wzero 
! is inout to eliminate components if the corresponding zzero is complex
    REAL(KIND=8),INTENT(OUT), DIMENSION(nroots) :: zzero 
    LOGICAL,INTENT(INOUT) :: sflag(6) ! solving system messages:
!         sflag(1) = true    OK
!         sflag(1) = false   there are two good solutions             
!         sflag(2) = true    OK                                               
!         sflag(2) = false   neither of the two evaluation is close to 0
!         sflag(3) = true    OK(the discriminant of the 2nd-deg poly is non-negative)
!         sflag(3) = false   the z-component of the solution is complex
!         sflag(4) = true    OK    
!         sflag(4) = false   leading coeff. of the 2nd degree pol. is small
!         sflag(5) = true    OK    
!         sflag(5) = false   1st and 2nd coeffs except of 2nd degree poly
!                            are small                       
!         sflag(6) = true    OK    
!         sflag(6) = false   2nd degree poly have all coefficients small
    LOGICAL,INTENT(INOUT) :: hwflag ! hwflag = .true.  OK!
                                    ! hwflag = .false. abs(root)>10^5
!   2nd degree poly smallness parameter
    REAL(KIND=8),PARAMETER :: eps_2deg=1.d-25
!   2nd degree poly bigness parameter
    REAL(KIND=8),PARAMETER :: big_root=1.d25
!   evaluation smallness parameter
    REAL(KIND=8),PARAMETER :: eps_eval=1.d-3  
!   --------- end interface ------------------------------------------
    REAL(KIND=8) :: svcomtil,cvcomtil,sVpltil,cVpltil
    REAL(KIND=8) :: sx,cx,sy,cy
    REAL(KIND=8) :: stmp(2)
    REAL(KIND=8) :: sf,cf,vcom,sfpl(2),cfpl(2),Vpl(2)
    REAL(KIND=8) :: sfplchk(N),cfplchk(N)
    REAL(KIND=8) :: sfplchktmp(N,N),cfplchktmp(N,N)
    REAL(KIND=8) :: sfchk(N),cfchk(N)
    REAL(KIND=8) :: vvcom(N),VVpl
    REAL(KIND=8) :: deg2z,deg1z,deg0z
!    REAL(KIND=8) :: deg2zold,deg1zold,deg0wold
    REAL(KIND=8) :: deg2ztmp,deg0wtmp
!                        
    REAL(KIND=8) :: z,w,evalst(2),evalstdeg2w(2)
!   choice flag ( choice = 0 ------> init value
!               ( choice = 1 ------> stmp(1)
!               ( choice = 2 ------> stmp(2)
    INTEGER :: choice
!   auxiliary
    REAL(KIND=8), DIMENSION(20) :: wzerotmp(N),zzerotmp(N) 
    LOGICAL :: discheck(N) ! check the discriminant of the 2nd-deg poly for each root
                           ! discheck(i) = .true. OK
                           ! discheck(i) = .false. negative discriminant 
    INTEGER :: nrealroots
!   loop indexes                                                      
    INTEGER :: i,j,k 
!   ==================================================================

!    verb_moid=20

! initialization
    discheck(1:N)=.true.
    sfplchktmp(1:N,1:N)=0.d0
    cfplchktmp(1:N,1:N)=0.d0
! initialization
    sfplchk(1:N)=0.d0
    cfplchk(1:N)=0.d0
    sfchk(1:N)=0.d0
    cfchk(1:N)=0.d0                        
    vvcom(1:N)=0.d0
    
    svcomtil = sin(vcomtil)
    cvcomtil = cos(vcomtil)
    sVpltil = sin(Vpltil)
    cVpltil = cos(Vpltil)
    
    DO i = 1,nroots     
       w=wzero(i) 
       write(*,*)'root number',i,' w=',w
! ==============================================                    
! deg2z(w)*z**2 + deg1z(w)*z + deg0z(w) = 0 
! (p2(w)*z**2 + p1(w)*z + p0(w) = 0)                     
       deg2z = p2(5)*w**4 + p2(4)*w**3 + p2(3)*w**2 + p2(2)*w + p2(1)
       deg1z = p1(5)*w**4 + p1(4)*w**3 + p1(3)*w**2 + p1(2)*w + p1(1)
       deg0z = p0(5)*w**4 + p0(4)*w**3 + p0(3)*w**2 + p0(2)*w + p0(1)
       write(*,*)'deg2z,deg1z,deg0z',deg2z,deg1z,deg0z

! CHECK SMALLNESS of the COEFFICIENTS
       IF(abs(deg2z).lt.eps_2deg)THEN 
          ! linear equation
          sflag(4) = .false.
          if(verb_moid.ge.20) then
             WRITE(*,*)'small leading coeff of 2nd degree poly',deg2z 
          endif
          IF(abs(deg1z).ge.eps_2deg) THEN
             zzero(i) = -deg0z/deg1z
             discheck(i) = .true. ! dummy value
             GOTO 133
          ELSEIF(abs(deg1z).lt.eps_2deg)THEN 
             sflag(5) = .false.
             if(verb_moid.ge.20) then
                WRITE(*,*)'small also the coeff of the linear term',deg1z 
             endif
             IF(abs(deg0z).lt.eps_2deg)THEN 
                sflag(6) = .false.
                if(verb_moid.ge.20) then
                   WRITE(*,*)'almost completely degenerate &
                        & 2nd degree poly',deg0z,deg1z,deg2z
                endif
             ENDIF
             ! to skip the computation of this root
             GOTO 124 !actually should compute the roots of the 4th-deg poly...
          ENDIF
       ENDIF

! NORMALIZING COEFFICIENTS (the order is IMPORTANT! deg2z must be the last)
       deg0z = deg0z/deg2z
       deg1z = deg1z/deg2z
       deg2z = deg2z/deg2z
!       write(*,*)'deg2z,deg1z,deg0z after normalization:',deg2z,deg1z,deg0z

! check the discriminant
       IF(deg1z**2 - 4.d0*deg0z*deg2z.lt.0.d0) THEN
          discheck(i) = .false.
          sflag(3) = .false.
          if (verb_moid.ge.20) then
             write(*,*)'negative discriminant',deg1z**2 - 4.d0*deg0z*deg2z
          endif
! *****************************************
          GOTO 133 ! skip this root
! *****************************************
       ELSE
          discheck(i) = .true.
          sflag(3) = .true.
       ENDIF

! solutions                                                         
       stmp(1) = (-deg1z + dsqrt(deg1z**2 - 4.d0*deg0z*deg2z))/       &
            &        (2.d0*deg2z)                                              
       stmp(2) = (-deg1z - dsqrt(deg1z**2 - 4.d0*deg0z*deg2z))/       &
            &        (2.d0*deg2z)                                              
                                                                        
       IF((abs(stmp(1)).gt.big_root).or.(abs(stmp(2)).gt.big_root)) THEN
          hwflag = .false.
          if(verb_moid.ge.20) then
             write(*,*)'solvesystem: large roots of 2nd deg poly!',abs(stmp(1:2))
          endif
          GOTO 124
       ENDIF

!     ======================================================            
!     deg4z(w)*z**4 + deg3z(w)*z**3 + deg2z(w)*z**2 + deg1z(w)*z + deg0z(w) = 0
!     (q4(w)*z**4 + q3(w)*z**3 + q2(w)*z**2 + q1(w)*z + q0(w) = 0)

!     conversion into true anomalies                               
       sx = 2.d0*wzero(i)/(1.d0+(wzero(i))**2) 
       cx = (1.d0-(wzero(i))**2)/(1.d0+(wzero(i))**2) 
       sf = sx*cvcomtil + cx*svcomtil
       cf = cx*cvcomtil - sx*svcomtil
       vcom = datan2(sf,cf) 
       vvcom(i) = vcom ! unshifted true anomaly of the comet   
       sfchk(i) = sx 
       cfchk(i) = cx 
!                                                                       
       DO j = 1,2 
          sy =  2.d0*stmp(j)/(1.d0+(stmp(j))**2) 
          cy = (1.d0-(stmp(j))**2)/(1.d0+(stmp(j))**2) 
          sfpl(j) = sy*cVpltil + cy*sVpltil
          cfpl(j) = cy*cVpltil - sy*sVpltil
          Vpl(j) = datan2(sfpl(j),cfpl(j)) 
          VVpl = Vpl(j) ! unshifted true anomaly of the planet
          sfplchktmp(i,j) =  sy 
          cfplchktmp(i,j) =  cy 

          evalstdeg2w(j) = -Ppl*(1.d0+ecom*cos(vvcom(i)))* &
               & ( -ecom*MM*cos(VVpl) - ecom*NN*sin(VVpl) + &
               & sin(vvcom(i))*KK*cos(VVpl) + sin(vvcom(i))*LL*sin(VVpl) - &
               & cos(vvcom(i))*MM*cos(VVpl) - cos(vvcom(i))*NN*sin(VVpl)) -&
               & ecom*pcom*sin(vvcom(i))*(1.d0+Epl*cos(VVpl))
!
          evalst(j) = -pcom*(1.d0+Epl*cos(VVpl))* &
               & ( Epl*LL*cos(vvcom(i)) + Epl*NN*sin(vvcom(i)) + &
               & cos(VVpl)*LL*cos(vvcom(i)) + cos(VVpl)*NN*sin(vvcom(i)) -&
               & sin(VVpl)*KK*cos(vvcom(i)) - sin(VVpl)*MM*sin(vvcom(i))) +&
               & Epl*Ppl*sin(VVpl)*(1.d0+ecom*cos(vvcom(i)))
       ENDDO

       write(*,*)'************* root(',i,')=',w
       WRITE(*,*)'EVALSTdeg2w(1)=',evalstdeg2w(1), &
            &'  EVALSTdeg2w(2)=',evalstdeg2w(2)
       WRITE(*,*)'EVALST(1)=',evalst(1),'  EVALST(2)=',evalst(2)
       write(*,*)'--------------------------------------------'

! -------------------------------------------
!     SELECTING the CORRESPONDING SOLUTION                              
       choice = 0 ! initialization 
!     ************ checking smallness of the evaluations *************  
       IF((abs(evalst(1)).lt.eps_eval).or.(abs(evalst(2)).lt.eps_eval))THEN 

!     choosing the one that gives the minimum value                     
          IF (abs(evalst(1)).le.abs(evalst(2))) THEN 
             choice = 1 
             zzero(i) = stmp(1) 
             sfplchk(i) = sfplchktmp(i,1) 
             cfplchk(i) = cfplchktmp(i,1)  

          ELSEIF (abs(evalst(1)).gt.abs(evalst(2))) THEN 
             choice = 2 
             zzero(i) = stmp(2) 
             sfplchk(i) = sfplchktmp(i,2) 
             cfplchk(i) = cfplchktmp(i,2)  

          ENDIF

!     if they are BOTH CLOSE TO ZERO                                    
          IF((abs(evalst(1)).lt.eps_eval).and.(abs(evalst(2)).lt.eps_eval))THEN
!     **************************                                         
             sflag(1) = .false.
!     **************************                                           
             if(verb_moid.ge.20) then
!     writing on screen                                                 
               WRITE(*,*)'SOLVING SYSTEM WARNING: TWO GOOD SOLUTIONS'        
               WRITE(*,*)'EVALST(1)=',evalst(1),'  EVALST(2)=',evalst(2)
            endif
         ENDIF
                                                                        
!     ================== DOUBLE CORRESPONDENCE CHECK ===================
!     if sflag(1)-.false., then do this check for all the following values of i
!     -------------------------------------------------------------     
!     CASE 1 (i=1)                
          IF((.not.sflag(1)).and.(i.eq.1)) THEN 
!     definitively choosing previously chosen value of z
!     -------------------------------------------------------------     
!     CASE 2 (i>1)                    
          ELSEIF((.not.sflag(1)).and.(i.gt.1))THEN 
!     at first choosing previously chosen value of z
             DO k = 1,i-1 
                IF( (abs(sfchk(i)-sfchk(k)).lt.1.d-1).and.     &
                     &  (abs(cfchk(i)-cfchk(k)).lt.1.d-1).and.         &
                     &  (abs(sfplchk(i)-sfplchk(k)).lt.1.d-1).and.     &
                     &  (abs(cfplchk(i)-cfplchk(k)).lt.1.d-1) )THEN      
!                   write(*,*)'OK, MODIFYING zzero to',stmp(2)
                   IF(choice.eq.1) THEN 
                      zzero(i) = stmp(2) 
                      sfplchk(i) = sfpl(2) 
                      cfplchk(i) = cfpl(2) 
                   ELSEIF(choice.eq.2) THEN 
                      zzero(i) = stmp(1) 
                      sfplchk(i) = sfpl(1) 
                      cfplchk(i) = cfpl(1) 
                   ELSE
                      WRITE(*,*)'ERROR: choice=',choice
                      STOP
                   ENDIF
                ELSE 
!     select previous value of z                                              
                ENDIF
             ENDDO
          ENDIF
          
!     neither of the evaluation is close to zero                        
       ELSE 
!     ***********************                                           
          sflag(2) = .false. 
!     ***********************                                           
          if(verb_moid.ge.20) then
             WRITE(*,*)' NEITHER OF THE EVALUATIONS IS CLOSE TO ZERO'   
             WRITE(*,*)' SELECTING THE ONE WITH THE SMALLEST VALUE '    
             WRITE(*,*)'EVALST(1)=',evalst(1),'EVALST(2)=',evalst(2)    
!            WRITE(*,*)'THE SMALLEST IS',min(abs(evalst(1)),            
!     *           abs(evalst(2)))                                       
          endif
!     choosing the one that gives the minimum value                     
          IF (abs(evalst(1)).le.abs(evalst(2))) THEN 
             zzero(i) = stmp(1) 
          ELSEIF (abs(evalst(1)).gt.abs(evalst(2))) THEN 
             zzero(i) = stmp(2) 
          ENDIF
                                                                        
       ENDIF

133    CONTINUE
       
    ENDDO

!     selecting roots with both real components
    nrealroots = 0
    DO i=1,nroots
       IF(discheck(i).eq..true.) THEN
          nrealroots=nrealroots+1
          wzerotmp(nrealroots)=wzero(i)
          zzerotmp(nrealroots)=zzero(i)
       ENDIF
    ENDDO

! renumbering wzero,zzero,nroots
    DO i = 1,nrealroots
       wzero(i) = wzerotmp(i)
       zzero(i) = zzerotmp(i)
    ENDDO
    nroots = nrealroots
                                                                        
124 CONTINUE ! to skip computations in case of degeneration of 2nd deg poly
             ! or in case of large roots

  END SUBROUTINE solvesystem_ta_shift

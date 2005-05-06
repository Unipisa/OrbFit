! ******************************************************************
! ******** HEART of the COMPUTATION of the STATIONARY POINTS *******
! ******************************************************************
! *********** written by GIOVANNI F. GRONCHI (2003) ****************
! ********** Department of Mathematics, UNIVERSITY of PISA *********
! last modified June 2004 (GFG)
! ==================================================================
  SUBROUTINE comp_heart_rot(u,upl,nroots,nummin,nummax, &
       & answer,warnflag,sflag,morse,weier,hzflag,hwflag,multfl)  
    USE fund_const                                      
    USE output_control
    IMPLICIT NONE                               
    DOUBLE PRECISION,INTENT(OUT) :: u(20),upl(20) ! eccentric anomalies  
    INTEGER,INTENT(OUT) :: nroots,nummin,nummax ! number of roots, 
                                                  !minima and maxima
    INTEGER,INTENT(OUT) :: answer(20) ! type of critical point:
!                                         answer =  1      maximum
!                                         answer =  0      saddle
!                                         answer = -1      minimum     
!                                         answer = -2      cannot decide
    LOGICAL,INTENT(INOUT) :: morse ! check with Morse theory 
!                        morse = true      OK   
!                        morse = false     ERROR  
    LOGICAL,INTENT(INOUT) :: weier ! check with Weierstrass theory:
!                        weier = true      OK 
!                        weier = false     ERROR
    LOGICAL,INTENT(INOUT) :: warnflag(3) ! program warning flags:
!      warnflag(1) = true    OK
!      warnflag(1) = false   leading coefficient of resultant is very small
!      warnflag(2) = true    OK                                          
!      warnflag(2) = false   higher degree terms in resultant are not small
!      warnflag(3) = true    OK                                           
!      warnflag(3) = false   low precision for resultant coefficients     
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
    LOGICAL,INTENT(INOUT) :: hzflag ! hzflag = .true.  OK!
                                    ! hzflag = .false. abs(root)>10^5
                                    ! (large root of the resultant poly)
    LOGICAL,INTENT(INOUT) :: hwflag ! hwflag = .true.  OK!
                                    ! hwflag = .false. abs(root)>10^5 
                                    ! (large root of 2nd deg poly)
    LOGICAL,INTENT(INOUT) :: multfl ! multfl = .true.  OK!
                                    ! multfl = .false. 0 has multiplicity > 4
!   =========== end interface ========================================
    INTEGER :: nsol
!     resultant evaluations (32) and coefficients
    DOUBLE PRECISION :: polycoe(32),evalpoly(32) 
!     for a test of the coefficient found                               
    DOUBLE PRECISION :: evalres(32),testevalpoly(32) 
    DOUBLE PRECISION :: zzero(20) ! roots of the resultant
!    DOUBLE PRECISION :: radius(20)! error estimate  
    DOUBLE PRECISION :: wzero(20) ! corresponding value of s 
!     for the eccentric anomalies                                       
    DOUBLE PRECISION :: su(20),cu(20) 
    DOUBLE PRECISION :: supl(20),cupl(20) 
!     for the angle translations
    DOUBLE PRECISION :: util,sutil,cutil,upltil,supltil,cupltil
    DOUBLE PRECISION :: sx,cx,sy,cy
!     Sylvester matrix elements                                         
    DOUBLE PRECISION :: alpha(32),beta(32),gamma(32) 
    DOUBLE PRECISION :: A(32),B(32),C(32),D(32),E(32) 
    DOUBLE PRECISION :: alphatil(32),betatil(32),gammatil(32) 
    DOUBLE PRECISION :: Atil(32),Btil(32),Ctil(32),Dtil(32),Etil(32) 
!     complex variable                                                  
    COMPLEX*16 :: complalpha(32),complbeta(32),complgamma(32) 
    COMPLEX*16 :: complA(32),complB(32),complC(32) 
    COMPLEX*16 :: complD(32),complE(32) 
    DOUBLE PRECISION :: SYLV(6,6,32) ! Sylvester matrix 
!     complex variables                                                 
    COMPLEX*16 :: evalSYLV(6,6,32),evalSYLVj(6,6) 
    COMPLEX*16 :: complSYLVline(6) 
    COMPLEX*16 :: detevalSYLV,deteval(32) ! for the determinant 
!    COMPLEX*16 :: zr1(20) ! complex zeros of the polynomial
    DOUBLE PRECISION :: s,t ! polynomial system variables 
    INTEGER :: count ! counter
    INTEGER :: j,l,nr ! loop indexes  
    INTEGER :: h,k ! matrix indexes 
    INTEGER :: ncoe,leftcoe ! indexes for the powers of the resultant
    INTEGER :: N,poldeg ! number of evaluations and degree of res(t)
    INTEGER :: ans !  ans = 1: maximum 
!                     ans = 0: saddle 
!                     ans= -1: minimum 
!                     ans= -2: cannot decide
!   ==================================================================
                                                                        
!   angular shift initialization
    util = 0.d0
    upltil = 0.d0

!   to force the rotation
!    util = 0.d0
!    upltil = dpig/4.d0
                                
    poldeg = 20 ! degree of the polynomial 
    N = 32 ! number of evaluations  
    count = 0 ! counter initialization

! to skip rotation
!    count = 10

10  CONTINUE

!   initialization of coefficients of polynomial res[P1(t,s),P2(t,s),t
!   (must be set to zero also from 22-th to 32-th for testing         
!   them again with DFT algorithm)                                    
    DO l = 1,N 
       polycoe(l)=0.d0 
       evalres(l)=0.d0 
       testevalpoly(l)=0.d0 
    ENDDO

    nummin = 0 
    nummax = 0 
                                                                        
!   ******************************************************************
!   ====================== READING MATRIX DATA =======================
!   ******************************************************************
    call matrixdat_rot(util,upltil,N,alpha,beta,gamma,A,B,C,D,E) 
                                                                        
!   ******************************************************************
!   DEFINING alphatil,betatil,gammatil,Atil,Btil,Ctil,Dtil,Etil       
!   ******************************************************************
    DO j = 1,N 
       alphatil(j) = ((-1.d0)**(j+1))*alpha(j) 
       betatil(j) = ((-1.d0)**(j+1))*beta(j) 
       gammatil(j) = ((-1.d0)**(j+1))*gamma(j) 
       Atil(j) = ((-1.d0)**(j+1))*A(j) 
       Btil(j) = ((-1.d0)**(j+1))*B(j) 
       Ctil(j) = ((-1.d0)**(j+1))*C(j) 
       Dtil(j) = ((-1.d0)**(j+1))*D(j) 
       Etil(j) = ((-1.d0)**(j+1))*E(j) 
    ENDDO
                                                                        
!   ******************************************************************
!   ========== EVALUATING the COEFFICIENTS of the MATRIX ===========  
!   ******************************************************************
    CALL rvfft(alpha,32,5) 
    CALL rvfft(alphatil,32,5) 
    CALL rvfft(beta,32,5) 
    CALL rvfft(betatil,32,5) 
    CALL rvfft(gamma,32,5) 
    CALL rvfft(gammatil,32,5) 
    CALL rvfft(A,32,5) 
    CALL rvfft(Atil,32,5) 
    CALL rvfft(B,32,5) 
    CALL rvfft(Btil,32,5) 
!   ========== C = 0 ===========                                      
    CALL rvfft(C,32,5) 
    CALL rvfft(Ctil,32,5) 
!   ============================                                      
    CALL rvfft(D,32,5) 
    CALL rvfft(Dtil,32,5) 
    CALL rvfft(E,32,5) 
    CALL rvfft(Etil,32,5) 
!   =================================================                 
                                                                        
!   ******************************************************************
!   ================ DECODING DFT OUTPUT  ============================
!   ******************************************************************
    CALL decode_out(N,alpha,alphatil,complalpha) 
    CALL decode_out(N,beta,betatil,complbeta) 
    CALL decode_out(N,gamma,gammatil,complgamma) 
    CALL decode_out(N,A,Atil,complA) 
    CALL decode_out(N,B,Btil,complB) 
    CALL decode_out(N,C,Ctil,complC) 
    CALL decode_out(N,D,Dtil,complD) 
    CALL decode_out(N,E,Etil,complE) 
                                                                        
!   ******************************************************************
!   ====== BUILDING EVALUATED MATRIX evalSYLV(h,k,j)  =========       
!   ******************************************************************
    DO j = 1,N 
       CALL clevcompsylv(complalpha(j),complbeta(j),complgamma(j),    &
            &        complA(j),complB(j),complC(j),complD(j),         &
            &        complE(j),evalSYLVj)                                      
       DO h = 1,6 
          DO k = 1,6 
!     building evaluated matrixes                                       
             evalSYLV(h,k,j) = evalSYLVj(h,k) 
          ENDDO
       ENDDO
    ENDDO
                                                                        
!   ******************************************************************
!   COMPUTATION OF THE EVALUATIONS OF THE RESULTANT IN THE N POINTS:  
!   ******************************************************************
    DO j = 1,N 
                                                                        
!   building matrixes whose determinant is to compute
       DO h = 1,6 
          DO k = 1,6 
             evalSYLVj(h,k)=evalSYLV(h,k,j) 
          ENDDO
       ENDDO
                                                                        
!   for a given value of j (1<j<32) compute det(evalSYLV(h,k,j))      
       CALL cdetcomp(evalSYLVj,detevalSYLV) 
                                                                        
!   evaluations of the resultant in the N-th roots of unity
       deteval(j)=detevalSYLV 

    ENDDO
                                                                        
!   ******************************************************************
!   =========== COMPUTE COEFFICIENTS of RESULTANT res(t) ==========   
!   ******************************************************************
!   coding input according to convention
    CALL code_inp(N,deteval,evalpoly) 
!   remember the evaluations of the resultant for test 2
    DO j = 1,N 
       evalres(j)=evalpoly(j) 
    ENDDO
!   compute coefficients
    CALL irvfft(evalpoly,32,5) 

! ******************************
!   CHECK LEADING COEFFICIENT
! ******************************
    IF (abs(evalpoly(poldeg+1)).le.1d-15) THEN 
       warnflag(1) = .false. 
!       if(verb_moid.ge.20) then
          WRITE(*,*)'comp_heart_rot: WARNING! very small leading coefficient &
               &  (<1D-15)!',evalpoly(poldeg+1),'apply rotation'
!       endif
       GOTO 13
    ELSE 
!   get a monic polynomial (polycoe(21) = 1)
       DO ncoe = 1,poldeg+1 
          polycoe(ncoe) = evalpoly(ncoe)/evalpoly(poldeg+1) 
!   remember the unnormalized coefficients of the resultant for test 1   
          testevalpoly(ncoe) = evalpoly(ncoe) 
       ENDDO
    ENDIF
                                                                        
!      check
!      write(*,*)'used not normalized coefficients:'
!      write(*,101) evalpoly(1:poldeg+1)
!      write(*,*)'left coefficients:'
!      write(*,101) evalpoly(poldeg+2:N)
!101   FORMAT (32(f15.5,1x))

!   ******************************************************************
!   ++++++++++++++++++++++++++ T E S T  1 ++++++++++++++++++++++++++++
!   ******************************************************************
!   CHECK THE NORMALIZED COEFFICIENTS OF THE HIGHER                 
!   THAN 20 DEGREE TERMS                                              
    DO ncoe = poldeg+2,N 
       IF (abs(evalpoly(ncoe)/evalpoly(poldeg+1)).gt.1.d-5) THEN 
          warnflag(2)=.false. 
!            WRITE(*,*)'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
!            WRITE(*,*)'$$$  WARNING: TEST 1 FAILED!      $$$'
!            WRITE(*,*)'$$$ HIGHER DEGREE TERMS NOT SMALL $$$'
!            WRITE(*,*)'  coeff(',ncoe,')=',evalpoly(ncoe)/evalpoly(poldeg+1)
!            WRITE(*,*)'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
       ENDIF
    ENDDO
                                                                        
!   ******************************************************************
!   ++++++++++++++++++++++++++ T E S T  2 ++++++++++++++++++++++++++++
!   ******************************************************************
!   EVALUATING THE POLYNOMIAL RESULTANT OBTAINED REMOVING THE 11      
!   COEFFICIENT FROM DEGREE 32 UP TO DEGREE 22 INCLUDED               
!   ******************************************************************
    CALL rvfft(testevalpoly,32,5) 
    DO j = 1,N 
!   =========================== CHECK ============================    
       IF ( abs(testevalpoly(j)-evalres(j)).gt.1.d-5) THEN 
!   HINT: polycoe(j) are coded, so we have to test them with evalres(j)
          warnflag(3) = .false. 
!          WRITE(*,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'   
!          WRITE(*,*)'ERROR: TEST 2 FAILED!'
!          WRITE(*,*)'DIFFERENCE IN THE EVALUATIONS OF'        
!          WRITE(*,*)'THE RESULTANT IS NOT SMALL (> 1.D-5):'                 
!          WRITE(*,*)'diff=',abs(testevalpoly(j)-evalres(j))          
!          WRITE(*,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'   
       ENDIF
    ENDDO
!   ******************************************************************
                                                                        
                                                                        
!solvpoly(poldeg,coef,roots,nroots,hzflag,multfl)
!   ******************************************************************
!   FIND REAL ROOTS OF THE POLYNOMIAL res[ P1(t,s),P2(t,s),t ]        
    CALL solvpoly(20,polycoe(1:21),zzero,nsol,hzflag,multfl) 
!   ******************************************************************
    if(.not.hzflag) then
       goto 13
    endif
 
!      check               
!      write(*,*)'nroots=',nsol
!      write(*,*)'zzero:',zzero(1:nsol)
                                                        
!   ************************************************************
!   COMPUTE FOR EACH t THE CORRESPONDING VALUE of s                   
    CALL solvesystem_rot(util,upltil,nsol,zzero(1:nsol),wzero,sflag,hwflag) 
!   ************************************************************
    if(.not.hwflag) then
       goto 13
    endif

    sutil = sin(util)
    cutil = cos(util)
    supltil = sin(upltil)
    cupltil = cos(upltil)
    
    DO nr = 1,nsol 

!   conversion from (t,s) to the angles (u,upl)
       sx = 2.d0*zzero(nr)/(1.d0+(zzero(nr))**2) 
       cx = (1.d0-(zzero(nr))**2)/(1.d0+(zzero(nr))**2) 
       su(nr) = sx*cutil + cx*sutil
       cu(nr) = cx*cutil - sx*sutil

       sy =  2.d0*wzero(nr)/(1.d0+(wzero(nr))**2) 
       cy = (1.d0-(wzero(nr))**2)/(1.d0+(wzero(nr))**2) 
       supl(nr) = sy*cupltil + cy*supltil
       cupl(nr) = cy*cupltil - sy*supltil

!   critical points (eccentric anomalies)                           
       u(nr) = datan2(su(nr),cu(nr)) 
       upl(nr) = datan2(supl(nr),cupl(nr)) 
                                                                        
!   compute type of critical points
       CALL hess(u(nr),upl(nr),ans) 

!       ans = -2 !dummy
       IF(ans.eq.-2) THEN
          WRITE(*,*)'comp_heart_rot: cannot decide type of critical point, calling int_eval'
          CALL int_eval(u(nr),upl(nr),ans)
       ENDIF
!          write(*,*)'u,upl',u(nr)*degrad,upl(nr)*degrad
!          write(*,*)'ans',ans

       IF (ans.eq.-1) THEN 
          nummin = nummin + 1 
       ELSEIF (ans.eq.1)THEN 
          nummax = nummax + 1 
       ENDIF
                                                                       
       answer(nr) = ans 
    
!   conversion into degrees                                           
       u(nr) =  u(nr)*degrad 
       upl(nr) = upl(nr)*degrad 
       
    ENDDO

!   *****************************************************             
!   ============= CHECK WITH MORSE THEORY ===============             
    IF (2*(nummax+nummin).ne.nsol) THEN 
       morse = .false.
    ENDIF
!   ============ CHECK WITH WEIERSTRASS THEOREM =========             
    IF((nummax.lt.1).or.(nummin.lt.1))THEN 
       weier = .false. 
    ENDIF
!   *****************************************************             

!   **********************************************************************
!   ******************  A N G U L A R   S H I F T  ***********************
13  IF((.not.morse).or. &
         & (.not.weier).or. &
         & (.not.warnflag(1)).or.&
!         & (.not.sflag(2)).or. &
         & (.not.sflag(3)).or. &
         & (.not.sflag(4)).or. &
         & (.not.sflag(5)).or. &
         & (.not.sflag(6)).or. &
         & (.not.hzflag).or. &
         & (.not.hwflag).or. &
         & (.not.multfl)) THEN
       count = count + 1
       IF(count.ge.10) THEN
          if(verb_moid.ge.20) then
!            WRITE(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' 
             WRITE(*,*)'comp_heart_rot: COMPUTATION FAILED!'
             WRITE(*,*) count,'ANGULAR SHIFT TRIED! &
                  & STOPPING CRITICAL POINTS COMPUTATION!'
!            WRITE(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' 
          endif
          nsol = -1 
          nummin = -1 
          nummax = -1 
          GOTO 1300
       ENDIF
          if(verb_moid.ge.20) then
             WRITE(*,*)'morse,weier,warnflag(1)',morse,weier,warnflag(1)
             WRITE(*,*)'sflag(2:6)',sflag(2:6)
             WRITE(*,*)'hzflag,hwflag,multfl',hzflag,hwflag,multfl
          endif
          util = util + 6.28318530717958648d0/15.d0
          upltil = upltil + 6.28318530717958648d0/16.d0
          if(verb_moid.ge.20) then
!            WRITE(*,*)'***********************************************'
             WRITE(*,*)'comp_heart_rot: APPLY ANGULAR SHIFT! count=',count
             WRITE(*,*)'               util,upltil:',util,upltil
!            WRITE(*,*)'***********************************************'
          endif
!   restoring flags
       morse = .true.
       weier = .true.
       warnflag(1) = .true.
       warnflag(2) = .true.
       warnflag(3) = .true.
       sflag(1) = .true.
       sflag(2) = .true.
       sflag(3) = .true.
       sflag(4) = .true.
       sflag(5) = .true.
       sflag(6) = .true.
       hzflag = .true.
       multfl = .true.
       GOTO 10 
    ELSE
!         WRITE(*,*)'computation OK!'
!         WRITE(*,*)'morse,weier,warnflag(1)',morse,weier,warnflag(1)
!         WRITE(*,*)'sflag(2:6)',sflag(2:6)
!         WRITE(*,*)'hzflag',hzflag
!     continue the computation
    ENDIF

!   maximal number of angular shifts
    IF(count.ge.10) THEN
       if(verb_moid.ge.20) then
          WRITE(*,*)'comp_heart_rot: COMPUTATION FAILED!'
          WRITE(*,*) count,'ANGULAR SHIFT TRIED! &
               & STOPPING CRITICAL POINTS COMPUTATION!'
       endif
       nsol = -1 
       nummin = -1 
       nummax = -1 
    ENDIF
! **********************************************************************

1300 CONTINUE ! to skip computation in case of failure
    
    nroots = nsol
    
    RETURN 
  END SUBROUTINE comp_heart_rot

! ******************************************************************
! ** MAIN SUBROUTINE for the COMPUTATION of the STATIONARY POINTS **
! ******************************************************************
! *********** written by GIOVANNI F. GRONCHI (Sept.2004) ***********
! ==================================================================
  SUBROUTINE comp_heart_ta_shift(fpl,fcom,nroots,nummin,nummax, &
       & answer,warnflag,sflag,morse,weier,hzflag,hwflag,multfl)  
    USE fund_const                                      
    USE output_control
    IMPLICIT NONE
    INTEGER, PARAMETER :: poldeg = 16 ! polynomial degree
    INTEGER, PARAMETER :: N = 16 ! number of evaluations of the 15th deg poly
    INTEGER, PARAMETER :: expo = 4 ! 2^expo = N
    REAL(KIND=8),INTENT(OUT) :: fpl(poldeg),fcom(poldeg)! true anomalies  
    INTEGER,INTENT(OUT) :: nroots,nummin,nummax ! number of roots, 
                                                  !minima and maxima
! ------- error flags ----------------------------------------------
    INTEGER,INTENT(OUT) :: answer(poldeg) ! type of critical point:
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
! warnflag(1) = true    OK
! warnflag(1) = false   leading coefficient of \tilde{r}(t) is very small
! warnflag(2) = true    OK                                          
! warnflag(2) = false   higher degree terms in resultant are not small
! warnflag(3) = true    OK                                           
! warnflag(3) = false   low precision for resultant coefficients     
    LOGICAL,INTENT(INOUT) :: sflag(6) ! solving system messages:
! sflag(1) = true    OK
! sflag(1) = false   there are two good solutions             
! sflag(2) = true    OK                                               
! sflag(2) = false   neither of the two evaluation is close to 0
! sflag(3) = true    OK                                               
! sflag(3) = false   the s-component of the solution is complex
! sflag(4) = true    OK    
! sflag(4) = false   leading coeff. of the 2nd degree pol. is small
! sflag(5) = true    OK    
! sflag(5) = false   1st and 2nd coeffs except of 2nd degree poly
!            are small                       
! sflag(6) = true    OK    
! sflag(6) = false   2nd degree poly have all coefficients small
    LOGICAL,INTENT(INOUT) :: hzflag ! hzflag = .true.  OK!
                                    ! hzflag = .false. abs(root)>10^5
    LOGICAL,INTENT(INOUT) :: hwflag ! hwflag = .true.  OK!
                                    ! hwflag = .false. abs(root)>10^5
    LOGICAL,INTENT(INOUT) :: multfl ! multfl = .true.  OK!
                                    ! multfl = .false. 0 has multiplicity > 4
! =========== end interface ========================================
    INTEGER :: nsol
! resultant evaluations (N) and coefficients
    REAL(KIND=8) :: polycoe(N),polycoe1(N+1),evalpoly(N),evalp0(N) 
!
    COMPLEX(KIND=8) :: evpol15(N)
! for a test of the coefficient found                               
    REAL(KIND=8) :: evalrt(N),testevalpoly(N) 
    REAL(KIND=8) :: zzero(poldeg) ! roots of the resultant
!    REAL(KIND=8) :: radius(poldeg)! error estimate  
    REAL(KIND=8) :: wzero(poldeg) ! corresponding value of s 
! for the eccentric anomalies                                       
    REAL(KIND=8) :: sfpl(poldeg),cfpl(poldeg) 
    REAL(KIND=8) :: sfcom(poldeg),cfcom(poldeg) 
! for the angle translations
    REAL(KIND=8) :: fpltil,sfpltil,cfpltil,fcomtil,sfcomtil,cfcomtil
    REAL(KIND=8) :: sx,cx,sy,cy
! Sylvester matrix elements                                         
    REAL(KIND=8) :: p0(N),p1(N),p2(N)
    REAL(KIND=8) :: q0(N),q1(N),q2(N),q3(N),q4(N)
    REAL(KIND=8) :: r31(N),r32(N),r33(N),r34(N),r35(N)
    REAL(KIND=8) :: r36(N),r41(N),r45(N),r46(N)

    COMPLEX(KIND=8) :: cp0,cp1,cp2 
    COMPLEX(KIND=8) :: cq0,cq1,cq2,cq3,cq4
    COMPLEX(KIND=8) :: cr31,cr32,cr33,cr34,cr35
    COMPLEX(KIND=8) :: cr36,cr41,cr45,cr46

! constant term of det(S_1)
    COMPLEX(KIND=8) :: S10(N,N),detS10
! complex variables                                                  
    COMPLEX(KIND=8) :: complp0(N),complp1(N),complp2(N)
    COMPLEX(KIND=8) :: complq0(N),complq1(N),complq2(N),complq3(N),complq4(N)
    COMPLEX(KIND=8) :: complr31(N),complr32(N),complr33(N),complr34(N),complr35(N)
    COMPLEX(KIND=8) :: complr36(N),complr41(N),complr45(N),complr46(N)
    REAL(KIND=8) :: SYLV(6,6,N) ! Sylvester matrix 
! complex variables                                                 
    COMPLEX(KIND=8) :: evalSYLV(6,6,N),evalSYLVj(6,6) 
    COMPLEX(KIND=8) :: complSYLVline(6) 
    COMPLEX(KIND=8) :: detevalSYLV,deteval(N) ! for the determinant 
! N-th roots of unity
    REAL(KIND=8) :: puroot(N),puroottil(N)
    COMPLEX(KIND=8) :: compluroot(N)
!    COMPLEX*16 :: zr1(poldeg) ! complex zeros of the polynomial
    REAL(KIND=8) :: s,t ! polynomial system variables 
    INTEGER :: count ! counter
    INTEGER :: j,l,nr ! loop indexes  
    INTEGER :: h,k ! matrix indexes 
    INTEGER :: ncoe,leftcoe ! indexes for the powers of the resultant
    INTEGER :: ans !  ans = 1: maximum 
!                     ans = 0: saddle 
!                     ans= -1: minimum 
!                     ans= -2: cannot decide
!   ==================================================================

! angular shift initialization
! dpig=6.28318530717958648d0
    fpltil = 0.d0
    fcomtil = 0.d0
    count = 0 ! counter initialization
!    verb_moid=21
    count = 9 ! to skip the shift

10  CONTINUE

! initialization of coefficients of polynomial r(t)
    polycoe(1:N)=0.d0 
    polycoe1(1:N+1) = 0.d0
    evalrt(1:N)=0.d0 
    testevalpoly(1:N)=0.d0 

    nummin = 0 
    nummax = 0 
                                                                        
! ******************
! READ MATRIX DATA
    call matrixdat_ta_shift(fpltil,fcomtil,N,p0,p1,p2,q0,q1,q2,q3,q4, &
         & r31,r32,r33,r34,r35,r36,r41,r45,r46)
!    write(*,*)'p0,p1,p2'
!    do j = 1,5
!       write(*,*) p0(j),p1(j),p2(j)
!    enddo

! *******************************
! COMPUTE CONSTANT TERM of r(t)
    cp0 = p0(1)
    cp1 = p1(1)
    cp2 = p2(1)
    cq0 = q0(1)
    cq1 = q1(1)
    cq2 = q2(1)
    cq3 = q3(1)
    cq4 = q4(1)
    cr31 = r31(1) 
    cr32 = r32(1) 
    cr33 = r33(1) 
    cr34 = r34(1) 
    cr35 = r35(1) 
    cr36 = r36(1) 
    cr41 = r41(1) 
    cr31 = r45(1) 
    cr31 = r46(1) 

!    write(*,*)'cp0',cp0
    CALL compmodsylv16_shift(cp0,cp1,cp2,cq0,cq1,cq2,cq3,cq4,cr31,cr32,cr33, &
         & cr34,cr35,cr36,cr41,cr45,cr46,S10)

    CALL cdetcomp(S10,detS10)
!    write(*,*)'r0=',detS10

! ****************************************
! p(x)=x POLYNOMIAL (for roots of unity)
    puroot(1) = 0.d0
    puroot(2) = 1.d0
    DO j = 3,N
       puroot(j) = 0.d0
    ENDDO

! *****************************************************
! EVALUATE the COEFFICIENTS of the MATRIX \tilde{S}_1
    CALL rvfft(p0,N,expo) 
    CALL rvfft(p1,N,expo) 
    CALL rvfft(p2,N,expo) 
    CALL rvfft(q0,N,expo) 
    CALL rvfft(q1,N,expo) 
    CALL rvfft(q3,N,expo) 
    CALL rvfft(q4,N,expo) 
    CALL rvfft(r31,N,expo) 
    CALL rvfft(r32,N,expo)
    CALL rvfft(r33,N,expo)
    CALL rvfft(r34,N,expo)
    CALL rvfft(r35,N,expo) ! e' costante, evitare la valutazione!
    CALL rvfft(r36,N,expo) ! e' costante, evitare la valutazione!
    CALL rvfft(r41,N,expo) 
    CALL rvfft(r45,N,expo) ! e' costante, evitare la valutazione!
    CALL rvfft(r46,N,expo) ! e' costante, evitare la valutazione!

    complp0(1) = DCMPLX(p0(1),0.d0)
    complp1(1) = DCMPLX(p1(1),0.d0)
    complp2(1) = DCMPLX(p2(1),0.d0)
    complq0(1) = DCMPLX(q0(1),0.d0)
    complq1(1) = DCMPLX(q1(1),0.d0)
    complq2(1) = DCMPLX(q2(1),0.d0)
    complq3(1) = DCMPLX(q3(1),0.d0)
    complq4(1) = DCMPLX(q4(1),0.d0) 
    complr31(1) = DCMPLX(r31(1),0.d0) 
    complr32(1) = DCMPLX(r32(1),0.d0) 
    complr33(1) = DCMPLX(r33(1),0.d0) 
    complr34(1) = DCMPLX(r34(1),0.d0) 
    complr35(1) = DCMPLX(r35(1),0.d0) 
    complr36(1) = DCMPLX(r36(1),0.d0) 
    complr41(1) = DCMPLX(r41(1),0.d0) 
    complr45(1) = DCMPLX(r45(1),0.d0) 
    complr46(1) = DCMPLX(r46(1),0.d0) 
!
    complp0(N/2+1) = DCMPLX(p0(N/2+1),0.d0)
    complp1(N/2+1) = DCMPLX(p1(N/2+1),0.d0)
    complp2(N/2+1) = DCMPLX(p2(N/2+1),0.d0)
    complq0(N/2+1) = DCMPLX(q0(N/2+1),0.d0)
    complq1(N/2+1) = DCMPLX(q1(N/2+1),0.d0)
    complq2(N/2+1) = DCMPLX(q2(N/2+1),0.d0)
    complq3(N/2+1) = DCMPLX(q3(N/2+1),0.d0)
    complq4(N/2+1) = DCMPLX(q4(N/2+1),0.d0)
    complr31(N/2+1) = DCMPLX(r31(N/2+1),0.d0) 
    complr32(N/2+1) = DCMPLX(r32(N/2+1),0.d0) 
    complr33(N/2+1) = DCMPLX(r33(N/2+1),0.d0) 
    complr34(N/2+1) = DCMPLX(r34(N/2+1),0.d0) 
    complr35(N/2+1) = DCMPLX(r35(N/2+1),0.d0) 
    complr36(N/2+1) = DCMPLX(r36(N/2+1),0.d0) 
    complr41(N/2+1) = DCMPLX(r41(N/2+1),0.d0) 
    complr45(N/2+1) = DCMPLX(r45(N/2+1),0.d0) 
    complr46(N/2+1) = DCMPLX(r46(N/2+1),0.d0) 
!
    DO j = 1,N/2-1
       complp0(j+1) = DCMPLX(p0(j+1),p0(N-j+1))
       complp0(N/2+j+1) = DCMPLX(p0(N/2-j+1),-p0(N/2+1+j))        

       complp1(j+1) = DCMPLX(p1(j+1),p1(N-j+1))
       complp1(N/2+j+1) = DCMPLX(p1(N/2-j+1),-p1(N/2+1+j))        

       complp2(j+1) = DCMPLX(p2(j+1),p2(N-j+1))
       complp2(N/2+j+1) = DCMPLX(p2(N/2-j+1),-p2(N/2+1+j))        

       complq0(j+1) = DCMPLX(q0(j+1),q0(N-j+1))
       complq0(N/2+j+1) = DCMPLX(q0(N/2-j+1),-q0(N/2+1+j))        

       complq1(j+1) = DCMPLX(q1(j+1),q1(N-j+1))
       complq1(N/2+j+1) = DCMPLX(q1(N/2-j+1),-q1(N/2+1+j))        

       complq2(j+1) = DCMPLX(q2(j+1),q2(N-j+1))
       complq2(N/2+j+1) = DCMPLX(q2(N/2-j+1),-q2(N/2+1+j))        

       complq3(j+1) = DCMPLX(q3(j+1),q3(N-j+1))
       complq3(N/2+j+1) = DCMPLX(q3(N/2-j+1),-q3(N/2+1+j))        

       complq4(j+1) = DCMPLX(q4(j+1),q4(N-j+1))
       complq4(N/2+j+1) = DCMPLX(q4(N/2-j+1),-q4(N/2+1+j))        

       complr31(j+1) = DCMPLX(r31(j+1),r31(N-j+1))
       complr31(N/2+j+1) = DCMPLX(r31(N/2-j+1),-r31(N/2+1+j))        

       complr32(j+1) = DCMPLX(r32(j+1),r32(N-j+1))
       complr32(N/2+j+1) = DCMPLX(r32(N/2-j+1),-r32(N/2+1+j))
        
       complr33(j+1) = DCMPLX(r33(j+1),r33(N-j+1))
       complr33(N/2+j+1) = DCMPLX(r33(N/2-j+1),-r33(N/2+1+j))        

       complr34(j+1) = DCMPLX(r34(j+1),r34(N-j+1))
       complr34(N/2+j+1) = DCMPLX(r34(N/2-j+1),-r34(N/2+1+j))        

       complr35(j+1) = DCMPLX(r35(j+1),r35(N-j+1))
       complr35(N/2+j+1) = DCMPLX(r35(N/2-j+1),-r35(N/2+1+j))        

       complr36(j+1) = DCMPLX(r36(j+1),r36(N-j+1))
       complr36(N/2+j+1) = DCMPLX(r36(N/2-j+1),-r36(N/2+1+j))        

       complr41(j+1) = DCMPLX(r41(j+1),r41(N-j+1))
       complr41(N/2+j+1) = DCMPLX(r41(N/2-j+1),-r41(N/2+1+j))        

       complr45(j+1) = DCMPLX(r45(j+1),r45(N-j+1))
       complr45(N/2+j+1) = DCMPLX(r45(N/2-j+1),-r45(N/2+1+j))        

       complr46(j+1) = DCMPLX(r46(j+1),r46(N-j+1))
       complr46(N/2+j+1) = DCMPLX(r46(N/2-j+1),-r46(N/2+1+j))        
    ENDDO

! === computing N-th roots of unity ===
    CALL rvfft(puroot,N,expo) 
      compluroot(1) = DCMPLX(puroot(1),0.d0)
      compluroot(N/2+1) = DCMPLX(puroot(N/2+1),0.d0)
      DO j = 1,N/2-1
         compluroot(j+1) = DCMPLX(puroot(j+1),puroot(N-j+1))
         compluroot(N/2+j+1) = DCMPLX(puroot(N/2-j+1),-puroot(N/2+1+j))        
      ENDDO

! check
!    write(*,*)'eval p0,p1,p2'
!    do j=1,N
!       write(*,*)complp0(j),complp1(j),complp2(j)
!    enddo
! check
!    write(*,*)'radici di 1' !selected clockwise
!    do j=1,N
!       write(*,*)compluroot(j)
!    enddo
! *TEST*: compute coefficients of \tilde{r}(t)
!    CALL code_input(N,complp0,evalp0) 
!    CALL irvfft(evalp0,N,expo) 
!    write(*,*)'coeffs of p0'
!    DO j=1,N
!       write(*,*)evalp0(j)
!    ENDDO

! **********************************
! EVALUATED MATRIX evalSYLV(h,k,j)
    DO j = 1,N 
       CALL compmodsylv16_shift(complp0(j),complp1(j),complp2(j), &
            & complq0(j),complq1(j),complq2(j),complq3(j),complq4(j), &
            & complr31(j),complr32(j),complr33(j),complr34(j),complr35(j), &
            & complr36(j),complr41(j),complr45(j),complr46(j),evalSYLVj)
       DO h = 1,6 
          DO k = 1,6 
!     defining evaluated matrixes                                       
             evalSYLV(h,k,j) = evalSYLVj(h,k) 
          ENDDO
       ENDDO
    ENDDO
                                                                        
! **********************************************
! evaluation of r(t) at the N-th roots o unity
    DO j = 1,N 
! defining matrixes whose determinant is to compute
       DO h = 1,6 
          DO k = 1,6 
             evalSYLVj(h,k)=evalSYLV(h,k,j) 
          ENDDO
       ENDDO
! for a given value of j (1<j<N) compute det(evalSYLV(h,k,j))
       CALL cdetcomp(evalSYLVj,detevalSYLV)
! evaluations of the resultant in the N-th roots of unity
       deteval(j)=detevalSYLV
    ENDDO

! ********************************************************
! compute evaluation \tilde{r}: (deteval(j)-r0)/uroot(j)
    DO j = 1,N
       evpol15(j) = (deteval(j)-detS10)/compluroot(j)
    ENDDO
!       write(*,*)'15 deg poly evals:',evpol15(1:N)

! **************************************
! compute coefficients of \tilde{r}(t)
    CALL code_input(N,evpol15,evalpoly) 
!   remember the evaluations of the resultant for test 2
    DO j = 1,N 
       evalrt(j)=evalpoly(j) 
    ENDDO
    CALL irvfft(evalpoly,N,expo) 
!   check leading coefficient (poldeg=16)
    IF (abs(evalpoly(poldeg)).le.1d-15) THEN 
       warnflag(1) = .false. 
       if(verb_moid.ge.20) then
          WRITE(*,*)'WARNING! very small leading coefficient (< 1D-15)'
          WRITE(*,*)'leading coe:',evalpoly(poldeg)
       endif
       GOTO 13
    ELSE 

!   get a monic polynomial (polycoe(poldeg) = 1)     
! !!!!!(QUESTA OPERAZIONE VA FATTA DOPO con il pol di grado 16) !!!!!
       DO ncoe = 1,poldeg 
!          polycoe(ncoe) = evalpoly(ncoe)/evalpoly(poldeg) 
          polycoe(ncoe) = evalpoly(ncoe) ! lasciamo stare i coeffs
       ENDDO
    ENDIF

! ******************************              
! compute coefficients of r(t)
       polycoe1(1) = DBLE(detS10)
       DO ncoe = 1,poldeg
          polycoe1(ncoe+1) = evalpoly(ncoe)
       ENDDO
! check
!    write(*,*)'r(t) coefficients'
!    write(*,101) polycoe1(1:poldeg+1)
101 FORMAT (32(es16.6,1x))
! eventualmente normalizzare ora i coeffs

! skip test 
    GOTO 111
! ****************************************************
! TEST:  EVALUATING THE OBTAINED POLYNOMIAL OBTAINED             
    CALL rvfft(testevalpoly,N,expo) 
    DO j = 1,N 
       IF ( abs(testevalpoly(j)-evalrt(j)).gt.1.d-5) THEN 
!   HINT: polycoe(j) are coded, so we have to test them with evalrt(j)
          warnflag(3) = .false. 
          if (verb_moid.ge.20)then
             WRITE(*,*)'WARNING: TEST FAILED!'
             WRITE(*,*)'DIFFERENCE IN THE EVALUATIONS OF'        
             WRITE(*,*)'THE RESULTANT IS NOT SMALL (> 1.D-5):'                 
             WRITE(*,*)'diff=',abs(testevalpoly(j)-evalrt(j))          
          endif
       ENDIF
    ENDDO
!   ******************************************************************
111 CONTINUE
                                                                        
!   ******************************************************************
!   FIND REAL ROOTS OF THE POLYNOMIAL res[ P1(t,s),P2(t,s),t ]        
    CALL solvpoly(poldeg,polycoe1(1:poldeg+1),zzero,nsol,hzflag,multfl) 
!   ******************************************************************
    if(.not.hzflag) then
       write(*,*)'warning: hzflag=',hzflag
       goto 13
    endif
 
! check               
    write(*,*)'nroots=',nsol
    write(*,*)'zzero:',zzero(1:nsol)

!   ************************************************************
!   COMPUTE FOR EACH t THE CORRESPONDING VALUE of s                   
    CALL solvesystem_ta_shift(fpltil,fcomtil,nsol,zzero(1:nsol),wzero,sflag,hwflag) 
!   ************************************************************
    if(.not.hwflag) then
       goto 13
    endif

    sfpltil = sin(fpltil)
    cfpltil = cos(fpltil)
    sfcomtil = sin(fcomtil)
    cfcomtil = cos(fcomtil)
    
    DO nr = 1,nsol 

!   conversion from (t,s) to the true anomalies (fcom,fpl)
       sy =  2.d0*wzero(nr)/(1.d0+(wzero(nr))**2) 
       cy = (1.d0-(wzero(nr))**2)/(1.d0+(wzero(nr))**2) 
       sfpl(nr) = sy*cfpltil + cy*sfpltil
       cfpl(nr) = cy*cfpltil - sy*sfpltil

       sx = 2.d0*zzero(nr)/(1.d0+(zzero(nr))**2) 
       cx = (1.d0-(zzero(nr))**2)/(1.d0+(zzero(nr))**2) 
       sfcom(nr) = sx*cfcomtil + cx*sfcomtil
       cfcom(nr) = cx*cfcomtil - sx*sfcomtil

!   critical points (true anomalies)                           
       fpl(nr) = datan2(sfpl(nr),cfpl(nr)) 
       fcom(nr) = datan2(sfcom(nr),cfcom(nr)) 

!   compute type of critical points
       CALL hess_ta(fpl(nr),fcom(nr),ans) 

!       ans = -2 !dummy to force int_eval
       IF(ans.eq.-2) THEN
          WRITE(*,*)'cannot decide type of critical point: calling int_eval'
          CALL int_eval_ta(fpl(nr),fcom(nr),ans)
       ENDIF

       IF (ans.eq.-1) THEN 
          nummin = nummin + 1 
       ELSEIF (ans.eq.1)THEN 
          nummax = nummax + 1 
       ENDIF
                                                                       
       answer(nr) = ans 

!   conversion into degrees                                           
       fpl(nr) = fpl(nr)*degrad 
       fcom(nr) =  fcom(nr)*degrad 
       
    ENDDO

!       write(*,*)'true anomalies of the critical points'
!       do nr = 1,nsol
!          write(*,*) fcom(nr), fpl(nr), answer(nr)        
!       enddo

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
       IF(count.eq.10) THEN
          if(verb_moid.ge.20) then
!            WRITE(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' 
             WRITE(*,*)'!! COMPUTATION FAILED:',count,'ANGULAR SHIFT TRIED !!'
             WRITE(*,*)'!!!     STOPPING CRITICAL POINTS COMPUTATION      !!!'
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
       fcomtil = fcomtil + 6.28318530717958648d0/15.d0
       fpltil = fpltil + 6.28318530717958648d0/16.d0
          if(verb_moid.ge.20) then
!            WRITE(*,*)'***********************************************'
             WRITE(*,*)'APPLY ANGULAR SHIFT: count=',count
             WRITE(*,*)'               fcomtil,fpltil:',fcomtil,fpltil
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
    IF(count.eq.10) THEN
       if(verb_moid.ge.20) then
          WRITE(*,*)'!!! COMPUTATION FAILED:',count,'ANGULAR SHIFT TRIED !!!'
          WRITE(*,*)'!!!      STOPPING CRITICAL POINTS COMPUTATION       !!!'
       endif
       nsol = -1 
       nummin = -1 
       nummax = -1 
    ENDIF
! **********************************************************************

1300 CONTINUE ! to skip computiation in case of failure
    
    nroots = nsol
    
    RETURN 
  END SUBROUTINE comp_heart_ta_shift

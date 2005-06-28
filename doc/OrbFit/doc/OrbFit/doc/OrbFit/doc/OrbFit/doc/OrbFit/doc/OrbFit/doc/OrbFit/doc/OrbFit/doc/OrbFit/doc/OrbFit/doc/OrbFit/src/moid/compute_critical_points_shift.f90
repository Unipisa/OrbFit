! ******************************************************************
! ** MAIN SUBROUTINE for the COMPUTATION of the STATIONARY POINTS **
! ******************************************************************
! *********** written by GIOVANNI F. GRONCHI (Oct.2004) ************
! last modified 22/06/2005 GFG
! ==================================================================
!    CALL compute_critical_points_shift(ecpl%coord(1:5),eccom%coord(1:5), &
!         & fpl,fcom,nstat,nummin,nummax, &
!         & answer,warnflag,sflag,morse,weier,hzflag,hwflag,multfl)

  SUBROUTINE compute_critical_points_shift(elpl,elcom,fpl,fcom,nroots, &
       & nummin,nummax,answer,warnflag,sflag,morse,weier,hzflag,hwflag, & 
       & multfl,hevalflag)
    USE critical_points_shift
    USE fund_const                                      
    USE output_control
    IMPLICIT NONE
!    INTEGER, PARAMETER :: poldeg = 16 ! polynomial degree
! cometary elements of the planet and of the comet
    REAL(KIND=8),INTENT(IN),DIMENSION(5) :: elpl
    REAL(KIND=8),INTENT(IN),DIMENSION(5) :: elcom 
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
! warnflag(2) = false   higher degree terms \tilde{r}(t) are not small
! warnflag(3) = true    OK                                           
! warnflag(3) = false   low precision for \tilde{r}(t) coefficients     
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
    LOGICAL,INTENT(INOUT) :: hevalflag ! hevalflag = .true. OK!
                                       ! hevalflag = .false. unsuccessful 
                                       !                     hessian evaluation
! =========== end interface ========================================
    REAL(KIND=8),DIMENSION(5) :: elc1,elc2 
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
    REAL(KIND=8) :: alpha,sa,ca,beta,sb,cb
!    REAL(KIND=8) :: fpltil,sfpltil,cfpltil,fcomtil,sfcomtil,cfcomtil
    REAL(KIND=8) :: sx,cx,sy,cy
! Sylvester matrix elements                                         
    REAL(KIND=8) :: evp0(N),evp1(N),evp2(N)
    REAL(KIND=8) :: evq0(N),evq1(N),evq2(N),evq3(N),evq4(N)
    REAL(KIND=8) :: evr31(N),evr32(N),evr33(N),evr34(N),evr35(N)
    REAL(KIND=8) :: evr36(N),evr41(N),evr45(N),evr46(N)
    COMPLEX(KIND=8) :: cp0,cp1,cp2 
    COMPLEX(KIND=8) :: cq0,cq1,cq2,cq3,cq4
    COMPLEX(KIND=8) :: cr31,cr32,cr33,cr34,cr35
    COMPLEX(KIND=8) :: cr36,cr41,cr45,cr46

! constant term of det(S_1)
    COMPLEX(KIND=8) :: S10(N,N),detS10
! complex variables                                                  
    COMPLEX(KIND=8) :: complp0(N),complp1(N),complp2(N)
    COMPLEX(KIND=8) :: complq0(N),complq1(N),complq2(N),complq3(N),complq4(N)
    COMPLEX(KIND=8) :: complr31(N),complr32(N),complr33(N)
    COMPLEX(KIND=8) :: complr34(N),complr35(N)
    COMPLEX(KIND=8) :: complr36(N),complr41(N),complr45(N),complr46(N)
! Sylvester matrix 
    REAL(KIND=8) :: SYLV(6,6,N)
! complex variables                                                 
    COMPLEX(KIND=8) :: evalSYLV(6,6,N),evalSYLVj(6,6) 
    COMPLEX(KIND=8) :: complSYLVline(6) 
    COMPLEX(KIND=8) :: detevalSYLV,deteval(N) ! for the determinant 
! N-th roots of unity
    REAL(KIND=8) :: puroot(N),puroottil(N)
    COMPLEX(KIND=8) :: compluroot(N)

    COMPLEX(KIND=8) :: z_1,z_2,z_3,z_4,z_5,z_6
    COMPLEX(KIND=8) :: w_1,w_2,w_3,w_4,w_5,w_6

    REAL(KIND=8) :: s,t ! polynomial system variables 
    INTEGER :: count ! counter
    INTEGER :: j,l,nr ! loop indexes  
    INTEGER :: h,k ! matrix indexes 
    INTEGER :: ncoe,leftcoe ! indexes for the powers of the resultant
    INTEGER :: ans !  ans = 1: maximum 
!                     ans = 0: saddle 
!                     ans= -1: minimum 
!                     ans= -2: cannot decide
    SAVE count,alpha,beta
    SAVE evp0,evp1,evp2,evq0,evq1,evq2,evq3,evq4
    SAVE evr31,evr32,evr33,evr34,evr35,evr36,evr41,evr45,evr46
    SAVE complp0,complp1,complp2,complq0,complq1,complq2,complq3,complq4
    SAVE complr31,complr32,complr33,complr34,complr35,complr36
    SAVE complr41,complr45,complr46
    SAVE detS10,compluroot
! ***** FINIRE DI CONTROLLARE IL SAVING DELLE VARIABILI LOCALI *****
! ***** OPPURE SCRIVERE 'SAVE' (per tutte) *****
!   ==================================================================

! angular shift initialization
!    alpha = 0.d0
!    beta = 0.d0
    alpha = 0.5d0
    beta = 0.5d0

    count = 0 ! counter initialization
!    count = 9 ! to skip the shift

    elc1(1:5)=elpl(1:5)
    elc2(1:5)=elcom(1:5)
!    write(*,*)'ec1 elements:',elc1(1:5)
!    write(*,*)'ec2 elements:',elc2(1:5)

    goto 10 ! jump writing the following roots
    if(verb_moid.ge.20) then
       write(*,*)'1.d0-elc1(2)*cos(alpha)=',1.d0-elc1(2)*cos(alpha)

       z_1= (elc1(2)*sin(alpha)+sqrt(elc1(2)**2-1.d0))/(1.d0-elc1(2)*cos(alpha))
       write(*,*)'z_1=',z_1
       write(*,*)'z_2=',(elc1(2)*sin(alpha)+sqrt(elc1(2)**2-1.d0))/(1.d0-elc1(2)*cos(alpha))
       write(*,*)'z_3=',(cos(alpha)-1.d0)/sin(alpha)
       write(*,*)'z_4=',(cos(alpha)+1.d0)/sin(alpha)
       write(*,*)'z_5=',(sin(alpha)+sqrt(1.d0-elc1(2)**2))/(elc1(2)+cos(alpha))
       write(*,*)'z_6=',(sin(alpha)-sqrt(1.d0-elc1(2)**2))/(elc1(2)+cos(alpha))
       
       write(*,*)'w_1=',(elc2(2)*sin(beta)+sqrt(elc2(2)**2-1.d0)) &
            & /(1.d0-elc2(2)*cos(beta))
       write(*,*)'w_2=',(elc2(2)*sin(beta)+sqrt(elc2(2)**2-1.d0)) &
            & /(1.d0-elc2(2)*cos(beta))
       write(*,*)'w_3=',(cos(beta)-1.d0)/sin(beta)
       write(*,*)'w_4=',(cos(beta)+1.d0)/sin(beta)
       write(*,*)'w_5=',(sin(beta)+sqrt(1.d0-elc2(2)**2))/(elc2(2)+cos(beta))
       write(*,*)'w_6=',(sin(beta)-sqrt(1.d0-elc2(2)**2))/(elc2(2)+cos(beta))
       
       write(*,*)'z_+=',(cos(alpha)+1.d0)/sin(alpha)
       write(*,*)'w_+=',(cos(beta)+1.d0)/sin(beta)
       
       write(*,*)'z_-=',(cos(alpha)-1.d0)/sin(alpha)
       write(*,*)'w_-=',(cos(beta)-1.d0)/sin(beta)
    endif

10  CONTINUE

    CALL orbitcoe_shift(alpha,beta,elc1,elc2)

! initialization
    evp0(1:N) = 0.d0
    evp1(1:N) = 0.d0
    evp2(1:N) = 0.d0
    evq0(1:N) = 0.d0
    evq1(1:N) = 0.d0
    evq2(1:N) = 0.d0
    evq3(1:N) = 0.d0
    evq4(1:N) = 0.d0
    evr31(1:N) = 0.d0
    evr32(1:N) = 0.d0
    evr33(1:N) = 0.d0
    evr34(1:N) = 0.d0
    evr35(1:N) = 0.d0
    evr36(1:N) = 0.d0
    evr41(1:N) = 0.d0
    evr45(1:N) = 0.d0
    evr46(1:N) = 0.d0

! initialization for polynomial r(t)
    polycoe(1:N)=0.d0 
    polycoe1(1:N+1) = 0.d0
    evalrt(1:N)=0.d0 
    testevalpoly(1:N)=0.d0 

    nummin = 0 
    nummax = 0 
                                                                        
! ******************
! READ MATRIX DATA
!    call matrixdat_ta_shift(alpha,beta,N,p0,p1,p2,q0,q1,q2,q3,q4, &
!         & r31,r32,r33,r34,r35,r36,r41,r45,r46)
    CALL matrixdat_ta_shift(alpha,beta)

!    write(*,*)'p0(1:5)',p0(1:5)
!    write(*,*)'p1(1:5)',p1(1:5)
!    write(*,*)'p2(1:5)',p2(1:5)
!    write(*,*)'q0(1:3)',q0(1:3)
!    write(*,*)'q1(1:3)',q1(1:3)
!    write(*,*)'q2(1:3)',q2(1:3)
!    write(*,*)'q3(1:3)',q3(1:3)
!    write(*,*)'q4(1:3)',q4(1:3)

! remember p0,p1,p2,q0,q1,q2,q3,q4 and use
! evp0,evp1,evp2,evq0,evq1,evq2,evq3,evq4 as
! in/out variables
    evp0(1:5) = p0(1:5)
    evp1(1:5) = p1(1:5)
    evp2(1:5) = p2(1:5)
    evq0(1:5) = q0(1:5)
    evq1(1:5) = q1(1:5)
    evq2(1:5) = q2(1:5)
    evq3(1:5) = q3(1:5)
    evq4(1:5) = q4(1:5)    
    evr31(1:5) = r31(1:5)
    evr32(1:5) = r32(1:5)
    evr33(1:5) = r33(1:5)
    evr34(1:5) = r34(1:5)
    evr35(1:5) = r35(1:5)
    evr36(1:5) = r36(1:5)
    evr41(1:5) = r41(1:5)
    evr45(1:5) = r45(1:5)
    evr46(1:5) = r46(1:5)

! *******************************
! COMPUTE CONSTANT TERM of r(t)
    cp0 = evp0(1)
    cp1 = evp1(1)
    cp2 = evp2(1)
    cq0 = evq0(1)
    cq1 = evq1(1)
    cq2 = evq2(1)
    cq3 = evq3(1)
    cq4 = evq4(1)
    cr31 = evr31(1) 
    cr32 = evr32(1) 
    cr33 = evr33(1) 
    cr34 = evr34(1) 
    cr35 = evr35(1) 
    cr36 = evr36(1) 
    cr41 = evr41(1) 
    cr45 = evr45(1) 
    cr46 = evr46(1) 

!    write(*,*)'cp0,cp1,cp2,cq0,cq1,cq2,cq3,cq4:',cp0,cp1, &
!    & cp2,cq0,cq1,cq2,cq3,cq4
!    write(*,*)'cr31,cr32,cr33,cr34,cr35,cr36,cr41,cr45,cr46:', &
!    & cr31,cr32,cr33,cr34,cr35,cr36,cr41,cr45,cr46

    CALL compmodsylv16_shift(cp0,cp1,cp2,cq0,cq1,cq2,cq3,cq4,cr31,cr32,cr33, &
         & cr34,cr35,cr36,cr41,cr45,cr46,S10)

    CALL cdetcomp(S10,detS10)
!    write(*,*)'r0=',detS10

! *****************************************************
! EVALUATE the COEFFICIENTS of the MATRIX \tilde{S}_1
    CALL rvfft(evp0,N,expo) 
    CALL rvfft(evp1,N,expo) 
    CALL rvfft(evp2,N,expo) 
    CALL rvfft(evq0,N,expo) 
    CALL rvfft(evq1,N,expo) 
    CALL rvfft(evq3,N,expo) 
    CALL rvfft(evq4,N,expo) 
    CALL rvfft(evr31,N,expo) 
    CALL rvfft(evr32,N,expo)
    CALL rvfft(evr33,N,expo)
    CALL rvfft(evr34,N,expo)
    CALL rvfft(evr35,N,expo) ! e' costante, evitare la valutazione!
    CALL rvfft(evr36,N,expo) ! e' costante, evitare la valutazione!
    CALL rvfft(evr41,N,expo) 
    CALL rvfft(evr45,N,expo) ! e' costante, evitare la valutazione!
    CALL rvfft(evr46,N,expo) ! e' costante, evitare la valutazione!

    complp0(1) = DCMPLX(evp0(1),0.d0)
    complp1(1) = DCMPLX(evp1(1),0.d0)
    complp2(1) = DCMPLX(evp2(1),0.d0)
    complq0(1) = DCMPLX(evq0(1),0.d0)
    complq1(1) = DCMPLX(evq1(1),0.d0)
    complq2(1) = DCMPLX(evq2(1),0.d0)
    complq3(1) = DCMPLX(evq3(1),0.d0)
    complq4(1) = DCMPLX(evq4(1),0.d0) 
    complr31(1) = DCMPLX(evr31(1),0.d0) 
    complr32(1) = DCMPLX(evr32(1),0.d0) 
    complr33(1) = DCMPLX(evr33(1),0.d0) 
    complr34(1) = DCMPLX(evr34(1),0.d0) 
    complr35(1) = DCMPLX(evr35(1),0.d0) 
    complr36(1) = DCMPLX(evr36(1),0.d0) 
    complr41(1) = DCMPLX(evr41(1),0.d0) 
    complr45(1) = DCMPLX(evr45(1),0.d0) 
    complr46(1) = DCMPLX(evr46(1),0.d0) 
!
    complp0(N/2+1) = DCMPLX(evp0(N/2+1),0.d0)
    complp1(N/2+1) = DCMPLX(evp1(N/2+1),0.d0)
    complp2(N/2+1) = DCMPLX(evp2(N/2+1),0.d0)
    complq0(N/2+1) = DCMPLX(evq0(N/2+1),0.d0)
    complq1(N/2+1) = DCMPLX(evq1(N/2+1),0.d0)
    complq2(N/2+1) = DCMPLX(evq2(N/2+1),0.d0)
    complq3(N/2+1) = DCMPLX(evq3(N/2+1),0.d0)
    complq4(N/2+1) = DCMPLX(evq4(N/2+1),0.d0)
    complr31(N/2+1) = DCMPLX(evr31(N/2+1),0.d0) 
    complr32(N/2+1) = DCMPLX(evr32(N/2+1),0.d0) 
    complr33(N/2+1) = DCMPLX(evr33(N/2+1),0.d0) 
    complr34(N/2+1) = DCMPLX(evr34(N/2+1),0.d0) 
    complr35(N/2+1) = DCMPLX(evr35(N/2+1),0.d0) 
    complr36(N/2+1) = DCMPLX(evr36(N/2+1),0.d0) 
    complr41(N/2+1) = DCMPLX(evr41(N/2+1),0.d0) 
    complr45(N/2+1) = DCMPLX(evr45(N/2+1),0.d0) 
    complr46(N/2+1) = DCMPLX(evr46(N/2+1),0.d0) 
!
    DO j = 1,N/2-1
       complp0(j+1) = DCMPLX(evp0(j+1),evp0(N-j+1))
       complp0(N/2+j+1) = DCMPLX(evp0(N/2-j+1),-evp0(N/2+1+j))        

       complp1(j+1) = DCMPLX(evp1(j+1),evp1(N-j+1))
       complp1(N/2+j+1) = DCMPLX(evp1(N/2-j+1),-evp1(N/2+1+j))        

       complp2(j+1) = DCMPLX(evp2(j+1),evp2(N-j+1))
       complp2(N/2+j+1) = DCMPLX(evp2(N/2-j+1),-evp2(N/2+1+j))        

       complq0(j+1) = DCMPLX(evq0(j+1),evq0(N-j+1))
       complq0(N/2+j+1) = DCMPLX(evq0(N/2-j+1),-evq0(N/2+1+j))        

       complq1(j+1) = DCMPLX(evq1(j+1),evq1(N-j+1))
       complq1(N/2+j+1) = DCMPLX(evq1(N/2-j+1),-evq1(N/2+1+j))        

       complq2(j+1) = DCMPLX(evq2(j+1),evq2(N-j+1))
       complq2(N/2+j+1) = DCMPLX(evq2(N/2-j+1),-evq2(N/2+1+j))        

       complq3(j+1) = DCMPLX(evq3(j+1),evq3(N-j+1))
       complq3(N/2+j+1) = DCMPLX(evq3(N/2-j+1),-evq3(N/2+1+j))        

       complq4(j+1) = DCMPLX(evq4(j+1),evq4(N-j+1))
       complq4(N/2+j+1) = DCMPLX(evq4(N/2-j+1),-evq4(N/2+1+j))        

       complr31(j+1) = DCMPLX(evr31(j+1),evr31(N-j+1))
       complr31(N/2+j+1) = DCMPLX(evr31(N/2-j+1),-evr31(N/2+1+j))        

       complr32(j+1) = DCMPLX(evr32(j+1),evr32(N-j+1))
       complr32(N/2+j+1) = DCMPLX(evr32(N/2-j+1),-evr32(N/2+1+j))
        
       complr33(j+1) = DCMPLX(evr33(j+1),evr33(N-j+1))
       complr33(N/2+j+1) = DCMPLX(evr33(N/2-j+1),-evr33(N/2+1+j))        

       complr34(j+1) = DCMPLX(evr34(j+1),evr34(N-j+1))
       complr34(N/2+j+1) = DCMPLX(evr34(N/2-j+1),-evr34(N/2+1+j))        

       complr35(j+1) = DCMPLX(evr35(j+1),evr35(N-j+1))
       complr35(N/2+j+1) = DCMPLX(evr35(N/2-j+1),-evr35(N/2+1+j))        

       complr36(j+1) = DCMPLX(evr36(j+1),evr36(N-j+1))
       complr36(N/2+j+1) = DCMPLX(evr36(N/2-j+1),-evr36(N/2+1+j))        

       complr41(j+1) = DCMPLX(evr41(j+1),evr41(N-j+1))
       complr41(N/2+j+1) = DCMPLX(evr41(N/2-j+1),-evr41(N/2+1+j))        

       complr45(j+1) = DCMPLX(evr45(j+1),evr45(N-j+1))
       complr45(N/2+j+1) = DCMPLX(evr45(N/2-j+1),-evr45(N/2+1+j))        

       complr46(j+1) = DCMPLX(evr46(j+1),evr46(N-j+1))
       complr46(N/2+j+1) = DCMPLX(evr46(N/2-j+1),-evr46(N/2+1+j))        
    ENDDO

! ****************************************
! p(x)=x POLYNOMIAL (for roots of unity)
    puroot(1) = 0.d0
    puroot(2) = 1.d0
    DO j = 3,N
       puroot(j) = 0.d0
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
!    write(*,*)'roots of 1' !selected clockwise
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
!    write(*,*)'detS10',detS10
!    write(*,*)'deteval:',deteval(1:N)
!    write(*,*)'compluroot:',compluroot(1:N)
!    write(*,*)'15 deg poly evals:',evpol15(1:N)

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
          WRITE(*,*)'WARNING! very small leading coefficient '
          WRITE(*,*)'of tilde{r} polynomial (< 1D-15)'
          WRITE(*,*)'leading coe:',evalpoly(poldeg)
       endif
       GOTO 13
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
!   get a monic polynomial (polycoe1(poldeg+1) = 1)     
    IF(abs(polycoe1(poldeg+1)).ge.1.d-5) THEN
       DO ncoe = 1,poldeg+1 
          polycoe1(ncoe) = polycoe1(ncoe)/polycoe1(poldeg+1) 
       ENDDO
! check
!       write(*,*)'r(t) coefficients after normalization'
!       write(*,101) polycoe1(1:poldeg+1)
    ELSE
! leave the coeficients unnormalized
       if(verb_moid.ge.20)then
          write(*,*)'leave the coefficients unnormalized'
          write(*,*)'leading coefficient:',polycoe1(poldeg+1)
       endif
    ENDIF

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
! ******************************************************************
111 CONTINUE
                                                                        
! ******************************************************************
! FIND REAL ROOTS OF THE POLYNOMIAL res[ P1(t,s),P2(t,s),t ]        
    CALL solvpoly(poldeg,polycoe1(1:poldeg+1),wzero,nsol,hzflag, &
         & multfl)
! ******************************************************************
    if(.not.hzflag) then
       if(verb_moid.ge.20) then
          write(*,*)'warning: hzflag=',hzflag
       endif
       goto 13
    endif
 
! check
    if(verb_moid.ge.20) then
       write(*,*)'nroots=',nsol
       write(*,*)'wzero:',wzero(1:nsol)
    endif

! ************************************************************
! COMPUTE FOR EACH w THE CORRESPONDING VALUE of z                   
    CALL solvesystem_ta_shift(alpha,beta,nsol,wzero(1:nsol),zzero(1:nsol),sflag,hwflag) 
! ************************************************************
    if(.not.hwflag) then
       goto 13
    endif
!    write(*,*)'zzero',zzero(1:nsol)
!    write(*,*)'alpha,beta',alpha,beta
    sa = sin(alpha)
    ca = cos(alpha)
    sb = sin(beta)
    cb = cos(beta)
    
    DO nr = 1,nsol 

! conversion from (z,w) to the true anomalies (Vpl,vcom)
       sx =  2.d0*zzero(nr)/(1.d0+(zzero(nr))**2) 
       cx = (1.d0-(zzero(nr))**2)/(1.d0+(zzero(nr))**2) 
       sfpl(nr) = sx*ca + cx*sa
       cfpl(nr) = cx*ca - sx*sa

       sy = 2.d0*wzero(nr)/(1.d0+(wzero(nr))**2) 
       cy = (1.d0-(wzero(nr))**2)/(1.d0+(wzero(nr))**2) 
       sfcom(nr) = sy*cb + cy*sb
       cfcom(nr) = cy*cb - sy*sb

! critical points (true anomalies)                           
       fpl(nr) = datan2(sfpl(nr),cfpl(nr)) 
       fcom(nr) = datan2(sfcom(nr),cfcom(nr)) 

! compute type of critical points
       CALL hess_ta(fpl(nr),fcom(nr),ans) 

!       ans = -2 !dummy to force int_eval
       IF(ans.eq.-2) THEN
          hevalflag = .false.
          WRITE(*,*)'cannot decide type of critical point: calling int_eval'
          CALL int_eval_ta(fpl(nr),fcom(nr),ans)
       ENDIF

       IF (ans.eq.-1) THEN 
          nummin = nummin + 1 
       ELSEIF (ans.eq.1)THEN 
          nummax = nummax + 1 
       ENDIF
                                                                       
       answer(nr) = ans 

! conversion into degrees                                           
       fpl(nr) = fpl(nr)*degrad 
       fcom(nr) =  fcom(nr)*degrad 
       
    ENDDO

! ***************************************************
! different checks if one orbit is cometary  ! STILL TO BE WRITTEN ....
    IF ((elc1(2).ge.1.d0).or.(elc2(2).ge.1.d0)) THEN
       GOTO 1300
    ENDIF

! *****************************************************             
! ============= CHECK WITH MORSE THEORY ===============             
    IF (2*(nummax+nummin).ne.nsol) THEN 
       morse = .false.
    ENDIF
! ============ CHECK WITH WEIERSTRASS THEOREM =========             
    IF((nummax.lt.1).or.(nummin.lt.1))THEN 
       weier = .false. 
    ENDIF
! *****************************************************             

! dummy: to force the shift
!    morse=.false.
! **********************************************************************
! ******************  A N G U L A R   S H I F T  ***********************
! ************** V = \Xi + \alpha  ;  v = \xi + \beta ******************
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
       alpha = alpha + dpig/16.d0
       beta = beta + dpig/15.d0
          if(verb_moid.ge.20) then
             WRITE(*,*)'APPLY ANGULAR SHIFT: count=',count
             WRITE(*,*)'               alpha,beta:',alpha,beta
          endif
! restoring flags
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
       hevalflag = .true.
       GOTO 10 
    ELSE
       if(verb_moid.ge.20) then
          WRITE(*,*)'computation OK! count = ',count
          WRITE(*,*)'morse,weier,warnflag(1)',morse,weier,warnflag(1)
          WRITE(*,*)'sflag(2:6)',sflag(2:6)
          WRITE(*,*)'hzflag',hzflag
       endif
! continue the computation
    ENDIF

! maximal number of angular shifts
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

1300 CONTINUE ! to skip computation in case of failure
    
    nroots = nsol
    
    RETURN 
  END SUBROUTINE compute_critical_points_shift

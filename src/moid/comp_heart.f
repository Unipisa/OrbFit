c     ********************************************************************
c     ******** HEART of the COMPUTATION of the STATIONARY POINTS *********
c     **************** (including HESSIAN COMPUTATION) *******************
c     ********************************************************************
c     *********** written by GIOVANNI F. GRONCHI (2001) ******************
c     ********** Department of Mathematics, UNIVERSITY of PISA ***********
c     ====================================================================
      SUBROUTINE comp_heart(u,upl,nroots,nummin,nummax,answer,
     *     warnflag,sflag)
      IMPLICIT NONE
c     eccentric anomalies
      DOUBLE PRECISION u(20),upl(20)
c     number of roots, minima and maxima
      INTEGER nroots,nummin,nummax
c     answer= 1:maximum 
c     answer= 0:saddle
c     answer=-1:minimum
      INTEGER answer(20)
c     flags:
c     general warning flags
c     warnflag(1) = 1: OK ; 
c     warnflag(1) = 0: leading coefficient of resultant is very small; 
c     warnflag(2) = 1: OK ; 
c     warnflag(2) = 0: higher degree terms in resultant are not small; 
c     warnflag(3) = 1: OK ; 
c     warnflag(3) = 0: low precision for resultant coefficients (fft); 
      INTEGER warnflag(3)
c     solving system flags
c     sflag(1) = 1: OK! ;  
c     sflag(1) = 0: WARNING! They are both good solutions!;
c     sflag(2) = 1: OK! ;  
c     sflag(2) = 0: WARNING! Neither of two is close to 0!
      INTEGER sflag(2)
c     -------------------------------------- end interface ---------------
c     RESULTANT EVALUATIONS (32) and COEFFICIENTS
      DOUBLE PRECISION polycoe(32),evalpoly(32) 
c     for a test of the coefficient found
      DOUBLE PRECISION evalres(32),testevalpoly(32)
c     ROOTS of the RESULTANT
      DOUBLE PRECISION tzero(20)
c     error estimate
      DOUBLE PRECISION radius(20)
c     corresponding value of s
      DOUBLE PRECISION szero(20)
c     for the eccentric anomalies
      DOUBLE PRECISION su(20),cu(20)
      DOUBLE PRECISION supl(20),cupl(20)
c     Sylvester matrix elements 
      DOUBLE PRECISION alpha(32),beta(32),gamma(32)
      DOUBLE PRECISION A(32),B(32),C(32),D(32),E(32)
      DOUBLE PRECISION alphatil(32),betatil(32),gammatil(32)
      DOUBLE PRECISION Atil(32),Btil(32),Ctil(32),Dtil(32),Etil(32)
c     complex variable
      COMPLEX*16 complalpha(32),complbeta(32),complgamma(32)
      COMPLEX*16 complA(32),complB(32),complC(32)
      COMPLEX*16 complD(32),complE(32)
c     Sylvester matrix
      DOUBLE PRECISION SYLV(6,6,32)
c     complex variables
      COMPLEX*16 evalSYLV(6,6,32),evalSYLVj(6,6)
      COMPLEX*16 complSYLVline(6)
c     for the determinant
      COMPLEX*16 detevalSYLV,deteval(32)
c     complex zeros of the polynomial
      COMPLEX*16 zr1(20)
c     polynomial system variables
      DOUBLE PRECISION s,t
c     loop indexes
      INTEGER j,l,nr
c     matrix indexes
      INTEGER h,k
c     index for the powers of the resultant
      INTEGER ncoe,leftcoe
c     number of evaluations and degree of res(t)     
      INTEGER N,poldeg
c     ans=1:maximum; ans=0:saddle; ans=-1:minimum
      INTEGER ans
c     trig constants
      INCLUDE 'trig.h'
c     ======================================================================
      
c     degree of the polynomial
      poldeg = 20
c     number of evaluations
      N = 32

c     ====================================================================
c     initialization of coefficients of polynomial res[P1(t,s),P2(t,s),t] 
c     (must be set to zero also from 22-th to 32-th for testing
c     them again with DFT algorithm)
c     ====================================================================
      DO l = 1,N
         polycoe(l)=0.d0
         evalres(l)=0.d0
         testevalpoly(l)=0.d0
      ENDDO
c     
      nummin = 0
      nummax = 0

c     ********************************************************************
c     ====================== READING MATRIX DATA =======================
c     ********************************************************************
      call matrixdat(N,alpha,beta,gamma,A,B,C,D,E)

c     ********************************************************************
c     DEFINING alphatil,betatil,gammatil,Atil,Btil,Ctil,Dtil,Etil
c     ********************************************************************
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
      
      
c     ********************************************************************
c     ========== EVALUATING the COEFFICIENTS of the MATRIX ===========
c     ********************************************************************
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
c     ========== C = 0 ===========
      CALL rvfft(C,32,5)
      CALL rvfft(Ctil,32,5)
c     ============================
      CALL rvfft(D,32,5)
      CALL rvfft(Dtil,32,5)
      CALL rvfft(E,32,5)
      CALL rvfft(Etil,32,5)
c     =================================================
      
c     *********************************************************************
c     ================ DECODING DFT OUTPUT  ===============================
c     *********************************************************************
      CALL decode_out(N,alpha,alphatil,complalpha)
      CALL decode_out(N,beta,betatil,complbeta)
      CALL decode_out(N,gamma,gammatil,complgamma)
      CALL decode_out(N,A,Atil,complA)
      CALL decode_out(N,B,Btil,complB)
      CALL decode_out(N,C,Ctil,complC)
      CALL decode_out(N,D,Dtil,complD)
      CALL decode_out(N,E,Etil,complE)
            
c     ******************************************************************
c     ====== BUILDING EVALUATED MATRIX evalSYLV(h,k,j)  =========
c     ******************************************************************
      DO j = 1,N
         CALL clevcompsylv(complalpha(j),complbeta(j),complgamma(j),
     *        complA(j),complB(j),complC(j),complD(j),
     *        complE(j),evalSYLVj)
         DO h = 1,6
            DO k = 1,6
c     building evaluated matrixes
               evalSYLV(h,k,j) = evalSYLVj(h,k)        
            ENDDO
         ENDDO
      ENDDO
      
c     **********************************************************************
c     COMPUTATION OF THE EVALUATIONS OF THE RESULTANT IN THE N POINTS:
c     **********************************************************************
c     BEGIN MAIN DO LOOP
      DO j = 1,N
         
c     ==== building matrixes whose determinant is to compute ===
         DO h = 1,6
            DO k = 1,6
               evalSYLVj(h,k)=evalSYLV(h,k,j)
            ENDDO
         ENDDO
c     ==========================================================
         
c     **********************************************************************
c     for a given value of j (1<j<32) compute det(evalSYLV(h,k,j))  
c     **********************************************************************
         CALL cdetcomp(evalSYLVj,detevalSYLV)
         
c     ===================================================
c     EVALUATIONS of RESULTANT in the N-th ROOTS of UNITY
c     ===================================================
         deteval(j)=detevalSYLV   
c     ===================================================
c     ===================================================

c     END MAIN DO LOOP
      ENDDO
      
c     *********************************************************************
c     =========== COMPUTE COEFFICIENTS of RESULTANT res(t) ==========
c     *********************************************************************
      
c     CODING INPUT ACCORDING TO CONVENTION
      CALL code_inp(N,deteval,evalpoly)

c     ======== keep in memory the evaluations of the =====
c     ======== resultant for the evaluation test ========= 
      DO j = 1,N
         evalres(j)=evalpoly(j)
      ENDDO
c     ==========================================
      
c     COMPUTE COEFFICIENTS
      CALL irvfft(evalpoly,32,5)
      
c     ============================================================
c     compute coefficient normalizing them (polycoe(21) = 1)
      IF (abs(evalpoly(poldeg+1)).le.1d-9) THEN         
         warnflag(1) = 0
         WRITE(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         WRITE(*,*)'!!! LEADING COEFFICIENT = 0 !!!'
         WRITE(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

C     DO NOT CONTINUE THE COMPUTATION: NROOTS = -1 NUMMIN = -1 NUMMAX = -1
            nroots = -1
            nummin = -1
            nummax = -1
            GOTO 1300

      ELSE
c     GET A MONIC POLYNOMIAL
         DO ncoe = 1,poldeg+1
            polycoe(ncoe) = evalpoly(ncoe)/evalpoly(poldeg+1) 
c     ===============================================================
c     keep in memory the unnormalized coefficients of the resultant
c     for the test
            testevalpoly(ncoe) = evalpoly(ncoe)
c     ===============================================================
         ENDDO
      ENDIF

c     *********************************************************************
c     +++++++++++++++++++++++++++++ T E S T +++++++++++++++++++++++++++++++
c     *********************************************************************
c     CHECK THE UNNORMALIZED COEFFICIENTS OF THE HIGHER 
c     THAN 20 DEGREE TERMS 
      DO ncoe = poldeg+2,N
         IF (abs(evalpoly(ncoe)).gt.1.d-5) THEN
            warnflag(2)=0
c            WRITE(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!'
c            WRITE(*,*)'!!!!!!!! WARNING !!!!!!!!!!'
c            WRITE(*,*)'!!! HIGHER DEGREE TERMS !!!'
c            WRITE(*,*)'!!!!!!!! NOT SMALL !!!!!!!!'
c            WRITE(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!'
c            WRITE(*,*)'coeff(',ncoe,')=',evalpoly(ncoe)
c            WRITE(*,*)'==============================================='
          ENDIF
      ENDDO

c     *********************************************************************
c     +++++++++++++++++++++++++++++ T E S T +++++++++++++++++++++++++++++++
c     *********************************************************************
c     EVALUATING THE POLYNOMIAL RESULTANT OBTAINED REMOVING THE 11
c     COEFFICIENT FROM DEGREE 32 UP TO DEGREE 22 INCLUDED 
c     *********************************************************************
      CALL rvfft(testevalpoly,32,5)
      DO j = 1,N
c     =========================== CHECK ============================
         IF ( abs(testevalpoly(j)-evalres(j)).gt.1.d-5) THEN
c     HINT: polycoe(j) are coded, so we have to test them with evalres(j)
            warnflag(3) = 0
c            WRITE(*,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
c            WRITE(*,*)'ERROR: DIFFERENCE IN THE EVALUATIONS OF' 
c            WRITE(*,*)'THE RESULTANT IS NOT SMALL!!!:'
c            WRITE(*,*)'diff=',abs(testevalpoly(j)-evalres(j))
c            WRITE(*,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
         ENDIF
      ENDDO
c     *********************************************************************
c     *********************************************************************
      
c     *********************************************************************
c     FIND REAL ROOTS OF THE POLYNOMIAL res[ P1(t,s),P2(t,s),t ] 
c     *********************************************************************
      CALL solv20(polycoe,tzero,nroots,radius,zr1)
      
c     check number of roots
c      write(*,*)'NUMBER of ROOTS=',nroots

c     *********************************************************************
c     COMPUTE FOR EACH t THE CORRESPONDING VALUE of s 
c     *********************************************************************
      CALL solvesystem(nroots,tzero,szero,sflag)

      DO nr = 1,nroots
         
c     *********************************************************************
c     CONVERSION OF THE SOLUTIONS (t,s) IN TERMS 
c     OF THE ECCENTRIC ANOMALIES (u,upl)
c     *********************************************************************
         su(nr) = 2.d0*tzero(nr)/(1.d0+(tzero(nr))**2) 
         cu(nr) = (1.d0-(tzero(nr))**2)/(1.d0+(tzero(nr))**2) 
         supl(nr) =  2.d0*szero(nr)/(1.d0+(szero(nr))**2) 
         cupl(nr) = (1.d0-(szero(nr))**2)/(1.d0+(szero(nr))**2) 
c     stationary points (eccentric anomalies)
         u(nr) = datan2(su(nr),cu(nr))
         upl(nr) = datan2(supl(nr),cupl(nr))
         
c     *********************************************************************
c     SELECT AMONG MINIMA, MAXIMA AND SADDLE
         CALL hess(u(nr),upl(nr),ans)
         IF (ans.eq.-1) THEN
            nummin = nummin + 1
         ELSEIF (ans.eq.1)THEN
            nummax = nummax + 1
         ENDIF
c
         answer(nr) = ans
c     conversion into degrees
         u(nr) =  u(nr)*degrad
         upl(nr) = upl(nr)*degrad
         
      ENDDO
c     ********************************************************************

c     DEGENERATE CASE
 1300 CONTINUE

      RETURN
      END

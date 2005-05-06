! ********************************************************************
! **************  PROGRAM   S T A T S E A R C H    *******************
! *********** written by GIOVANNI F. GRONCHI (May 2004) **************
! *************** E-MAIL gronchi@dm.unipi.it *************************
! ********************************************************************
  PROGRAM stat_search
    USE fund_const
    USE output_control
    IMPLICIT NONE
! ==================== 1st ELLIPSE =======================
    DOUBLE PRECISION :: apl,epl,Ipl,ompl,Omnodpl,effepl
! ==================== 2nd ELLIPSE =======================
    DOUBLE PRECISION :: a,e,Inc,om,Omnod,effe
! auxiliary variables
    DOUBLE PRECISION :: ap,ep,Ip,opp,Onp,aa,ea,Ia,opa,Ona
! for compute_statpts
! elements of the Earth and of the asteroid (angles in deg.)
    DOUBLE PRECISION :: ekpl(6),elkep(6) 
! eccentric anomalies 
    DOUBLE PRECISION, DIMENSION(20) :: u,upl
! SQUARED DISTANCE function 
    DOUBLE PRECISION :: D2(20)
! number of stationary points (maximum = 20)
    INTEGER :: nstat
! number of relative minima/maxima found
    INTEGER :: nummin,nummax
! type of singular point : answer(j) =  1    MAXIMUM      
!                          answer(j) = -1    MINIMUM
!                          answer(j) =  0    SADDLE
!                          answer(j) =  2    CANNOT DECIDE
    INTEGER :: answer(20)
    CHARACTER(LEN=13) :: char
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

! treshold for writing data
    INTEGER :: treshold
! loop indexes and counters
    INTEGER :: h,i,j,k,l,m,n,o,p,q,count,ii,nfail
! ==================================================================

!    verb_moid=21

! CRITICAL POINTS FILE (for matlab)
    OPEN(4,file='critpts',status='unknown')

! MORSE-WEIERSTRASS CHECK FILE
    OPEN(3,file='morse_weier',status='unknown')

! HIGH NUMBER of STATIONARY POINTS  FILE
    OPEN(2,file='hnstat.out',status='unknown')

    WRITE(2,*)'+++++++++++++++++++++++++++++++++++++++++++++++++++++++&
         &++++++++++++'
    WRITE(2,*)'apl,epl,Ipl,ompl,Omnodpl'
    WRITE(2,*)'a,e,Inc,om,Omnod'
    WRITE(2,*)'nstat,nummin,nummax'
    WRITE(2,*)'+++++++++++++++++++++++++++++++++++++++++++++++++++++++&
         &++++++++++++'

! flags initialization                                              
    weier = .true.
! weier = .false. ! dummy
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

! initialization
    count = 0
    nfail = 0

    write(*,*) 'write integer treshold for stat points to write:'
    read(*,*) treshold
    
! READING ORBITAL DATA
    OPEN(1,file='orbelems.dat',status='old')
    READ(1,*)
    READ(1,*)a,e,Inc,om,Omnod,effe
    READ(1,*)
    READ(1,*)apl,epl,Ipl,ompl,Omnodpl,effepl
    CLOSE(1)
    
! 1st ELLIPSE
    aa = a
    ea = e
    Ia = Inc*radeg
    opa = om*radeg
    Ona = Omnod*radeg

! 2nd ELLIPSE
    ap = apl
    ep = epl
    Ip = Ipl*radeg
    opp = ompl*radeg
    Onp = Omnodpl*radeg


! =======  ENTER MAIN LOOP =======

! 1st SEMIMAJOR AXIS
    DO h = 0,0
!         WRITE(*,*)'h=', h
       elkep(1)   = aa + 1*h

! 1st ECCENTRICITY
       DO i = 0,99
!            WRITE(*,*)'i=', i
          elkep(2) = ea + 0.001*i

! 1st INCLINATION
          DO j = 0,0
!                WRITE(*,*)'j=', j
             elkep(3) = Ia + 0.1*j*radeg

! 1st OMNOD
             DO k = 0,0      
!                  WRITE(*,*)'k=', k
                elkep(4) = Ona + 4*k*radeg
! 1st OMEGA
                
                DO l = 0,0
!                   WRITE(*,*)'l=', l
                   elkep(5) = opa + 4*l*radeg
     
! 2nd SEMIMAJOR AXIS
                   DO m = 0,0
!                     WRITE(*,*)'m=', m
                      ekpl(1)   = ap + 1*m

! 2nd ECCENTRICITY
                      DO n = 0,0
!                        WRITE(*,*)'n=', n
                         ekpl(2) = ep + 0.1*n

! 2nd INCLINATION
                         DO o = 0,0
!                           WRITE(*,*)'o=',o
                            ekpl(3) = Ip + 0.1*o*radeg

! 2nd OMNOD
                            DO p = 0,0       
!                              WRITE(*,*)'p=', p
                               ekpl(4) = Onp + 4*p*radeg
! 2nd OMEGA
                
                               DO q = 0,0
!                                 WRITE(*,*)'l=', l
                                  ekpl(5) = opp + 4*q*radeg
                      
                                  count = count+1
                                  
                                  elkep(6) = 0.d0
                                  ekpl(6) = 0.d0

                                  write(*,*)'elkep'
                                  write(*,100)elkep(1:2),elkep(3:5)*degrad
                                  write(*,*)'ekpl'
                                  write(*,100)ekpl(1:2),ekpl(3:5)*degrad

! ==================================================================
! ******************************************************************
                      CALL compute_statpts(ekpl,elkep,u,upl,D2,nstat,&
                           & nummin,nummax,answer, &
                           & morse,weier,warnflag,sflag)
! ******************************************************************
! ==================================================================
                        
                      IF(nstat.ge.treshold) THEN
                         write(*,*)'nstat,nummax,nummin = ',nstat,&
                              & nummax,nummin 
                         WRITE(2,100)ekpl(1),ekpl(2),ekpl(3)*degrad,ekpl(4)*degrad,ekpl(5)*degrad
                         WRITE(2,100)elkep(1),elkep(2),elkep(3)*degrad,elkep(4)*degrad,elkep(5)*degrad
                         WRITE(2,101)nstat,nummin,nummax
!                         WRITE(4,*) nstat,0,0
                         DO ii=1,nstat
                            IF(answer(ii).eq.1) THEN   
                               char='MAXIMUM'      
                            ELSEIF(answer(ii).eq.-1) THEN 
                               char='MINIMUM'
                            ELSEIF(answer(ii).eq.0) THEN
                               char='SADDLE'
                            ELSEIF(answer(ii).eq.2) THEN
                               char='CANNOT DECIDE'
                            ENDIF
                            WRITE(2,102)u(ii),upl(ii),D2(ii),char
                            WRITE(*,102)u(ii),upl(ii),D2(ii),char
                            WRITE(4,103)u(ii),upl(ii),D2(ii),answer(ii)
                         ENDDO
                         WRITE(2,*)'NUMBER',count
                         WRITE(2,*)'*************************************'
100                      FORMAT(1x,5(f10.3,2x))
101                      FORMAT(i3,1x,i3,1x,i3)
102                      FORMAT(1x,3(f10.5,4x),a13)
103                      FORMAT(1x,3(f10.5,4x),i3)
                      ENDIF
                        
! ============= CHECK WITH MORSE THEORY ===============
                      IF (.not.morse) THEN
                         WRITE(*,*)'stat_search: ERROR! MORSE THEORY VIOLATION!'
                         WRITE(3,*)'stat_search: ERROR! MORSE THEORY VIOLATION!'   
                         WRITE(3,*)'FAILED COMPUTATION FOR DATA'
                         WRITE(3,100)ap,ep,Ip*degrad,opp*degrad, &
                              & Onp*degrad
                         WRITE(3,100) aa,ea,Ia*degrad,opa*degrad, &
                              & Ona*degrad
                         WRITE(3,*)'NUMSTAT=',nstat
                         WRITE(3,*)'NUMMIN=',nummin
                         WRITE(3,*)'NUMMAX=',nummax
                         nfail = nfail+1
                      ENDIF
! =====================================================
                        
! ============ CHECK WITH WEIERSTRASS THEOREM =========
                      IF(.not.weier)THEN
                         WRITE(*,*)'stat_search: ERROR! WEIERSTRASS THEOREM VIOLATION!'
                         WRITE(3,*)'stat_search: ERROR! WEIERSTRASS THEOREM VIOLATION!'
                         WRITE(3,*)'FAILED COMPUTATION FOR DATA'
                         WRITE(3,100)ap,ep,Ip*degrad,opp*degrad,Onp*degrad
                         WRITE(3,100) aa,ea,Ia*degrad,opa*degrad,Ona*degrad
                         IF(morse)THEN ! add new failure only if 
                                       ! not counted my 'Morse check' 
                            nfail = nfail + 1
                         ENDIF
                      ENDIF
                      
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO
ENDDO
ENDDO
ENDDO
ENDDO

! closing morse_weier
    CLOSE(3)
! closing hnstat.out
    CLOSE(2)
    write(*,*)'Number of processed orbits:',count
    write(*,*)'Number of failed computations:',nfail

! closing matlab file
    CLOSE(4)
    
  END PROGRAM stat_search

! ********************************************************************
! *******  PROGRAM   S T A T _ S E A R C H _ T A _ S H I F T  ********
! *********** written by GIOVANNI F. GRONCHI (Nov.2004) **************
! ********************************************************************
! last modified 24/06/2005 GFG
! ********************************************************************
  PROGRAM stat_search_ta_shift
    USE fund_const
    USE output_control
    USE critical_points_shift
!    USE ouput_control
    IMPLICIT NONE
    DOUBLE PRECISION :: perih1,e1,i1,om1,Omn1 !1st ELLIPSE 
    DOUBLE PRECISION :: perih2,e2,i2,om2,Omn2 !2nd ELLIPSE
!    INTEGER, PARAMETER :: poldeg=16
! cometary elements of the first and second orbit
    DOUBLE PRECISION :: ec1(6),ec2(6)
    DOUBLE PRECISION :: ph1step,e1step,i1step,Omn1step,om1step
    INTEGER :: ph1iter,e1iter,i1iter,Omn1iter,om1iter
    DOUBLE PRECISION :: ph2step,e2step,i2step,Omn2step,om2step
    INTEGER :: ph2iter,e2iter,i2iter,Omn2iter,om2iter
    CHARACTER*1 :: ansiter
! true anomalies
    DOUBLE PRECISION :: v1(poldeg),v2(poldeg),vv1,vv2,chv1,chv2
    DOUBLE PRECISION :: D2(poldeg),DD2 !squared distance function
    INTEGER :: D2_ord(poldeg)
    INTEGER :: nstat! number of critical points (maximum = poldeg)
    INTEGER :: nummin,nummax! number of relative minima/maxima
! type of singular point : answer(j) =  1    MAXIMUM      
!                          answer(j) = -1    MINIMUM
!                          answer(j) =  0    SADDLE
!                          answer(j) =  2    CANNOT DECIDE
    INTEGER :: answer(poldeg)
! error flags
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
    LOGICAL :: hevalflag ! hevalflag = .true. OK!
                         ! hevalflag = .false. unsuccessful hessian evaluation
    LOGICAL :: circ_copl_orb ! circ_copl_orb = .false. OK!
                             ! circ_copl_orb = .true. circular coplanar orbits
    LOGICAL :: overlap_ell ! overlap_ell = .false. OK!
                           ! overlap_ell = .true. overlapping ellipses
    LOGICAL :: write_on_screen ! write_on_screen = .false. do not write
                               ! write_on_screen = .true.  write
    LOGICAL :: write_for_mlab ! write_for_mlab = .false. do not write
                               ! write_for_mlab = .true.  write
    LOGICAL :: select_one_branch ! to select one branch of hyperbolas
! treshold for writing data
    INTEGER :: threshold
! auxiliary
    CHARACTER*1 :: wrimat,wriscreen,wricompl
    DOUBLE PRECISION :: cosupl,sinupl,cosu,sinu,uupl,uu,chupl,chu
    DOUBLE PRECISION ,DIMENSION(poldeg) :: zeros 
    INTEGER :: nstattmp,anstmp(poldeg)
    DOUBLE PRECISION :: v1tmp(poldeg),v2tmp(poldeg)
! loop indexes
    INTEGER h,i,j,k,l,m,enn,o,p,q,count,ii,nfail,jj
    SAVE
! ==================================================================

!    verb_moid=20
! initialization
    count = 0
    nfail = 0
    zeros(1:poldeg)=0.d0
! default values (dummy)
    cosupl=1.d0
    sinupl=0.d0
    cosu=1.d0
    sinu=0.d0
    uupl=0.d0
    uu=0.d0

! reading orbital data
    OPEN(1,file='orbit_cometary.dat',status='old')
    READ(1,*)
    READ(1,*)perih1,e1,i1,Omn1,om1
    READ(1,*)
    READ(1,*)perih2,e2,i2,Omn2,om2
    CLOSE(1)
! 1st orbit
    i1 = i1*radeg
    Omn1 = Omn1*radeg
    om1 = om1*radeg
! 2nd orbit
    i2 = i2*radeg
    Omn2 = Omn2*radeg
    om2 = om2*radeg

! options file
    OPEN(10,file='stat_search.opt',status='old')
    READ(10,*)
    READ(10,*) ph1step,e1step,i1step,Omn1step,om1step
    READ(10,*)
    READ(10,*)
    READ(10,*) ph1iter,e1iter,i1iter,Omn1iter,om1iter
    READ(10,*)
    READ(10,*)
    READ(10,*) ph2step,e2step,i2step,Omn2step,om2step
    READ(10,*)
    READ(10,*)
    READ(10,*) ph2iter,e2iter,i2iter,Omn2iter,om2iter
    READ(10,*)
    READ(10,*)
    READ(10,*) select_one_branch
    CLOSE(10)
! options requested at the terminal
    WRITE(*,*)'only one iteration? (y/n)'
    READ(*,*)ansiter
    IF(ansiter.eq.'y')THEN
       ph1iter=0
       e1iter=0
       i1iter=0
       Omn1iter=0
       om1iter=0
       ph2iter=0
       e2iter=0
       i2iter=0
       Omn2iter=0
       om2iter=0
    ELSE
       ansiter='n'
    ENDIF
    WRITE(*,*)'write on file for matlab? (y/n)'
    READ(*,*)wrimat
    IF(wrimat.eq.'y')THEN
       write_for_mlab = .true.
    ELSE
       write_for_mlab = .false.
    ENDIF
    WRITE(*,*)'write output on screen? (y/n)'
    READ(*,*)wriscreen
    IF(wriscreen.eq.'y')THEN
       write_on_screen = .true.
    ELSE
       write_on_screen = .false.
    ENDIF
    WRITE(*,*)'threshold for writing output on hnstat.out &
         &(a positive integer):'
    READ(*,*) threshold
    WRITE(*,*)'verbosity? (>20=yes)'
    READ(*,*) verb_moid
    WRITE(*,*)'write complex roots of r(t) on file? (y/n)'
    READ(*,*)wricompl
    IF(wricompl.eq.'y')THEN
       follow_roots = .true.
    ELSE
       follow_roots = .false.
    ENDIF

! opening files
! for FOLLOWROOTS
    IF(follow_roots) THEN
       OPEN(12,file='complex_roots',status='unknown')
    ENDIF
! for MATLAB
    OPEN(11,file='comet_elements',status='unknown')
    OPEN(5,file='mlabstatpts',status='unknown')
! ERRORs
    OPEN(4,file='error_file',status='unknown')
! check with MORSE-WEIERSTRASS
    OPEN(3,file='morse_weier',status='unknown')
! HIGH NUMBER of STATIONARY POINTS
    OPEN(2,file='hnstat.out',status='unknown')
    WRITE(2,*)'threshold=',threshold
    WRITE(2,*)'perih1,e1,i1,Omn1,om1'
    WRITE(2,*)'perih2,e2,i2,Omn2,om2'
    WRITE(2,*)'nstat,nummin,nummax'

! MAIN LOOP
! 1st PERIHELION DISTANCE
    DO h = 0,ph1iter
       WRITE(*,*)'h=', h
       ec1(1)   = perih1 + ph1step*h
! 1st ECCENTRICITY
    DO i = 0,e1iter
       WRITE(*,*)'i=', i
       ec1(2) = e1 + e1step*i
! 1st INCLINATION
    DO j = 0,i1iter
          WRITE(*,*)'j=', j
       ec1(3) = i1 + i1step*j*radeg
! 1st OMNOD
    DO k = 0,Omn1iter    
          WRITE(*,*)'k=', k
       ec1(4) = Omn1 + Omn1step*k*radeg
! 1st OMEGA
    DO l = 0,om1iter
          WRITE(*,*)'l=', l
       ec1(5) = om1 + om1step*l*radeg           
! -------------------------------------------------------              
! 2nd PERIHELION DISTANCE
    DO m = 0,ph2iter
!          WRITE(*,*)'m=', m
       ec2(1)   = perih2 + ph2step*m               
! 2nd ECCENTRICITY
    DO enn = 0,e2iter
!          WRITE(*,*)'enn=', enn
       ec2(2) = e2 + e2step*enn
! 2nd INCLINATION
    DO o = 0,i2iter
!          WRITE(*,*)'o=',o
       ec2(3) = i2 + i2step*o*radeg
! 2nd OMNOD
    DO p = 0,Omn2iter      
!          WRITE(*,*)'p=', p
       ec2(4) = Omn2 + Omn2step*p*radeg
! 2nd OMEGA
    DO q = 0,om2iter
!          WRITE(*,*)'l=', l
       ec2(5) = om2 + om2step*q*radeg
! dummy variables
       ec1(6) = 0.d0
       ec2(6) = 0.d0

       count = count+1     
                     
!           write(*,113)'ec1',ec1(1:5)
!           write(*,113)'ec2',ec2(1:5)
113    FORMAT(a3,1x,5(f10.5,1x))

       circ_copl_orb=.false.
! check for circular coplanar orbits
       IF((ec1(2).eq.0.d0).and. &
            & (ec2(2).eq.0.d0).and. &
            & (ec1(3).eq.0.d0).and. &
            & ec2(3).eq.0.d0) THEN
          write(*,*)'circular coplanar orbits: skip computation'
          write(4,*)'circular coplanar orbits: skip computation'
          write(4,113)'ec1',ec1(1:2),ec1(3:5)*degrad
          write(4,113)'ec2',ec2(1:2),ec2(3:5)*degrad
          nfail = nfail+1
          circ_copl_orb=.true.
          GOTO 13
       ENDIF
       overlap_ell=.false.
! check for overlapping ellipses
       IF((ec1(1).eq.ec2(1)).and. &
&          (ec1(2).eq.ec2(2)).and. &
&          (ec1(3).eq.ec2(3)).and. &
&          (ec1(4).eq.ec2(4)).and. &
&          (ec1(5).eq.ec2(5))) THEN
          write(*,*)'overlapping orbits: skip computation'
          write(4,*)'overlapping orbits: skip computation'
          write(4,113)'ec1',ec1(1:2),ec1(3:5)*degrad
          write(4,113)'ec2',ec2(1:2),ec2(3:5)*degrad
          nfail = nfail+1
          overlap_ell=.true.
          GOTO 13
       ENDIF
          
!          IF(follow_roots) THEN
! write the values of the roots to be discarded
!                         WRITE(12,102) sqrt((e2+1)/(1-e2)),zeros(2:poldeg)
!                         WRITE(*,102) sqrt((e2+1)/(1-e2)),zeros(2:poldeg)
!102          FORMAT(f20.8,15(2x,f20.8))
!          ENDIF

! flags initialization                                              
       weier = .true.
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
       hevalflag = .true.

! ****************************************************
       CALL compute_critical_points_shift(ec1(1:5),&
            & ec2(1:5),v1,v2,nstat,nummin,nummax,&
            & answer,warnflag,sflag,morse,weier,&
            & hzflag,hwflag,multfl,hevalflag)
! ****************************************************
       
! write in error_file
       if(nstat.eq.-1) then
          nfail = nfail+1
          write(*,*)'failed computation:nstat',nstat
          write(*,*)
          write(4,*)'nstat=',nstat
          write(4,113)'ec1',ec1(1:2),ec1(3:5)*degrad
          write(4,113)'ec2',ec2(1:2),ec2(3:5)*degrad
!             stop
          cycle
       endif
       if(.not.hevalflag) then
          write(4,*)'hessian evaluation failed',hevalflag
          write(4,113)'ec1',ec1(1:2),ec1(3:5)*degrad
          write(4,113)'ec2',ec2(1:2),ec2(3:5)*degrad
       endif

       IF(select_one_branch) THEN
! in case of one or two hyperbolas select only one branch
          IF((ec1(2).gt.1).or.(ec2(2).gt.1)) THEN
             nstattmp=0
             DO jj = 1,nstat
                IF((1.d0+ec1(2)*cos(radeg*v1(jj)).gt.0.d0).and. &
                     & (1.d0+ec2(2)*cos(radeg*v2(jj)).gt.0.d0)) THEN
                   ! accept critical point
                   nstattmp = nstattmp+1
                   v1tmp(nstattmp)=v1(jj)
                   v2tmp(nstattmp)=v2(jj)
                   anstmp(nstattmp)=answer(jj)
                ENDIF
             ENDDO
!             write(*,*)'nstat,nstattmp',nstat,nstattmp
             nstat=nstattmp
             nummin=0
             nummax=0
             DO jj=1,nstat
                v1(jj)=v1tmp(jj)
                v2(jj)=v2tmp(jj)
                answer(jj)=anstmp(jj)
!                write(*,*)'answer',answer(jj)
                IF(answer(jj).eq.-1) THEN
                   nummin = nummin+1
!                   write(*,*)'nummin',nummin
                ELSEIF(answer(jj).eq.1) THEN
                   nummax = nummax+1
                ENDIF
             ENDDO
          ENDIF
       ENDIF

! write in files for matlab
       IF (write_for_mlab) THEN
          write(11,114) ec1(1:2),ec1(3:5)*degrad
          write(11,114) ec2(1:2),ec2(3:5)*degrad
114       FORMAT(5(f10.5,1x))
          write(5,*) nstat,nummin,nummax
          DO jj = 1,nstat 
             write(5,104) v1(jj),v2(jj),answer(jj)
          ENDDO
104       FORMAT(f14.7,2x,f14.7,2x,i2)
       ENDIF
       
! compute d^2 at critical (fpl,fcom)
       DO jj = 1,nstat 
          vv1=v1(jj)*radeg ! conversion in radians for D2eval 
          vv2=v2(jj)*radeg 
!             write(*,*)'v1,v2 in radians',vv1,vv2
          CALL d2eval_ta(vv1,vv2,DD2) 
          D2(jj)=DD2 
! angles between 0 and 360 degrees                                  
!                                     CALL choosedeg(v1(jj),chv1) 
!                                     CALL choosedeg(v2(jj),chv2)        
!                                     v1(jj) = chv1 
!                                     v2(jj) = chv2 
       ENDDO

! sorting 3-uples (v1,v2,D2) according to value of D2               
       CALL heapsort(D2,nstat,D2_ord)

! ***************************************************************************
       if (write_on_screen) then
          WRITE(*,*)'#########################################&
               &########################################'
          WRITE(*,*)'######### S T A T I O N A R Y   P O I N T S &
               &#####################################'
          WRITE(*,*)'#########################################&
               &########################################'
          WRITE(*,*)'       u              upl         v2 &
               &           v1             DIST        TYPE'
          WRITE(*,*)'========================================&
               &========================================='
          DO jj = 1,nstat       
             IF(ec1(2).lt.1.d0) THEN
                cosupl = (ec1(2) + cos(radeg*v1(D2_ord(jj))))/ &
                     & (1.d0+ec1(2)*cos(radeg*v1(D2_ord(jj))))
                sinupl = sin(radeg*v1(D2_ord(jj)))*sqrt(1.d0-ec1(2)**2)/&
                     &(1.d0+ec1(2)*cos(radeg*v1(D2_ord(jj))))
                ! write(*,*)'cosupl, sinupl',cosupl,sinupl
                uupl = datan2(sinupl,cosupl)
                CALL choosedeg(degrad*uupl,chupl) 
                uupl=chupl
             ENDIF
             IF(ec2(2).lt.1.d0) THEN
                cosu = (ec2(2) + cos(radeg*v2(D2_ord(jj))))/ &
                     & (1.d0+ec2(2)*cos(radeg*v2(D2_ord(jj))))
                sinu = sin(radeg*v2(D2_ord(jj)))*sqrt(1.d0-ec2(2)**2)/&
                     &(1.d0+ec2(2)*cos(radeg*v2(D2_ord(jj))))
                ! write(*,*)'cosu, sinu',cosu,sinu
                uu = datan2(sinu,cosu)
                CALL choosedeg(degrad*uu,chu) 
                uu=chu
             ENDIF
             WRITE(*,108)uu,uupl,v2(D2_ord(jj)),v1(D2_ord(jj)), &
                  & dsqrt(D2(D2_ord(jj))),answer(D2_ord(jj))
          ENDDO
108       FORMAT(4(2x,f13.7),2x,f12.8,2x,i2) 
       ELSE
       ENDIF
                                  
! if nstat >= threshold then write on file                       
       IF(nstat.ge.threshold) THEN
          write(*,*)'nstat,nummin,nummax = ',nstat,&
               & nummin,nummax 
          WRITE(*,100)ec1(1),ec1(2),ec1(3)*degrad,&
               & ec1(4)*degrad,ec1(5)*degrad
          WRITE(*,100)ec2(1),ec2(2),ec2(3)*degrad,&
               & ec2(4)*degrad,ec2(5)*degrad
          WRITE(2,100)ec1(1),ec1(2),ec1(3)*degrad,&
               & ec1(4)*degrad,ec1(5)*degrad
          WRITE(2,100)ec2(1),ec2(2),ec2(3)*degrad,&
               & ec2(4)*degrad,ec2(5)*degrad
          WRITE(2,101)nstat,nummin,nummax
          WRITE(2,*) 'hevalflag=',hevalflag
          DO ii=1,nstat
             WRITE(2,*)v1(D2_ord(ii)),v2(D2_ord(ii)),D2(D2_ord(ii)), &
                  & answer(D2_ord(ii))
          ENDDO
          WRITE(2,*)'NUMBER',count
          WRITE(2,*)'*************************************'
100       FORMAT(1x,5(f10.3,2x))
101       FORMAT(i3,1x,i3,1x,i3)
       ENDIF
       
! ============= CHECK WITH MORSE THEORY ===============
       IF (.not.morse) THEN
          if(verb_moid.ge.20)then
             WRITE(*,*)'ERROR: MORSE THEORY VIOLATION!!!'
          endif
          WRITE(3,*)'ERROR: MORSE THEORY VIOLATION!!!'   
          WRITE(3,*)'FAILED COMPUTATION FOR DATA'
          WRITE(3,100)ec1(1),ec1(2),ec1(3)*degrad,&
               & ec1(4)*degrad,ec1(5)*degrad
          WRITE(3,100)ec2(1),ec2(2),ec2(3)*degrad,&
               & ec2(4)*degrad,ec2(5)*degrad
          WRITE(3,*)'NUMSTAT=',nstat
          WRITE(3,*)'NUMMIN=',nummin
          WRITE(3,*)'NUMMAX=',nummax
          nfail = nfail+1
       ENDIF
! ============ CHECK WITH WEIERSTRASS THEOREM =========
       IF(.not.weier)THEN
          if(verb_moid.ge.20)then
             WRITE(*,*)'ERROR: WEIERSTRASS THEOREM VIOLATION!!!'
          endif
          WRITE(3,*)'ERROR: WEIERSTRASS THEOREM VIOLATION!!!'
          WRITE(3,*)'FAILED COMPUTATION FOR DATA'
          WRITE(3,100)ec1(1),ec1(2),ec1(3)*degrad,&
               & ec1(4)*degrad,ec1(5)*degrad
          WRITE(3,100)ec2(1),ec2(2),ec2(3)*degrad,&
               & ec2(4)*degrad,ec2(5)*degrad
          IF(morse)THEN !add failure only if 
             !morse check is ok
             nfail = nfail + 1
          ENDIF
       ENDIF
       
13     CONTINUE ! to skip the computation
          
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

! closing files
    CLOSE(2) ! hnstat.out
    CLOSE(3) ! morse_weier
    CLOSE(4) ! error_file
    CLOSE(5) ! matlab files
    CLOSE(11)
    IF(follow_roots) THEN
       CLOSE(12) ! follow_roots file
    ENDIF

    write(*,*)'Number of processed orbits:',count
    write(*,*)'Number of failed computations:',nfail
    
  END PROGRAM stat_search_ta_shift

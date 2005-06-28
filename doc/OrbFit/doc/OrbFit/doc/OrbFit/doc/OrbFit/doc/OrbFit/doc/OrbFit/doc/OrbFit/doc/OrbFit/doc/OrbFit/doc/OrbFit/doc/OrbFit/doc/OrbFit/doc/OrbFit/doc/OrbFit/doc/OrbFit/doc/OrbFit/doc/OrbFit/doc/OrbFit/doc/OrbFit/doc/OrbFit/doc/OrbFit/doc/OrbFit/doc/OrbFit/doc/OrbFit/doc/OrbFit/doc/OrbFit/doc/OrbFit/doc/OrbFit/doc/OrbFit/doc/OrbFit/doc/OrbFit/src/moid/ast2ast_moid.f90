! ******************************************************
! *******  PROGRAM   A S T 2 A S T _ M O I D ***********
! ****** written by GIOVANNI F. GRONCHI (Mar.2005) *****
! ******************************************************
  PROGRAM ast2ast_moid
    USE fund_const
    USE output_control
    USE critical_points_shift
    USE orbit_elements
!    USE ouput_control
    IMPLICIT NONE
    INTEGER, PARAMETER :: norbx=1000000
! for read_elems
    REAL(KIND=8) :: sax
    CHARACTER*80 :: catnam
    INTEGER :: iunin
    TYPE(orbit_elem) :: elcom(norbx) !cometary elements
    CHARACTER*9 :: name(norbx)
    LOGICAL :: eof
! for heapsort
    INTEGER :: hmagx
    INTEGER :: norb !number of orbits in the catalog
    INTEGER,DIMENSION(norbx) :: srtnum
! for coo_cha
    INTEGER :: fail_flag
!
!    INTEGER, PARAMETER :: poldeg=16
! cometary elements of the first and second orbit, temporary keplerian
    TYPE(orbit_elem) :: ec1,ec2,elkep
! true anomalies
    REAL(KIND=8) :: v1(poldeg),v2(poldeg),vv1,vv2,chv1,chv2
    REAL(KIND=8) :: D2(poldeg),DD2 !squared distance function
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
    LOGICAL,PARAMETER :: select_one_branch=.false. ! to select one branch of 
                                                   ! hyperbolic orbits
! treshold for writing data
   INTEGER :: threshold
   REAL(KIND=8) :: eps_moid
! auxiliary
    CHARACTER*1 :: wriscreen
    REAL(KIND=8) :: cosupl,sinupl,cosu,sinu,uupl,uu,chupl,chu
    REAL(KIND=8) ,DIMENSION(poldeg) :: zeros 
    INTEGER :: nstattmp,anstmp(poldeg)
    REAL(KIND=8) :: v1tmp(poldeg),v2tmp(poldeg)
! loop indexes
    INTEGER h,k,ii,jj,count,nfail,countbig
    SAVE
! ==================================================================

! initialization
    threshold=12
    verb_moid=0
    count = 0
    countbig=0
    nfail = 0
    zeros(1:poldeg)=0.d0
! default values (dummy)
    cosupl=1.d0
    sinupl=0.d0
    cosu=1.d0
    sinu=0.d0
    uupl=0.d0
    uu=0.d0

! option
    WRITE(*,*)'maximal semimajor axis'
    READ(*,*) sax
! reading orbital data
    catnam='allnum.cat'
    CALL filopn(iunin,catnam,'OLD')
    CALL oporbf(catnam,iunin)
    norb=0
    DO h=1,norbx
       name(h)='h'
       CALL read_elems(elkep,name(h),eof,'allnum.cat',iunin)
       IF(eof) EXIT
! 1st orbit
       IF(elkep%coord(1).gt.sax)CYCLE
       norb=norb+1
       CALL coo_cha(elkep,'COM',elcom(norb),fail_flag)!conversion into cometary
!       write(*,*)'elements:',elkep(h)%coord(1:6)
    ENDDO
    WRITE(*,*)' number of orbits in input ', norb
    CALL clorbf
! sorting by absolute magnitude 
    srtnum(1:norbx)=0 !initialization
    CALL heapsort(elcom(1:norb)%h_mag,norb,srtnum)
! check
!    write(*,*) 'elements sorted by abs. magnitude'
!    DO h =1,10
!       write(*,*) elcom(srtnum(h))%h_mag
!    ENDDO

! other options requasted at the terminal
    WRITE(*,*)'maximal absolute magnitude for more massive asteroids:'
    READ(*,*) hmagx
    WRITE(*,*)'write output on screen? (y/n)'
    READ(*,*)wriscreen
    IF(wriscreen.eq.'y')THEN
       write_on_screen = .true.
    ELSE
       write_on_screen = .false.
    ENDIF
!    WRITE(*,*)'nstatmin (threshold to write in hnstat.out): &
!         &(a positive integer):'
!    READ(*,*) threshold
    WRITE(*,*)'eps_moid max (threshold to write in low_moid.out):'
    READ(*,*) eps_moid
    WRITE(*,*)'verbosity? (>20=yes)'
    READ(*,*) verb_moid

! output files
! ERRORs
    OPEN(4,file='error_file',status='unknown')
! check with MORSE-WEIERSTRASS
    OPEN(3,file='morse_weier',status='unknown')
! HIGH NUMBER of STATIONARY POINTS
    OPEN(2,file='hnstat.out',status='unknown')
    WRITE(2,*)'perih1,e1,i1,Omn1,om1'
    WRITE(2,*)'perih2,e2,i2,Omn2,om2'
    WRITE(2,*)'nstat,nummin,nummax'
    OPEN(1,file='low_moid.out',status='unknown')
    WRITE(1,*) 'eps_moid=',eps_moid
    WRITE(1,*) 'number 1    H1      number 2    H2          MOID'

113 FORMAT(a3,1x,5(f10.5,1x))

    DO h =1,norb-1
       IF(elcom(srtnum(h))%h_mag.gt.hmagx) THEN
          EXIT
       ENDIF
       countbig=countbig+1
       DO k = h+1,norb
          IF(elcom(srtnum(k))%h_mag.gt.hmagx) THEN
             EXIT
          ENDIF
          count = count+1
! 1st orbit
          ec1=elcom(srtnum(h))
! 2nd orbit
          ec2=elcom(srtnum(k))
 
          circ_copl_orb=.false.
! check for circular coplanar orbits
          IF((ec1%coord(2).eq.0.d0).and. &
               & (ec2%coord(2).eq.0.d0).and. &
               & (ec1%coord(3).eq.0.d0).and. &
               & ec2%coord(3).eq.0.d0) THEN
             write(*,*)'circular coplanar orbits: skip computation'
             write(4,*)'circular coplanar orbits: skip computation'
             write(4,113)'ec1',ec1%coord(1:2),ec1%coord(3:5)*degrad
             write(4,113)'ec2',ec2%coord(1:2),ec2%coord(3:5)*degrad
             nfail = nfail+1
             circ_copl_orb=.true.
             GOTO 13
          ENDIF
          overlap_ell=.false.
! check for overlapping ellipses
          IF((ec1%coord(1).eq.ec2%coord(1)).and. &
               &          (ec1%coord(2).eq.ec2%coord(2)).and. &
               &          (ec1%coord(3).eq.ec2%coord(3)).and. &
               &          (ec1%coord(4).eq.ec2%coord(4)).and. &
               &          (ec1%coord(5).eq.ec2%coord(5))) THEN
             write(*,*)'overlapping orbits: skip computation'
             write(4,*)'overlapping orbits: skip computation'
             write(4,113)'ec1',ec1%coord(1:2),ec1%coord(3:5)*degrad
             write(4,113)'ec2',ec2%coord(1:2),ec2%coord(3:5)*degrad
             nfail = nfail+1
             overlap_ell=.true.
             GOTO 13
          ENDIF
          
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

       ! *********************************************
          CALL compute_critical_points_shift(ec1%coord(1:5),&
               & ec2%coord(1:5),v1,v2,nstat,nummin,nummax,&
               & answer,warnflag,sflag,morse,weier,&
               & hzflag,hwflag,multfl,hevalflag)
       ! *********************************************
          
! write in error_file
          if(nstat.eq.-1) then
             nfail = nfail+1
             write(*,*)'failed computation:nstat',nstat
             write(*,*)
             write(4,*)'nstat=',nstat
             write(4,113)'ec1',ec1%coord(1:2),ec1%coord(3:5)*degrad
             write(4,113)'ec2',ec2%coord(1:2),ec2%coord(3:5)*degrad
!             stop
             cycle
          endif
          if(.not.hevalflag) then
             write(4,*)'hessian evaluation failed',hevalflag
             write(4,113)'ec1',ec1%coord(1:2),ec1%coord(3:5)*degrad
             write(4,113)'ec2',ec2%coord(1:2),ec2%coord(3:5)*degrad
          endif

          IF(select_one_branch) THEN
! in case of one or two hyperbolas select only one branch
             IF((ec1%coord(2).gt.1).or.(ec2%coord(2).gt.1)) THEN
                nstattmp=0
                DO jj = 1,nstat
                   IF((1.d0+ec1%coord(2)*cos(radeg*v1(jj)).gt.0.d0).and. &
                        & (1.d0+ec2%coord(2)*cos(radeg*v2(jj)).gt.0.d0)) THEN
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
                   ! write(*,*)'answer',answer(jj)
                   IF(answer(jj).eq.-1) THEN
                      nummin = nummin+1
!                   write(*,*)'nummin',nummin
                   ELSEIF(answer(jj).eq.1) THEN
                      nummax = nummax+1
                   ENDIF
                ENDDO
             ENDIF
          ENDIF

! compute d^2 at critical (fpl,fcom)
          DO jj = 1,nstat 
             vv1=v1(jj)*radeg! conversion in radians for D2eval 
             vv2=v2(jj)*radeg 
             !write(*,*)'v1,v2 in radians',vv1,vv2
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
             WRITE(*,*)'###### S T A T I O N A R Y   P O I N T S &
                  & for PAIR ',name(srtnum(h)),name(srtnum(k)),'#####'
             WRITE(*,*)'#########################################&
                  &########################################'
             WRITE(*,*)'       u              upl         v2 &
                  &           v1             DIST        TYPE'
             WRITE(*,*)'========================================&
                  &========================================='
             DO jj = 1,nstat       
                IF(ec1%coord(2).lt.1.d0) THEN
                   cosupl = (ec1%coord(2) + cos(radeg*v1(D2_ord(jj))))/ &
                        & (1.d0+ec1%coord(2)*cos(radeg*v1(D2_ord(jj))))
                   sinupl = sin(radeg*v1(D2_ord(jj)))* &
                        & sqrt(1.d0-ec1%coord(2)**2)/&
                        &(1.d0+ec1%coord(2)*cos(radeg*v1(D2_ord(jj))))
                   ! write(*,*)'cosupl, sinupl',cosupl,sinupl
                   uupl = datan2(sinupl,cosupl)
                   CALL choosedeg(degrad*uupl,chupl) 
                   uupl=chupl
                   ! uupl=uupl*degrad
                ENDIF
                IF(ec2%coord(2).lt.1.d0) THEN
                   cosu = (ec2%coord(2) + cos(radeg*v2(D2_ord(jj))))/ &
                        & (1.d0+ec2%coord(2)*cos(radeg*v2(D2_ord(jj))))
                   sinu = sin(radeg*v2(D2_ord(jj)))* &
                        & sqrt(1.d0-ec2%coord(2)**2)/&
                        &(1.d0+ec2%coord(2)*cos(radeg*v2(D2_ord(jj))))
                   ! write(*,*)'cosu, sinu',cosu,sinu
                   uu = datan2(sinu,cosu)
                   CALL choosedeg(degrad*uu,chu) 
                   uu=chu
                   ! uu=uu*degrad
                ENDIF
                WRITE(*,108)uu,uupl,v2(D2_ord(jj)),v1(D2_ord(jj)), &
                     & dsqrt(D2(D2_ord(jj))),answer(D2_ord(jj))
             ENDDO
108          FORMAT(4(2x,f13.7),2x,f12.8,2x,i2) 
          ELSE
          ENDIF
          
          ! if nstat >= threshold then write on file                       
          IF(nstat.ge.threshold) THEN
             write(*,*)'nstat,nummin,nummax = ',nstat,&
                  & nummin,nummax 
             WRITE(*,100)ec1%coord(1),ec1%coord(2),ec1%coord(3)*degrad,&
                  & ec1%coord(4)*degrad,ec1%coord(5)*degrad
             WRITE(*,100)ec2%coord(1),ec2%coord(2),ec2%coord(3)*degrad,&
                  & ec2%coord(4)*degrad,ec2%coord(5)*degrad
             WRITE(2,100)ec1%coord(1),ec1%coord(2),ec1%coord(3)*degrad,&
                  & ec1%coord(4)*degrad,ec1%coord(5)*degrad
             WRITE(2,100)ec2%coord(1),ec2%coord(2),ec2%coord(3)*degrad,&
                  & ec2%coord(4)*degrad,ec2%coord(5)*degrad
             WRITE(2,101)nstat,nummin,nummax
             WRITE(2,*) 'hevalflag=',hevalflag
             DO ii=1,nstat
                WRITE(2,*)v1(D2_ord(ii)),v2(D2_ord(ii)),D2(D2_ord(ii)), &
                     & answer(D2_ord(ii))
             ENDDO
             WRITE(2,*)'NUMBER',count
             WRITE(2,*)'*************************************'
100          FORMAT(1x,5(f10.3,2x))
101          FORMAT(i3,1x,i3,1x,i3)
          ENDIF
          
          ! if moid <= eps_moid then write on file                       
          IF(sqrt(D2(D2_ord(1))).le.eps_moid) THEN
             WRITE(1,114) name(srtnum(h)),'&',elcom(srtnum(h))%h_mag, &
                  & '&',name(srtnum(k)),'&',elcom(srtnum(k))%h_mag, &
                  &'&','&',sqrt(D2(D2_ord(1)))
114          FORMAT(a7,2x,a1,f6.3,2x,a1,a7,2x,a1,f6.3,4x,a1,2x,a1,f10.7)
!             WRITE(1,100)ec1%coord(1),ec1%coord(2),ec1%coord(3)*degrad,&
!                  & ec1%coord(4)*degrad,ec1%coord(5)*degrad
!             WRITE(1,100)ec2%coord(1),ec2%coord(2),ec2%coord(3)*degrad,&
!                  & ec2%coord(4)*degrad,ec2%coord(5)*degrad
!             WRITE(1,101)nstat,nummin,nummax
!             WRITE(1,*) 'hevalflag=',hevalflag
!             DO ii=1,nstat
!                WRITE(1,*)v1(D2_ord(ii)),v2(D2_ord(ii)),D2(D2_ord(ii)), &
!                     & answer(D2_ord(ii))
!             ENDDO
!             WRITE(1,*)'NUMBER',count
!             WRITE(1,*)'*************************************'
          ENDIF
          
          ! ============= CHECK WITH MORSE THEORY ===============
          IF (.not.morse) THEN
             if(verb_moid.ge.20)then
                WRITE(*,*)'ERROR: MORSE THEORY VIOLATION!!!'
             endif
             WRITE(3,*)'ERROR: MORSE THEORY VIOLATION!!!'   
             WRITE(3,*)'FAILED COMPUTATION FOR DATA'
             WRITE(3,100)ec1%coord(1),ec1%coord(2),ec1%coord(3)*degrad,&
                  & ec1%coord(4)*degrad,ec1%coord(5)*degrad
             WRITE(3,100)ec2%coord(1),ec2%coord(2),ec2%coord(3)*degrad,&
                  & ec2%coord(4)*degrad,ec2%coord(5)*degrad
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
             WRITE(3,100)ec1%coord(1),ec1%coord(2),ec1%coord(3)*degrad,&
                  & ec1%coord(4)*degrad,ec1%coord(5)*degrad
             WRITE(3,100)ec2%coord(1),ec2%coord(2),ec2%coord(3)*degrad,&
                  & ec2%coord(4)*degrad,ec2%coord(5)*degrad
             IF(morse)THEN !add failure only if 
                !morse check is ok
                nfail = nfail + 1
             ENDIF
          ENDIF
          
13        CONTINUE ! to skip the computation
       ENDDO
    ENDDO
    
    ! closing files
    CLOSE(1) ! low_moid.out
    CLOSE(2) ! hnstat.out
    CLOSE(3) ! morse_weier
    CLOSE(4) ! error_file
    
    write(*,*)'Number of processed couples of orbits:',count
    write(*,*)'Number of big asteroids', countbig
    write(*,*)'Number of failed computations:',nfail
    
  END PROGRAM ast2ast_moid

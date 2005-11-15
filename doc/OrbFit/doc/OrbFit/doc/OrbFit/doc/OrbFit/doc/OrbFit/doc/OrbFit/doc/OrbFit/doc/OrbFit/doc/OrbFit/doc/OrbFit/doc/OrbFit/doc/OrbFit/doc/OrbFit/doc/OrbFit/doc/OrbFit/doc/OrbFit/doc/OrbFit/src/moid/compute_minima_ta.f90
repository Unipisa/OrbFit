
! ************************************************************      
! ***** C O M P U T I N G    M I N I M U M   P O I N T S *****      
! ************************************************************      
! ******* written by GIOVANNI F. GRONCHI (April 2005) ********
! ************************************************************      
! last modified 24/06/2005 GFG
! ==================================================================      
  SUBROUTINE compute_minima_ta(car,carpl,iplam,cmin,cplmin,D2min,nummin)
    USE critical_points_shift
    USE fund_const 
    USE orbit_elements
    USE planet_masses
    IMPLICIT NONE 
! cartesian coordinates of the asteroid and of the planet:
    DOUBLE PRECISION,INTENT(IN) :: car(6)
    DOUBLE PRECISION,INTENT(IN) :: carpl(6)
    INTEGER,INTENT(IN) :: iplam ! planet number (3 = Earth) 
! cartesian coordinates of the asteroid and of the planet 
! at minimum points:
!    INTEGER, PARAMETER :: poldeg=16 !polynomial degree
! cartesian position of the critical points on the two orbits
    DOUBLE PRECISION,INTENT(OUT) :: cmin(3,poldeg),cplmin(3,poldeg) 
    DOUBLE PRECISION,INTENT(OUT) :: D2min(poldeg) ! squared distance function 
                                                  ! at minimum points 
    INTEGER,INTENT(OUT) :: nummin ! number of relative minima found 
! ========= end interface ==========================================
    TYPE(orbit_elem) :: ecarpl ! cartesian elements of the planet
    TYPE(orbit_elem) :: ecar   ! cartesian elements of the asteroid/comet
! ---------- for compute_critical_points_shift ---------------------
    TYPE(orbit_elem) :: ecpl   ! cometary elements of the planet
    TYPE(orbit_elem) :: ec     ! cometary elements of the asteroid/comet
    REAL(KIND=8):: fpl(poldeg),fcom(poldeg)! true anomalies  
    INTEGER :: nstat,nummax ! number of statpts,minima and maxima
! error flags
    INTEGER :: answer(poldeg) ! type of critical point:
!                                         answer =  1      maximum
!                                         answer =  0      saddle
!                                         answer = -1      minimum     
!                                         answer = -2      cannot decide
    LOGICAL :: morse ! check with Morse theory 
!                        morse = true      OK   
!                        morse = false     ERROR  
    LOGICAL :: weier ! check with Weierstrass theory:
!                        weier = true      OK 
!                        weier = false     ERROR
    LOGICAL :: warnflag(3) ! program warning flags:
! warnflag(1) = true    OK
! warnflag(1) = false   leading coefficient of \tilde{r}(t) is very small
! warnflag(2) = true    OK                                          
! warnflag(2) = false   higher degree terms \tilde{r}(t) are not small
! warnflag(3) = true    OK                                           
! warnflag(3) = false   low precision for \tilde{r}(t) coefficients     
    LOGICAL :: sflag(6) ! solving system messages:
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
    LOGICAL :: hzflag ! hzflag = .true.  OK!
                      ! hzflag = .false. abs(root)>10^5
    LOGICAL :: hwflag ! hwflag = .true.  OK!
                      ! hwflag = .false. abs(root)>10^5
    LOGICAL :: multfl ! multfl = .true.  OK!
                      ! multfl = .false. 0 has multiplicity > 4
    LOGICAL :: hevalflag ! hevalflag = .true. OK!
                         ! hevalflag = .false. unsuccessful 
                         !                     hessian evaluation
! -------------------------------------------------------------------
!   auxiliary
    DOUBLE PRECISION :: cosn,sinn,coso,sino,cosi,sini
!   pos and vel unit vectors at pericenter
    DOUBLE PRECISION,DIMENSION(3) :: x0_ast,y0_ast,x0_pl,y0_pl
    DOUBLE PRECISION,DIMENSION(3) :: hatC_ast,hatC_pl ! ang.mom. unit vector
    DOUBLE PRECISION,DIMENSION(3,3) :: inv_chb_ast,inv_chb_pl ! for basis change
    DOUBLE PRECISION,DIMENSION(3,3) :: chb_ast,chb_pl ! for basis change
    DOUBLE PRECISION,DIMENSION(3,3) :: rtmp_ast,rtmp_pl ! auxiliary matrixes
    DOUBLE PRECISION,DIMENSION(3,3) :: rotm_ast,rotm_pl ! rotation matrixes
    DOUBLE PRECISION,DIMENSION(3,poldeg) :: hatP_ast,hatP_pl !versors to crit.pts
    DOUBLE PRECISION,DIMENSION(poldeg) :: r_ast,r_pl ! distance at crit.pts
    DOUBLE PRECISION :: ffcom,ffpl,chfcom,chfpl 
    DOUBLE PRECISION :: D2(poldeg) ! squared distance function 
    DOUBLE PRECISION :: DD2 
    DOUBLE PRECISION :: fplmin(poldeg),fmin(poldeg)
! -------------------------------------------------
    INTEGER,DIMENSION(poldeg) :: srtnum !ordered 3-uples, for sorting 
    INTEGER :: fail_flag ! for coo_cha
    INTEGER :: k,j,kmin ! loop index 
! ==================================================================

!   initialization
    cmin = 0.d0
    cplmin = 0.d0
    D2min(1:poldeg) = 0.d0
    srtnum(1:poldeg)=0 

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

! conversion into COMETARY elements
    ecar=undefined_orbit_elem
    ecar%coord(1:6)=car(1:6)
    ecar%coo='CAR'
    ecar%t=0.d0 !dummy time
    CALL coo_cha(ecar,'COM',ec,fail_flag) !conversion into cometary elems
!    write(*,*)'compute_minima; ast elements(com):',ec%coord(1:6)
    IF(fail_flag.gt.5)THEN !no conversion provided
       write(*,*)'compute_minima_ta: failed conversion into cometary elements &
            & for the asteroid'
       nummin=0
       RETURN              
    ENDIF

    ecarpl=undefined_orbit_elem
    ecarpl%coord(1:6)=carpl(1:6)
    ecarpl%coo='CAR'
    ecarpl%t=0.d0 !dummy time
    CALL coo_cha(ecarpl,'COM',ecpl,fail_flag) !conversion into cometary elems
!    write(*,*)'compute_minima; pl elements(com):',ecpl%coord(1:6)
    IF(fail_flag.gt.5)THEN !no conversion provided
       write(*,*)'compute_minima_ta: failed conversion into cometary elements&
            & for the planet'
       nummin=0
       RETURN              
    ENDIF

! *********************************************************************
    CALL compute_critical_points_shift(ecpl%coord(1:5),ec%coord(1:5),fpl,&
         & fcom,nstat,nummin,nummax,answer,warnflag,sflag,morse,weier, &
         & hzflag,hwflag,multfl,hevalflag)
! *********************************************************************

! check of the solution
    if(.not.sflag(3)) then
       write(*,*)'negative discriminant!  ecc=',ec%coord(2)
    endif

    IF(nstat.le.0) THEN
       WRITE(*,*)'compute_minima_ta: error! nstat=',nstat
       nummin=0
       RETURN
    ENDIF

! loop on number of stationary points                               
    DO j = 1,nstat 
! conversion in radians for D2eval                                  
       ffpl=fpl(j)*radeg 
       ffcom=fcom(j)*radeg 
! compute d^2 at critical (fpl,fcom)
       CALL d2eval_ta(ffpl,ffcom,DD2) 
       D2(j)=DD2 
! angles between 0 and 360 degrees                                  
       CALL choosedeg(fpl(j),chfpl) 
       CALL choosedeg(fcom(j),chfcom)        
       fpl(j) = chfpl 
       fcom(j) = chfcom 
    ENDDO

! sorting 3-uples (fcom,fpl,D2) according to value of D2               
    CALL heapsort(D2,nstat,srtnum)
    kmin = 0 !counting minimum points
    DO k = 1,nstat 
!   write(*,*)'fpl,fcom,d2',fpl(srtnum(k)),fcom(srtnum(k)),sqrt(D2(srtnum(k)))
       IF(answer(srtnum(k)).eq.-1) THEN
          kmin=kmin+1
          fplmin(kmin) = fpl(srtnum(k))*radeg
          fmin(kmin) = fcom(srtnum(k))*radeg
          D2min(kmin) = D2(srtnum(k))
       ENDIF
    ENDDO
!    write(*,*)'kmin,nummin',kmin,nummin

! compute Cartesian coords of the critical points 
! on the asteroid and planet orbits 
    cosn=cos(ec%coord(4))
    sinn=sin(ec%coord(4))
    coso=cos(ec%coord(5))
    sino=sin(ec%coord(5))
    cosi=cos(ec%coord(3))
    sini=sin(ec%coord(3))
    x0_ast=(/coso*cosn-sino*sinn*cosi, coso*sinn+sino*cosn*cosi, sino*sini/)
    y0_ast=(/-sino*cosn-coso*sinn*cosi, -sino*sinn+coso*cosn*cosi, coso*sini/)
    CALL prvec(x0_ast,y0_ast,hatC_ast)
    inv_chb_ast(1:3,1) = x0_ast(1:3)
    inv_chb_ast(1:3,2) = y0_ast(1:3)
    inv_chb_ast(1:3,3) = hatC_ast(1:3) 
    chb_ast = TRANSPOSE(inv_chb_ast) ! basis change matrix

    cosn=cos(ecpl%coord(4))
    sinn=sin(ecpl%coord(4))
    coso=cos(ecpl%coord(5))
    sino=sin(ecpl%coord(5))
    cosi=cos(ecpl%coord(3))
    sini=sin(ecpl%coord(3))
    x0_pl=(/coso*cosn-sino*sinn*cosi, coso*sinn+sino*cosn*cosi, sino*sini/)
    y0_pl=(/-sino*cosn-coso*sinn*cosi, -sino*sinn+coso*cosn*cosi, coso*sini/)
    CALL prvec(x0_pl,y0_pl,hatC_pl)
    inv_chb_pl(1:3,1) = x0_pl(1:3)
    inv_chb_pl(1:3,2) = y0_pl(1:3)
    inv_chb_pl(1:3,3) = hatC_pl(1:3) 
    chb_pl = TRANSPOSE(inv_chb_pl) ! basis change matrix

    DO k=1,kmin
       CALL rotmt(-fmin(k),rtmp_ast,3)
       CALL rotmt(-fplmin(k),rtmp_pl,3)
       rotm_ast = MATMUL(inv_chb_ast,MATMUL(rtmp_ast,chb_ast))
       rotm_pl = MATMUL(inv_chb_pl,MATMUL(rtmp_pl,chb_pl))
       hatP_ast(1:3,k) = MATMUL(rotm_ast,X0_ast)
       hatP_pl(1:3,k) = MATMUL(rotm_pl,X0_pl)
       r_ast(k) = ec%coord(1)*(1.d0+ec%coord(2))/(1.d0+ec%coord(2)*cos(fmin(k)))
       cmin(1:3,k) = r_ast(k)*hatP_ast(1:3,k) 
       r_pl = ecpl%coord(1)*(1.d0+ecpl%coord(2))/(1.d0+ecpl%coord(2)* &
            & cos(fplmin(k)))
       cplmin(1:3,k) = r_pl(k)*hatP_pl(1:3,k) 
    ENDDO

  END SUBROUTINE compute_minima_ta
  

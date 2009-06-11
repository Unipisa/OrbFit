! *** M O D U L E  kepint ***
! *** written by  G.F. Gronchi and L. Dimare ***
! *** December 2008 ***
MODULE kepint
  USE fund_const
  USE force_sat, ONLY: eamoon_mass
  USE planet_masses, ONLY: gmearth
  PRIVATE
  INTEGER, PARAMETER, PUBLIC  :: nprelx=48 ! maximum number of prelim. orbits
  INTEGER,PARAMETER :: Nev=64,expo=6,poldeg=48
  INTEGER :: obscode1,obscode2
! heliocentric pos/vel of the observer
  REAL(KIND=dkind),DIMENSION(3,2) :: pe,pde 

  PUBLIC twobodyint

  REAL(KIND=dkind),DIMENSION(3,2) :: Dvec,Evec,Fvec,Gvec !ang mom coeffs
  REAL(KIND=dkind),DIMENSION(2) :: c0,c1,c2,c3,c31,c32,c4,c5 !energy coeffs
  REAL(KIND=dkind),DIMENSION(2) :: cosa,sina,cosd,sind 
  REAL(KIND=dkind),DIMENSION(3,2) :: rhohat,rhohat_a,rhohat_d 
  REAL(KIND=dkind) :: gmvar

  CONTAINS

SUBROUTINE twobodyint(att1,att2,maxnorm,nprel,orbprel, &
         & qobs1,qpobs1,qobs2,qpobs2,orbunc,deltaoml,norm,nsol_step,checkder)
  USE orbit_elements
  USE attributable
  USE output_control
  USE station_coordinates, ONLY: statcode
  IMPLICIT NONE
  TYPE(attrib),INTENT(IN) :: att1,att2
  REAL(KIND=dkind), INTENT(IN) :: maxnorm
  INTEGER, INTENT(OUT):: nprel ! number of preliminary orbits accepted
  TYPE(orbit_elem),INTENT(OUT) :: orbprel(nprelx)
! OPTIONAL arguments
! interpolated geocentric position and velocity of the observer (equatorial);
  REAL(KIND=dkind), INTENT(IN), DIMENSION(3),OPTIONAL :: qobs1,qpobs1
  REAL(KIND=dkind), INTENT(IN), DIMENSION(3),OPTIONAL :: qobs2,qpobs2
  TYPE(orb_uncert),INTENT(OUT),OPTIONAL  :: orbunc(nprelx) !covariance
  REAL(KIND=dkind),INTENT(OUT),OPTIONAL :: deltaoml(nprelx,2) 
  REAL(KIND=dkind),INTENT(OUT),OPTIONAL :: norm(nprelx) ! identification norm
! num. of solutions passing different controls
!  1) num. real sols of q=p=0 after discarding those with negative components
!     and applying N-R, discarding not-convergent 
!     (hint! after N-R some components could have become negative)
!  2) num. sols with no very small component, after removing duplicates
!  3) num. sols after elim. spurious
!  4) num. sols after elim. hyperbolic orbits
  INTEGER,INTENT(OUT),OPTIONAL :: nsol_step(4) 
  LOGICAL,INTENT(IN),OPTIONAL ::checkder
! ----------- end interface --------------------------------------------
  INTEGER :: nsol,nsol_tmp ! number of positive real roots
  REAL(KIND=dkind) :: idnorm(nprelx) ! identification norm

! Keplerian elements of the asteroid
  TYPE(orbit_elem),DIMENSION(nprelx) :: elkep1,elkep2
  TYPE(orbit_elem),DIMENSION(nprelx) :: elkep1_aux,elkep2_aux
  TYPE(orbit_elem),DIMENSION(nprelx) :: elcar1,elcar2 ! asteroid CAR elems
  TYPE(orbit_elem) :: elatt1,elatt2 ! attributable elements of the asteroid
  TYPE(orbit_elem) :: elcar1_aux, elcar2_aux
  INTEGER :: fail_flag

  REAL(KIND=dkind) :: n1(nprelx),n1_h ! mean motion

  REAL(KIND=dkind), DIMENSION(nprelx,2) :: rrdot1, rrdot2 !radial dist/vel
  REAL(KIND=dkind), DIMENSION(nprelx,2) :: rrdot1_aux, rrdot2_aux
  REAL(KIND=dkind), DIMENSION(nprelx,2) :: rrdot1_h, rrdot2_h 
!  REAL(KIND=qkind), DIMENSION(2) :: rrd1, rrd2 !radial dist/vel (for i=isel)
!  REAL(KIND=qkind), DIMENSION(2) :: rrd1_h, rrd2_h !auxiliary
  REAL(KIND=dkind), DIMENSION(3,2) :: qqd1,qqd2 !obs pos/vel
  REAL(KIND=dkind), DIMENSION(4,4) :: der_phi_r
!derivatives w.r.t. att
  REAL(KIND=dkind), DIMENSION(4,8) :: der_phi_att
  REAL(KIND=dkind), DIMENSION(4,8) :: der_rrd_att 
  REAL(KIND=dkind), DIMENSION(4) :: phi_atth
  REAL(KIND=dkind), DIMENSION(nprelx,4) :: phi

  REAL(KIND=dkind), DIMENSION(6,6) :: der1_att, der2_att
  REAL(KIND=dkind), DIMENSION(poldeg,6,6) :: der1_car , der2_car
  REAL(KIND=dkind), DIMENSION(nprelx,6,6) :: dKEP_dATT_1
  REAL(KIND=dkind), DIMENSION(6,6) :: der_car_att1, der_car_att2 
  REAL(KIND=dkind), DIMENSION(4) :: deltaw_att1, deltaw_att2
  REAL(KIND=dkind), DIMENSION(4) :: deltaell_att1, deltaell_att2
!  REAL(KIND=qkind), DIMENSION(1,2) :: delta_well
  REAL(KIND=dkind), DIMENSION(2,8) :: deltawl_att
  REAL(KIND=dkind) :: deltaoml_aux(nprelx,2) 
! for checkder
  REAL(KIND=dkind), DIMENSION(nprelx,4,4) :: der_phi_r_isel
  REAL(KIND=dkind), DIMENSION(nprelx,4,8) :: der_phi_att_isel
  REAL(KIND=dkind), DIMENSION(nprelx,4,8) :: der_rrd_att_isel
  REAL(KIND=dkind), DIMENSION(nprelx,2,8) :: deltawl_att_isel
!
  REAL(KIND=dkind), DIMENSION(8,2) :: mat2x8
  REAL(KIND=dkind), DIMENSION(8,2) :: tdeltawl_att
  REAL(KIND=dkind), DIMENSION(8,8) :: gamma_att
  REAL(KIND=dkind), DIMENSION(2,2) :: gamma_delta, c_delta

  TYPE(orb_uncert), DIMENSION(nprelx) :: orbunc_aux,orbunc_kep 
  REAL(KIND=dkind), DIMENSION(2,2) :: gamma_R1
  REAL(KIND=dkind), DIMENSION(4,2) :: gamma_a1R1
! for tchinv subroutine
  INTEGER :: ierr
  REAL(KIND=dkind) :: ws2(2),ws6(6) 

  REAL(KIND=dkind) :: tmpsum(2)
  REAL(KIND=dkind) :: gammatmp(8,2)
  REAL(KIND=dkind) :: det_phi, det_gamma
  REAL(KIND=dkind) :: incr(8)!increments:ha1,hd1,had1,hdd1,ha2,hd2,had2,hdd2
  REAL(KIND=dkind),DIMENSION(nprelx) :: deltaom,deltaell
  REAL(KIND=dkind) :: deltaom_h,deltaell_h
  REAL(KIND=dkind) :: pridif,primea,ell2_to_t1,ell2_to_t1_h
  TYPE(attrib) :: att1h,att2h
  TYPE(orbit_elem),DIMENSION(nprelx) :: elcar1_h,elcar2_h
  TYPE(orbit_elem) :: elkep1_h , elkep2_h
  INTEGER :: ising,nprel_h

  INTEGER,DIMENSION(poldeg) :: orb_srt
  INTEGER,DIMENSION(poldeg) :: indsrt
  REAL(KIND=dkind) :: orb_d
  REAL(KIND=dkind), DIMENSION(poldeg) :: orb_dist
  INTEGER :: i,j,k,isrt,isel ! loop indexes
  LOGICAL checkd
! ==================================================================

  IF(rhs.eq.1)THEN
     gmvar=gms
  ELSEIF(rhs.eq.2)THEN
     CALL eamoon_mass
     gmvar=gmearth
  ELSE
     WRITE(*,*)'twobodyint: value not allowed for rhs:',rhs
     STOP
  ENDIF

  IF(PRESENT(checkder))THEN
     checkd=checkder
  ELSE
     checkd=.false.
  ENDIF
!
  IF(PRESENT(nsol_step))THEN
     IF(PRESENT(qobs1).and.PRESENT(qpobs1).and. &
          & PRESENT(qobs2).and.PRESENT(qpobs2))THEN
        CALL orbprelim_2b(att1,att2,nsol,elcar1,elcar2,rrdot1,rrdot2, &
             & nsol_step(1:2),qobs1,qpobs1,qobs2,qpobs2)
     ELSE
        CALL orbprelim_2b(att1,att2,nsol,elcar1,elcar2,rrdot1,rrdot2, &
             & nsol_step(1:2))
     ENDIF
     nsol_step(3)=nsol
  ELSE
     IF(PRESENT(qobs1).and.PRESENT(qpobs1).and. &
          & PRESENT(qobs2).and.PRESENT(qpobs2))THEN
        CALL orbprelim_2b(att1,att2,nsol,elcar1,elcar2,rrdot1,rrdot2, &
             & QOBS1=qobs1,QPOBS1=qpobs1,QOBS2=qobs2,QPOBS2=qpobs2)
     ELSE
        CALL orbprelim_2b(att1,att2,nsol,elcar1,elcar2,rrdot1,rrdot2)
     ENDIF
  ENDIF

  IF(nsol.gt.0)THEN
     nsol_tmp=0
     DO isel=1,nsol
!        write(*,*)'center=',elcar1(isel)%center
        CALL coo_cha(elcar1(isel),'KEP',elkep1(isel),fail_flag, &
             & der1_car(isel,1:6,1:6))
        IF(fail_flag.lt.5)THEN
           IF(verb_2body.gt.10)THEN
!             write(*,*)'elkep1%center=',elkep1(isel)%center
              WRITE(*,*)'----------------------------------------------'
              WRITE(*,101)'elkep1:',elkep1(isel)%coord(1:2), &
                   & elkep1(isel)%coord(3:6)*180.d0/pig,elkep1(isel)%t
           ENDIF
           CALL coo_cha(elcar2(isel),'KEP',elkep2(isel),fail_flag,&
                & der2_car(isel,1:6,1:6))
           IF(fail_flag.lt.5)THEN
              IF(verb_2body.gt.10)THEN
                 WRITE(*,101)'elkep2:',elkep2(isel)%coord(1:2), &
                      & elkep2(isel)%coord(3:6)*180.d0/pig,elkep2(isel)%t
              ENDIF
              nsol_tmp=nsol_tmp+1
              elkep1_aux(nsol_tmp)= elkep1(isel)
              elkep2_aux(nsol_tmp)= elkep2(isel)
              rrdot1_aux(nsol_tmp,1:2) = rrdot1(isel,1:2)
              rrdot2_aux(nsol_tmp,1:2) = rrdot2(isel,1:2)
           ELSEIF(fail_flag.eq.5)THEN
              IF(verb_2body.gt.10)THEN
                 write(*,*)'unbounded orbit: ecc=',elkep2(isel)%coord(2)
              ENDIF
           ELSE
              write(*,*)'ERROR! fail_flag=',fail_flag
           ENDIF
        ELSEIF(fail_flag.eq.5)THEN
              IF(verb_2body.gt.10)THEN
                 write(*,*)'fail_flag:',fail_flag
                 write(*,*)'unbounded orbit: ecc=',elkep1(isel)%coord(2)
              ENDIF
        ELSE
           write(*,*)'ERROR! fail_flag=',fail_flag
        ENDIF
     ENDDO
101  FORMAT(a7,1x,6(f10.6,2x),f12.6)
  ELSE
     IF(verb_2body.gt.10) THEN
        write(*,*)'no preliminary orbit: nsol=',nsol
     ENDIF
     nprel=nsol
     RETURN
  ENDIF
  nprel=nsol_tmp
  IF(nprel.eq.0)RETURN

  IF(PRESENT(nsol_step))THEN
     nsol_step(4)=nprel
  ENDIF

!  write(*,*)'nprel=',nprel
  DO isel=1,nprel
     elkep1(isel)=elkep1_aux(isel)
     elkep2(isel)=elkep2_aux(isel)
     rrdot1(isel,1:2)=rrdot1_aux(isel,1:2)
     rrdot2(isel,1:2)=rrdot2_aux(isel,1:2)
  ENDDO
  qqd1(1:3,1) = pe(1:3,1)
  qqd1(1:3,2) = pde(1:3,1)
  qqd2(1:3,1) = pe(1:3,2)
  qqd2(1:3,2) = pde(1:3,2)
  deltaom=0.d0
  deltaell=0.d0
  DO isel=1,nprel
     deltaom(isel)=pridif(elkep1(isel)%coord(5),elkep2(isel)%coord(5))
     n1(isel)= sqrt(gmvar)*(elkep1(isel)%coord(1)**(-1.5d0))
     ell2_to_t1=elkep2(isel)%coord(6) - &
          & n1(isel)*(elkep2(isel)%t-elkep1(isel)%t)
     deltaell(isel)=pridif(elkep1(isel)%coord(6),ell2_to_t1)
  ENDDO
!  WRITE(*,*) '========propagation of covariance============'
  DO isel=1,nprel
     CALL derphi(att1,att2,rrdot1(isel,1:2),rrdot2(isel,1:2), &
          & qqd1,qqd2,der_phi_att,der_phi_r)
! storing values for CHECKDER 
     IF(checkd)THEN
        der_phi_att_isel(isel,1:4,1:8) = der_phi_att(1:4,1:8) 
        der_phi_r_isel(isel,1:4,1:4) = der_phi_r(1:4,1:4)
     ENDIF
! ------------------------------------------------------
     CALL matin(der_phi_r,det_phi,4,0,4,ising,1)
     IF (ising.eq.1) THEN 
        WRITE(*,*)'ERROR! singular matrix'
     ELSE
        der_rrd_att=0.d0
        DO i=1,4
           DO j=1,8
              DO k=1,4
                 der_rrd_att(i,j)=der_rrd_att(i,j)- &
                      & der_phi_r(i,k)*der_phi_att(k,j)
              ENDDO
           ENDDO
        ENDDO
     ENDIF

! storing values for CHECKDER 
     IF(checkd)THEN
        der_rrd_att_isel(isel,1:4,1:8) = der_rrd_att(1:4,1:8) 
     ENDIF

! ****************************************************************
! compute the jacobian matrix of KEP with rispect to ATT
! ****************************************************************
! Attributable  elements
    elatt1=undefined_orbit_elem
    elatt1%coo='ATT'
    elatt1%center=3
    elatt1%t=att1%tdtobs- rrdot1(isel,1)/vlight
    elatt1%coord(1:4)=att1%angles
    elatt1%coord(5)= rrdot1(isel,1)
    elatt1%coord(6)= rrdot1(isel,2)
    elatt1%obscode=obscode1 
 !
    elatt2=undefined_orbit_elem
    elatt2%coo='ATT'
    elatt2%center=3
    elatt2%t=att2%tdtobs- rrdot2(isel,1)/vlight
    elatt2%coord(1:4)=att2%angles
    elatt2%coord(5)= rrdot2(isel,1)
    elatt2%coord(6)= rrdot2(isel,2)
    elatt2%obscode=obscode2 

! compute the derivatives of CAR with respect to ATT
    CALL coo_cha(elatt1,'CAR',elcar1_aux,fail_flag,der_car_att1)
    CALL coo_cha(elatt2,'CAR',elcar2_aux,fail_flag,der_car_att2)

! write(*,*)'der1_car(isel,:,:)',der1_car(isel,1:6,1:6)
    der1_att = 0.d0
    der2_att = 0.d0
    DO i=1,6
       DO j=1,6
          DO k=1,6
             der1_att(i,j) =  der1_att(i,j) + &
                  & der1_car(isel,i,k)*der_car_att1(k,j)
             der2_att(i,j) =  der2_att(i,j) + &
                  & der2_car(isel,i,k)*der_car_att2(k,j)
          ENDDO
       ENDDO
    ENDDO
    
! computation of partial derivatives of Delta_1,2
    deltaw_att1= der1_att(5,1:4)+MATMUL(der1_att(5,5:6),der_rrd_att(1:2,1:4)) &
         & -MATMUL(der2_att(5,5:6),der_rrd_att(3:4,1:4))
    deltaw_att2= MATMUL(der1_att(5,5:6),der_rrd_att(1:2,5:8))- &
         & der2_att(5,1:4) &
         & -MATMUL(der2_att(5,5:6),der_rrd_att(3:4,5:8))
    deltaell_att1= der1_att(6,1:4)+ &
         & MATMUL(der1_att(6,5:6),der_rrd_att(1:2,1:4)) &
         & -MATMUL(der2_att(6,5:6),der_rrd_att(3:4,1:4)) & 
         & -1.5d0*(n1(isel)/elkep1(isel)%coord(1))*(der1_att(1,1:4)+ &
         & MATMUL(der1_att(1,5:6),der_rrd_att(1:2,1:4)))*(elatt2%t-elatt1%t)+&
         & (n1(isel)/vlight)*(der_rrd_att(1,1:4)-der_rrd_att(3,1:4))
    deltaell_att2= MATMUL(der1_att(6,5:6),der_rrd_att(1:2,5:8))- &
         & der2_att(6,1:4)&
         & -MATMUL(der2_att(6,5:6),der_rrd_att(3:4,5:8)) & 
         & -1.5d0*(n1(isel)/elkep1(isel)%coord(1))* &
         & MATMUL(der1_att(1,5:6),der_rrd_att(1:2,5:8))*(elatt2%t-elatt1%t)+&
         & (n1(isel)/vlight)*(der_rrd_att(1,5:8)-der_rrd_att(3,5:8))
    
    deltawl_att(1,1:4)= deltaw_att1
    deltawl_att(1,5:8)= deltaw_att2
    deltawl_att(2,1:4)= deltaell_att1
    deltawl_att(2,5:8)= deltaell_att2

    IF(checkd)THEN
       deltawl_att_isel(isel,1:2,1:8) = deltawl_att(1:2,1:8)
    ENDIF

! ***************************************************************
! covariance matrix of the attributables
    gamma_att=0.d0                    
    gamma_att(1:4,1:4)=att1%g
    gamma_att(5:8,5:8)=att2%g
! marginal covariance matrix
    DO i=1,8
       DO j=1,2
          tdeltawl_att(i,j)=deltawl_att(j,i)
       ENDDO
    ENDDO

    CALL matrix_mult(deltawl_att,gamma_att,2,8,8,mat2x8)
    CALL matrix_mult(mat2x8,tdeltawl_att,2,8,2,gamma_delta)
!    gamma_delta=MATMUL(deltawl_att,MATMUL(gamma_att,tdeltawl_att))

    CALL tchinv(gamma_delta,2,c_delta,ws2,ierr)
    IF(ierr.gt.0)THEN
       WRITE(*,*)'ERROR: matrix inversion unsuccesful!,ierr=',ierr
    ENDIF
!    CALL MATIN_qp(c_delta,det_gamma,2,0,2,ising,1) !wrong

    IF(PRESENT(deltaoml))THEN
       deltaoml_aux(isel,1) = deltaom(isel)
       deltaoml_aux(isel,2) = deltaell(isel)
       IF(PRESENT(norm))THEN
! identification norm
          idnorm(isel)=0.d0
          tmpsum=0.d0
          DO j=1,2
             DO k=1,2
                tmpsum(j)= tmpsum(j) + c_delta(j,k)*deltaoml_aux(isel,k)
             ENDDO
             idnorm(isel)=idnorm(isel) + deltaoml_aux(isel,j)*tmpsum(j)
          ENDDO
          idnorm(isel)=sqrt(idnorm(isel))
          IF(verb_2body.gt.10)THEN    
             write(*,*)'norm of Delta_1,2=',idnorm(isel),', isel=', isel
          ENDIF
       ENDIF
    ENDIF   

    IF(PRESENT(orbunc).and.PRESENT(norm).and.idnorm(isel).lt.maxnorm)THEN
! covariance matrices of the orbits at time t_1
       gammatmp=0.d0
       CALL matrix_mult(gamma_att,TRANSPOSE(der_rrd_att(1:2,1:8)), &
            & 8,8,2,gammatmp)
       gamma_R1=0.d0
       CALL matrix_mult(der_rrd_att(1:2,1:8),gammatmp,2,8,2,gamma_R1)
       gamma_a1R1=0.d0
       CALL matrix_mult(gamma_att(1:4,1:4),TRANSPOSE(der_rrd_att(1:2,1:4)),&
            & 4,4,2,gamma_a1R1)
       
       orbunc_aux(isel)%g(1:4,1:4)=gamma_att(1:4,1:4)
       orbunc_aux(isel)%g(5:6,5:6)=gamma_R1
       orbunc_aux(isel)%g(1:4,5:6)=gamma_a1R1
       orbunc_aux(isel)%g(5:6,1:4)=TRANSPOSE(gamma_a1R1)
       
! matrix inversion
       CALL tchinv(orbunc_aux(isel)%g,6,orbunc_aux(isel)%c,ws6,ierr)
       orbunc_aux(isel)%ndim=6
       IF(ierr.eq.0)THEN
          orbunc_aux(isel)%succ=.true.
       ELSE
          WRITE(*,*)'ERROR: matrix inversion unsuccesful!,ierr=',ierr
       ENDIF

! derivative matrix: dKEP/dATT
       CALL matrix_mult(der1_car(isel,1:6,1:6),&
            & der_car_att1,6,6,6, dKEP_dATT_1(isel,1:6,1:6))

! uncertainty for Keplerian elements
       CALL convertunc(orbunc_aux(isel),dKEP_dATT_1(isel,1:6,1:6),&
            & orbunc_kep(isel))
    ENDIF

 ENDDO

! ****************************
! CHECK DERIVATIVES 
 IF(checkd)THEN
    DO isel=1,nprel
 ! values of the Phi map at Att
       phi(isel,1:3)=Dvec(1:3,1)*rrdot1(isel,2)-Dvec(1:3,2)*rrdot2(isel,2)+ &
            & Evec(1:3,1)*rrdot1(isel,1)**2-&
            & Evec(1:3,2)*rrdot2(isel,1)**2+Fvec(1:3,1)*rrdot1(isel,1)- &
            & Fvec(1:3,2)*rrdot2(isel,1)+&
            & Gvec(1:3,1)-Gvec(1:3,2) 

       phi(isel,4)=rrdot1(isel,2)**2+c1(1)*rrdot1(isel,2)+ &
            & c2(1)*rrdot1(isel,1)**2+ &
            & c3(1)*rrdot1(isel,1)+c4(1)-&
            & 2*(gmvar/sqrt(rrdot1(isel,1)**2+c5(1)*rrdot1(isel,1)+c0(1)))-&
            & (rrdot2(isel,2)**2+c1(2)*rrdot2(isel,2)+ &
            & c2(2)*rrdot2(isel,1)**2+ c3(2)*rrdot2(isel,1)+c4(2)-&
            & 2*(gmvar/sqrt(rrdot2(isel,1)**2+c5(2)*rrdot2(isel,1)+c0(2))))

    ENDDO
! increments
    incr=0.d0 ! initialization
    WRITE(*,*)'***  CHECK DERIVATIVES  ***'
    OPEN(1,file='increments',status='old')
    READ(1,*)
    READ(1,*)incr
    CLOSE(1)
    write(*,*)'increments:',incr
    att1h = att1
    att2h = att2
    DO i=1,8
       IF(incr(i).ne.0.d0)THEN
          IF(i.le.4)THEN
             att1h%angles(i) = att1%angles(i)+incr(i)
          ELSEIF(i.gt.4)THEN
             att2h%angles(i-4) = att2%angles(i-4)+incr(i)
          ENDIF
       ENDIF
    ENDDO
! incremented values of rho,rhod
    CALL orbprelim_2b(att1h,att2h,nprel_h,elcar1_h,elcar2_h,&
         & rrdot1_h,rrdot2_h,QOBS1=qobs1, &
         & QPOBS1=qpobs1,QOBS2=qobs2,QPOBS2=qpobs2)
    IF (nprel .ne. nprel_h) THEN
       WRITE(*,*) 'Check error: nprel has changed.'
    ENDIF
    DO isel=1,nprel
! values of the Phi map at Atth
       phi_atth(1:3)=Dvec(1:3,1)*rrdot1(isel,2)- &
            & Dvec(1:3,2)*rrdot2(isel,2)+ Evec(1:3,1)*rrdot1(isel,1)**2-&
            & Evec(1:3,2)*rrdot2(isel,1)**2+Fvec(1:3,1)*rrdot1(isel,1)- &
            & Fvec(1:3,2)*rrdot2(isel,1)+ Gvec(1:3,1)-Gvec(1:3,2) 

       phi_atth(4)=rrdot1(isel,2)**2 + c1(1)*rrdot1(isel,2)+ &
            & c2(1)*rrdot1(isel,1)**2 + c3(1)*rrdot1(isel,1)+&
            & c4(1)-2*(gmvar/sqrt(rrdot1(isel,1)**2 + &
            & c5(1)*rrdot1(isel,1)+c0(1)))-&
            & (rrdot2(isel,2)**2 + c1(2)*rrdot2(isel,2)+ &
            & c2(2)*rrdot2(isel,1)**2 + c3(2)*rrdot2(isel,1)+&
            & c4(2)-2*(gmvar/sqrt(rrdot2(isel,1)**2 + &
            & c5(2)*rrdot2(isel,1)+c0(2))))
              
       DO i=1,8
          IF(incr(i).ne.0.d0)THEN
             DO j=1,4
                write(*,113) 'phi(j)(atth)-phi(j)(att)/incr(i) ',&
                     & (phi_atth(j)-phi(isel,j))/incr(i) 
                write(*,113)'(d phi(j))/(d att(j)) ',der_phi_att_isel(isel,j,i)
             ENDDO
             EXIT
          ENDIF
       ENDDO
113    FORMAT(a35,3x,f18.10)
       CALL coo_cha(elcar1_h(isel),'KEP',elkep1_h,fail_flag)
       CALL coo_cha(elcar2_h(isel),'KEP',elkep2_h,fail_flag)       
!       rrd1_h(1:2) = rrdot1_h(isel,1:2)
!       rrd2_h(1:2) = rrdot2_h(isel,1:2)
       deltaom_h=pridif(elkep1_h%coord(5),elkep2_h%coord(5))
       n1_h= sqrt(gmvar)*(elkep1_h%coord(1)**(-1.5d0))
       ell2_to_t1_h=elkep2_h%coord(6) - n1_h*(elkep2_h%t-elkep1_h%t)
       deltaell_h=pridif(elkep1_h%coord(6),ell2_to_t1_h)

! *************************************************************
       DO i=1,8
          IF(incr(i).ne.0.d0)THEN
             write(*,113)' (rho1(atth)-rho1(att))/incr(i)',&
                  &(rrdot1_h(isel,1)-rrdot1(isel,1))/incr(i)
             write(*,113)'(d rho1)/(d att_i)=',der_rrd_att_isel(isel,1,i)
             
             write(*,113)' (rhodot1(atth)-rho1dot(att))/incr(i)',&
                  &(rrdot1_h(isel,2)-rrdot1(isel,2))/incr(i)
             write(*,113)'(d rhodot1)/(d att_i)=',der_rrd_att_isel(isel,2,i)
             
             write(*,113)' (deltaom(atth)-deltaom(att))/incr(i)',&
                  & pridif(deltaom_h,deltaom(isel))/incr(i)
             write(*,113)'(d deltaom)/(d alpha1)=',deltawl_att_isel(isel,1,i)
             
             write(*,113)' (deltaell(atth)-deltaell(att))/incr(i)',&
                  & pridif(deltaell_h,deltaell(isel))/incr(i)
             write(*,113)'(d deltaell)/(d alpha1)=',deltawl_att_isel(isel,2,i)
          ENDIF
       ENDDO
    ENDDO
 ENDIF

! =================================================================
 IF(PRESENT(norm))THEN
    orb_srt(1:nprel)=0
    IF(verb_2body.gt.10)THEN
       WRITE(*,*)'... selecting good orbits and reordering ... '
       WRITE(*,*)''
    ENDIF
    CALL heapsort(idnorm(1:nprel),nprel,indsrt(1:nprel)) 
!    DO isel=1,nprel
!       write(*,*)'indsrt,dist',indsrt(isel),idnorm(indsrt(isel))
!    ENDDO

! preliminary orbits
    nsol_tmp=0
    DO isel=1,nprel
       isrt=indsrt(isel)
       IF(idnorm(isrt).le.maxnorm)THEN
          nsol_tmp=nsol_tmp+1
          ! preliminary orbits
          orbprel(isel)=undefined_orbit_elem
          IF(rhs.EQ.2)THEN
             orbprel(isel)%center=3
          ENDIF
          orbprel(isel)%coord(1:6)= elkep1(isrt)%coord(1:6) 
          orbprel(isel)%t= elkep1(isrt)%t
          orbprel(isel)%coo='KEP'
          ! orbit uncertainty
          orbunc(isel)=orbunc_kep(isrt)
          ! deltaoml
          deltaoml(isel,1:2)=deltaoml_aux(isrt,1:2)
          ! identification norm
          norm(isel)=idnorm(isrt)
       ENDIF
    ENDDO
    nprel=nsol_tmp   
    
 ENDIF
 

END SUBROUTINE twobodyint

! ========================================================================
SUBROUTINE orbprelim_2b(att1,att2,nprel,elcar1,elcar2,rrdot1,rrdot2, &
     &nsolstep,qobs1,qpobs1,qobs2,qpobs2)
!  USE fund_const
  USE orbit_elements
  USE attributable
  USE output_control
  USE reference_systems, ONLY: observer_position,pvobs
  USE station_coordinates, ONLY: statcode
  IMPLICIT NONE
  TYPE(attrib),INTENT(IN) :: att1,att2
  INTEGER, INTENT(OUT):: nprel ! number of preliminary orbits accepted
  TYPE(orbit_elem),INTENT(OUT),DIMENSION(nprelx):: elcar1,elcar2 !CAR elems
  REAL(KIND=dkind),INTENT(OUT),DIMENSION(nprelx,2):: rrdot1,rrdot2!rdist/vel
  INTEGER,INTENT(OUT),DIMENSION(2),OPTIONAL:: nsolstep
! interpolated geocentric position and velocity of the observer (equatorial);
  REAL(KIND=dkind), INTENT(IN), DIMENSION(3),OPTIONAL :: qobs1,qpobs1
  REAL(KIND=dkind), INTENT(IN), DIMENSION(3),OPTIONAL :: qobs2,qpobs2
! ----------- end interface --------------------------------------------
!  REAL(KIND=dkind), PARAMETER :: mindist=2.q-2
  REAL(KIND=dkind) :: mindist
  REAL(KIND=dkind) :: t1,t2 ! time of the attributables MJD
  REAL(KIND=dkind) :: norm2
  REAL(KIND=dkind),PARAMETER :: domx=2.d0,norm2x=1.d-3
  LOGICAL :: duplfl
! orbital elements
  REAL(KIND=dkind),DIMENSION(6) :: eqtmp1,eqtmp2 ! auxiliary
  TYPE(orbit_elem) :: elem1,elem2   ! EQU elements of the Earth
  TYPE(orbit_elem) :: ecpl1,ecpl2   ! CAR elements of the Earth
  TYPE(orbit_elem) :: elatt1,elatt2 ! attributable elements of the asteroid
  REAL(KIND=dkind), DIMENSION(poldeg,6) :: cfr_el1,cfr_el2 ! auxiliary, Keplerian
!
  REAL(KIND=dkind),DIMENSION(2) :: alpha,delta,ad,dd ! attributables
  REAL(KIND=dkind), DIMENSION(3) :: pos1,pos2     ! observer position vectors
  REAL(KIND=dkind),DIMENSION(3) :: vel1,vel2     ! observer velocity vectors

! auxiliary --------------------------------------------------------
  REAL(KIND=dkind),DIMENSION(3,2) :: pe_pde,pe_rhat
  REAL(KIND=dkind),DIMENSION(3,2) :: pe_rhata,pe_rhatd,rhat_rhata,rhat_rhatd
  REAL(KIND=dkind),DIMENSION(3,2) :: rhat_pde
!  REAL(KIND=qkind) :: q,ecc,inc,omnod,omeg,tperi
  INTEGER :: j1,j2,j3
! -------------------------------------------------------------------
!  REAL(KIND=qkind),DIMENSION(3,2) :: A,B,C,D
  REAL(KIND=dkind) :: det12,det13,det23,detD,maxdet
  REAL(KIND=qkind) :: coe_f(0:2,0:2),coe_g(0:2,0:2)
  REAL(KIND=qkind) :: coe_rd1(0:2,0:2),coe_rd2(0:2,0:2)
  REAL(KIND=qkind) :: coe_rd1_quad(0:4,0:4),coe_rd2_quad(0:4,0:4)
  REAL(KIND=qkind) :: coe_F1(0:4,0:4),coe_F2(0:4,0:4)
  REAL(KIND=qkind) :: coe_G1(0:2,0:2),coe_G2(0:2,0:2)
  REAL(KIND=qkind) :: coe_G1G2(0:4,0:4),coe_fact1(0:4,0:4)
  REAL(KIND=qkind) :: coe_diffF1F2(0:4,0:4),coe_dF1F2_quad(0:8,0:8)
  REAL(KIND=qkind) :: coe_fact2tmp(0:12,0:12),coe_addG1G2(0:2,0:2)
  REAL(KIND=qkind) :: coe_fact2(0:12,0:12),coe_fact2_quad(0:24,0:24)
  REAL(KIND=qkind) :: coe_p(0:24,0:24),coe_q(0:2,0:2)
  REAL(KIND=qkind) :: ncoe_p(0:24,0:24),ncoe_q(0:2,0:2)!normalized, qkind
  REAL(KIND=dkind) :: ncoe_p_DP(0:24,0:24),ncoe_q_DP(0:2,0:2)!normalized,dkind
! evaluations of the coefficients a(rho2), b(rho2) in the 64th roots of unity
  REAL(KIND=qkind) :: aa(1:21,1:Nev),bb(1:3,1:Nev) !auxiliary 
  REAL(KIND=qkind) :: ev_a(1:21,1:Nev),ev_b(1:3,1:Nev),ev_x(1:Nev)
  COMPLEX(KIND=qkind),DIMENSION(21,Nev) :: compl_eva
  COMPLEX(KIND=qkind),DIMENSION(3,Nev) :: compl_evb
  COMPLEX(KIND=qkind),DIMENSION(Nev) :: compl_evx
  COMPLEX(KIND=qkind) :: ev_Sylv_j(22,22),ev_Sylv(22,22,Nev)
  COMPLEX(KIND=qkind) :: deteval(Nev)
  REAL(KIND=qkind) :: maxdeteval,norm_res,prea,pimg ! auxiliary

! computation of determinants
  COMPLEX(KIND=qkind) :: det_ev_Sylv_j,det_ev_Sylv(Nev)
  REAL(KIND=qkind) ::    ev_detSylv(Nev)
! interpolation
  INTEGER :: ncoe,ncoe1
  REAL(KIND=qkind) :: detS_coe(0:Nev-1)
! for computing roots
  INTEGER :: pdeg
  REAL(KIND=qkind) :: roots(poldeg) !real roots 
  REAL(KIND=qkind) :: rho2sol(poldeg) !real roots 
  INTEGER :: nroots ! number of real roots
  REAL(KIND=qkind) :: rho1pos(poldeg),rho2pos(poldeg) !positive real roots
  REAL(KIND=qkind) :: rho1_tmp(poldeg),rho2_tmp(poldeg) !positive real roots
  INTEGER :: nsol,nsol_tmp ! number of positive real roots
!
  LOGICAL :: hzflag ! hzflag = .true.  OK!
                    ! hzflag = .false. abs(root)>10^5
  LOGICAL :: multfl ! multfl = .true.  OK!
                    ! multfl = .false. 0 has multiplicity > 4
  REAL(KIND=dkind) :: rho1dot(poldeg),rho2dot(poldeg)
  REAL(KIND=dkind) :: deltat,enne2
!
  REAL(KIND=qkind) :: rho1,rho2,drho1,drho2 ! auxiliary, for newton-raphson
  REAL(KIND=qkind) :: pk,qk ! auxiliary, for newton-raphson
  REAL(KIND=qkind) :: mincoep,maxcoep,pnorm
  REAL(KIND=qkind) :: mincoeq,maxcoeq,qnorm
!
  INTEGER :: rd1quamod,rd2quamod
  INTEGER :: G1G2mod,dF1F2mod,dF1F2quamod,f2tmpmod,addG1G2mod
  INTEGER :: fact2mod,f2quamod
  INTEGER :: pmod,qmod ! degree of the polynomials p,q
! 
  REAL(KIND=qkind) :: F1eval,F2eval,fact2eval! to discard spurious solutions
  REAL(KIND=qkind) :: G1eval,G2eval,dF1F2eval
!
  INTEGER :: isel ! selected index
  INTEGER :: i,j,h,k,n ! loop indexes
  INTEGER, PARAMETER :: nitmax=10
  INTEGER :: fail_flag
  REAL(KIND=dkind) :: pridif,primea,ell2_to_t1

  REAL(KIND=dkind) :: xo1(3),vo1(3)
  REAL(KIND=dkind) :: xea1(6),xoec1(3),voec1(3),xeaeq1(6),xogeo1(3)
  REAL(KIND=dkind) :: xo2(3),vo2(3)
  REAL(KIND=dkind) :: xea2(6),xoec2(3),voec2(3),xeaeq2(6),xogeo2(3)
! conversion factor and vectors
  REAL(KIND=dkind) :: convdist,gmvarc
!  REAL(KIND=dkind) :: qobs1c(3),qpobs1c(3),qobs2c(3),qpobs2c(3)
!scaled ang mom coeffs
  REAL(KIND=dkind),DIMENSION(3,2) :: Dvecc,Evecc,Fvecc,Gvecc
!scaled energy coeffs !
  REAL(KIND=dkind),DIMENSION(2) :: c0c,c1c,c2c,c3c,c31c,c32c,c4c,c5c
  REAL(KIND=dkind) :: maxv,minv,rho1store,rho2store

! =============================================================

  IF(rhs.EQ.1)THEN
     convdist=1.d0
  ELSEIF(rhs.EQ.2)THEN
     convdist=150.d6/42.d3 ! unit 42000 km
!     convdist=150.d6/35.d3
!     convdist=150.d6/6378.388d0
!     convdist=100.d0
  ELSE
     WRITE(*,*)'twobodyint: value not allowed for rhs:',rhs
     STOP
  ENDIF

  gmvarc=gmvar*convdist**3

  IF(rhs.eq.1)THEN
     mindist=2.q-2*convdist
  ELSEIF(rhs.eq.2)THEN
     mindist=1.q-10*convdist
  ELSE
     WRITE(*,*)'orbprelim_2b: rhs=', rhs
     STOP
  ENDIF

  verb_2body=11
  WRITE(*,*)'verb_2body=',verb_2body

  alpha(1) = att1%angles(1)
  delta(1) = att1%angles(2)
  ad(1) = att1%angles(3)
  dd(1) = att1%angles(4)
  CALL statcode(att1%obscod,obscode1)
  alpha(2) = att2%angles(1)
  delta(2) = att2%angles(2)
  ad(2) = att2%angles(3)
  dd(2) = att2%angles(4)
  CALL statcode(att2%obscod,obscode2)
  IF(verb_2body.gt.10) THEN
     WRITE(*,*)'obscode1=',obscode1
     WRITE(*,*)'obscode2=',obscode2
  ENDIF
  t1 = att1%tdtobs
  t2 = att2%tdtobs

  DO i=1,2
     cosa(i)=cos(alpha(i))
     sina(i)=sin(alpha(i))
     cosd(i)=cos(delta(i))
     sind(i)=sin(delta(i))
  ENDDO
  
  DO i=1,2
     rhohat(1,i) = cosa(i)*cosd(i)
     rhohat(2,i) = sina(i)*cosd(i)
     rhohat(3,i) = sind(i)
     rhohat_a(1,i) = -sina(i)*cosd(i)
     rhohat_a(2,i) = cosa(i)*cosd(i)
     rhohat_a(3,i) = 0.d0
     rhohat_d(1,i) = -cosa(i)*sind(i)
     rhohat_d(2,i) = -sina(i)*sind(i)
     rhohat_d(3,i) = cosd(i)
  ENDDO

  IF(rhs.eq.1)THEN
! Earth elements at time t1
     CALL earcar(t1,xea1,1) ! earth coord, heliocentric ecliptic
  ELSEIF(rhs.eq.2)THEN
     xea1=0.d0
  ENDIF

  IF(PRESENT(qobs1).and.PRESENT(qpobs1))THEN
! interpolated obs pos. and vel. 
! (heliocentric/geocentric (dep. on rhs)equatorial)
     pe(1:3,1) = qobs1+ MATMUL(roteceq,xea1(1:3)) 
     pde(1:3,1)= qpobs1+MATMUL(roteceq,xea1(4:6))
  ELSE
     IF(rhs.EQ.1)THEN
        CALL pvobs(t1,obscode1,xoec1,voec1) ! obs. coord, geocentric ecliptic
        xoec1=xoec1+xea1(1:3) ! obs coord, heliocentric ecliptic
        voec1=voec1+xea1(4:6) ! id velocity
        pe(1:3,1)=MATMUL(roteceq,xoec1) ! obs coord, heliocentric equatorial
        pde(1:3,1)=MATMUL(roteceq,voec1) ! id velocity
     ELSEIF(rhs.EQ.2)THEN
        CALL observer_position(t1,pe(1:3,1),pde(1:3,1),obscode1)
     ENDIF
  ENDIF
  WRITE(*,*)'Equatorial observer position at t1:', pe(1:3,1)
  WRITE(*,*)'Equatorial observer velocity at t1:', pde(1:3,1)

  IF(rhs.eq.1)THEN
! Earth elements at time t2
     CALL earcar(t2,xea2,1) ! earth coord, heliocentric ecliptic
  ELSEIF(rhs.eq.2)THEN
     xea2=0.d0
  ENDIF
  IF(PRESENT(qobs2).and.PRESENT(qpobs2))THEN
! interpolated obs pos. and vel. (heliocentric equatorial)
     pe(1:3,2) = qobs2+ MATMUL(roteceq,xea2(1:3)) 
     pde(1:3,2)= qpobs2+MATMUL(roteceq,xea2(4:6)) 
  ELSE
     IF(rhs.EQ.1)THEN
        CALL pvobs(t2,obscode2,xoec2,voec2) ! obs. coord, geocentric ecliptic
        xoec2=xoec2+xea2(1:3) ! obs coord, heliocentric ecliptic
        voec2=voec2+xea2(4:6) ! id velocity
        pe(1:3,2)=MATMUL(roteceq,xoec2) ! obs coord, heliocentric equatorial
        pde(1:3,2)=MATMUL(roteceq,voec2) ! id velocity
     ELSEIF(rhs.EQ.2)THEN
        CALL observer_position(t2,pe(1:3,2),pde(1:3,2),obscode2)
     ENDIF
  ENDIF
  WRITE(*,*)'Equatorial observer position at t2:', pe(1:3,2)
  WRITE(*,*)'Equatorial observer velocity at t2:', pde(1:3,2)

  DO j=1,2
     CALL cross_prod(pe(1:3,j),pde(1:3,j),pe_pde(1:3,j))
     CALL cross_prod(pe(1:3,j),rhohat_a(1:3,j),pe_rhata(1:3,j))
     CALL cross_prod(pe(1:3,j),rhohat_d(1:3,j),pe_rhatd(1:3,j))
     CALL cross_prod(rhohat(1:3,j),pde(1:3,j),rhat_pde(1:3,j))
     CALL cross_prod(rhohat(1:3,j),rhohat_a(1:3,j),rhat_rhata(1:3,j))
     CALL cross_prod(rhohat(1:3,j),rhohat_d(1:3,j),rhat_rhatd(1:3,j))
     CALL cross_prod(pe(1:3,j),rhohat(1:3,j),pe_rhat(1:3,j))
  ENDDO

! Admissible Region defining coefficients
!  OPEN(4,file='AR_coeff',status='unknown') 
  DO j=1,2
     c0(j) = DOT_PRODUCT(pe(1:3,j),pe(1:3,j))
     c0c(j)=c0(j)*convdist**2
     c1(j) = 2.q0*DOT_PRODUCT(pde(1:3,j),rhohat(1:3,j))
     c1c(j)=c1(j)*convdist
     c2(j) = (ad(j)*cos(delta(j)))**2 + dd(j)**2
     c2c(j)=c2(j)
     c31(j) = 2.q0*DOT_PRODUCT(pde(1:3,j),rhohat_a(1:3,j))
     c31c(j)=c31(j)*convdist
     c32(j) = 2.q0*DOT_PRODUCT(pde(1:3,j),rhohat_d(1:3,j))
     c32c(j)=c32(j)*convdist
     c3(j) = ad(j)*c31(j) + dd(j)*c32(j);
     c3c(j) = ad(j)*c31c(j) + dd(j)*c32c(j);
     c4(j) = DOT_PRODUCT(pde(1:3,j),pde(1:3,j))
     c4c(j)=c4(j)*convdist**2
     c5(j) = 2.q0*DOT_PRODUCT(pe(1:3,j),rhohat(1:3,j))
     c5c(j)=c5(j)*convdist
!     WRITE(4,112) c0(j),c1(j),c2(j),c3(j),c4(j),c5(j),gmvar
     IF(verb_2body.gt.10)THEN
        WRITE(*,*)'cj coeff, gmvarc'
        WRITE(*,112) c0c(j),c1c(j),c2c(j),c3c(j),c4c(j),c5c(j),gmvarc
!        WRITE(4,112) c0(j),c1(j),c2(j),c3(j),c4(j),c5(j),gmvar
     ENDIF
112  FORMAT(7(e25.20,1x))
  ENDDO
!  CLOSE(4)

! Angular Momentum curve
  DO j = 1,2
     Gvec(1:3,j)= pe_pde(1:3,j)
     Gvecc(1:3,j)=Gvec(1:3,j)*convdist**2
     Fvec(1:3,j)= ad(j)*pe_rhata(1:3,j)+dd(j)*pe_rhatd(1:3,j)+rhat_pde(1:3,j)
     Fvecc(1:3,j)=Fvec(1:3,j)*convdist
     Evec(1:3,j)= ad(j)*rhat_rhata(1:3,j)+dd(j)*rhat_rhatd(1:3,j)
     Evecc(1:3,j)=Evec(1:3,j)
     Dvec(1:3,j)= pe_rhat(1:3,j)
     Dvecc(1:3,j)=Dvec(1:3,j)*convdist
  ENDDO
  det12= -Dvecc(1,1)*Dvecc(2,2)+Dvecc(2,1)*Dvecc(1,2)
  det13= -Dvecc(1,1)*Dvecc(3,2)+Dvecc(3,1)*Dvecc(1,2)
  det23= -Dvecc(2,1)*Dvecc(3,2)+Dvecc(3,1)*Dvecc(2,2)
  maxdet=MAX(abs(det12),abs(det13),abs(det23))
  IF(maxdet.eq.abs(det12))THEN
     j1 = 1; j2 = 2; j3 = 3
     detD=det12
  ELSEIF(maxdet.eq.abs(det13))THEN
     j1 = 1; j2 = 3; j3 = 2
     detD=det13
  ELSEIF(maxdet.eq.abs(det23))THEN
     j1 = 2; j2 = 3; j3 = 1
     detD=det23
  ENDIF

! coefficients of rd1,rd2 as functions of rho1,rho2
  coe_f=0.q0; coe_g=0.q0
  coe_f(0,0) = (Gvecc(j1,2)-Gvecc(j1,1)) 
  coe_f(1,0) = -Fvecc(j1,1)
  coe_f(2,0) = -Evecc(j1,1)
  coe_f(0,1) =  Fvecc(j1,2)
  coe_f(0,2) =  Evecc(j1,2)

  coe_g(0,0) = (Gvecc(j2,2)-Gvecc(j2,1)) 
  coe_g(1,0) = -Fvecc(j2,1)
  coe_g(2,0) = -Evecc(j2,1)
  coe_g(0,1) =  Fvecc(j2,2)
  coe_g(0,2) =  Evecc(j2,2)

  coe_rd1=0.q0; coe_rd2=0.q0
  DO h=0,2
     DO k=0,2 
        IF(h.eq.0.or.k.eq.0)THEN
           coe_rd1(h,k) = (-coe_f(h,k)*Dvecc(j2,2) + coe_g(h,k)*Dvecc(j1,2))/detD
           coe_rd2(h,k) = ( Dvecc(j1,1)*coe_g(h,k) - Dvecc(j2,1)*coe_f(h,k))/detD
        ENDIF
     ENDDO   
  ENDDO

  CALL poly_product_QP(2,2,coe_rd1,coe_rd1,rd1quamod,coe_rd1_quad)
  CALL poly_product_QP(2,2,coe_rd2,coe_rd2,rd2quamod,coe_rd2_quad)

  coe_F1=coe_rd1_quad
  coe_F2=coe_rd2_quad

  coe_F1(0,0) = coe_F1(0,0) + c1c(1)*coe_rd1(0,0) + c4c(1)
  coe_F1(1,0) = coe_F1(1,0) + c1c(1)*coe_rd1(1,0) + c3c(1)
  coe_F1(2,0) = coe_F1(2,0) + c1c(1)*coe_rd1(2,0) + c2c(1)
  coe_F1(0,1) = coe_F1(0,1) + c1c(1)*coe_rd1(0,1)
  coe_F1(0,2) = coe_F1(0,2) + c1c(1)*coe_rd1(0,2)

  coe_F2(0,0) = coe_F2(0,0) + c1c(2)*coe_rd2(0,0) + c4c(2)
  coe_F2(1,0) = coe_F2(1,0) + c1c(2)*coe_rd2(1,0) 
  coe_F2(2,0) = coe_F2(2,0) + c1c(2)*coe_rd2(2,0) 
  coe_F2(0,1) = coe_F2(0,1) + c1c(2)*coe_rd2(0,1) + c3c(2)
  coe_F2(0,2) = coe_F2(0,2) + c1c(2)*coe_rd2(0,2) + c2c(2)

  coe_G1=0.q0; coe_G2=0.q0
  coe_G1(0,0) = c0c(1)
  coe_G1(1,0) = c5c(1)
  coe_G1(2,0) = 1.q0

  coe_G2(0,0) = c0c(2)
  coe_G2(0,1) = c5c(2)
  coe_G2(0,2) = 1.q0

  CALL poly_product_QP(2,2,coe_G1,coe_G2,G1G2mod,coe_G1G2)

  coe_fact1 = 64.q0*gmvarc**4.q0*coe_G1G2

  CALL poly_sum_QP(4,4,coe_F1,-coe_F2,dF1F2mod,coe_diffF1F2)

  CALL poly_product_QP(dF1F2mod,dF1F2mod,coe_diffF1F2,coe_diffF1F2, &
       & dF1F2quamod,coe_dF1F2_quad)

  CALL poly_product_QP(dF1F2quamod,G1G2mod,coe_dF1F2_quad,coe_G1G2, &
       & f2tmpmod,coe_fact2tmp)

  CALL poly_sum_QP(2,2,coe_G1,coe_G2,addG1G2mod,coe_addG1G2)

  coe_addG1G2 = -4.q0*gmvarc**2.q0*coe_addG1G2

  CALL poly_sum_QP(f2tmpmod,addG1G2mod,coe_fact2tmp,coe_addG1G2, &
       & fact2mod,coe_fact2)

  CALL poly_product_QP(fact2mod,fact2mod,coe_fact2,coe_fact2, &
       & f2quamod,coe_fact2_quad)

! ********************************************************************
! p(rho1,rho2) polynomial
  CALL poly_sum_QP(G1G2mod,f2quamod,coe_fact1,-coe_fact2_quad,pmod,coe_p)
! ********************************************************************

! normalization of coe_p
  pnorm = coe_p(0,0)
!  pnorm=1.d0 !dummy  
  ncoe_p =  coe_p/pnorm
  ncoe_p_DP = ncoe_p 

  maxcoep=maxval(abs(ncoe_p(0:20,0:20)))
  mincoep=1.q20 ! initialization
  DO h=0,20
     DO k=0,20
        IF(ncoe_p(h,k).ne.0.q0)THEN
!           WRITE(*,*)'ncoe_p(',h,k,')=',ncoe_p(h,k)
           WRITE(37,*)  h,k,ncoe_p(h,k)
           IF(abs(ncoe_p(h,k)).lt.mincoep)THEN
              mincoep=abs(ncoe_p(h,k))
           ENDIF
        ENDIF
     ENDDO
  ENDDO
100 FORMAT(i2,2x,i2,e28.20)
  write(*,*)'maxcoep=',maxcoep,'mincoep=',mincoep

! *********************************************************************
! q(rho1,rho2) polynomial
  coe_q=0.q0
  qmod = 2
  coe_q(0,0)= Dvecc(j3,1)*coe_rd1(0,0)-Dvecc(j3,2)*coe_rd2(0,0)+ &
       & Gvecc(j3,1)-Gvecc(j3,2) 
  coe_q(1,0)= Dvecc(j3,1)*coe_rd1(1,0)-Dvecc(j3,2)*coe_rd2(1,0) + Fvecc(j3,1) 
  coe_q(2,0)= Dvecc(j3,1)*coe_rd1(2,0)-Dvecc(j3,2)*coe_rd2(2,0) + Evecc(j3,1) 
  coe_q(0,1)= Dvecc(j3,1)*coe_rd1(0,1)-Dvecc(j3,2)*coe_rd2(0,1) - Fvecc(j3,2)
  coe_q(0,2)= Dvecc(j3,1)*coe_rd1(0,2)-Dvecc(j3,2)*coe_rd2(0,2) - Evecc(j3,2)
! *********************************************************************

! normalization of coe_q
  qnorm = coe_q(0,0)
!  qnorm=1.d0 !dummy
  ncoe_q =  coe_q/qnorm
  ncoe_q_DP = ncoe_q 

  maxcoeq=maxval(abs(ncoe_q(0:2,0:2)))
  mincoeq=1.q20 
  DO h=0,2
     DO k=0,2
        IF(ncoe_q(h,k).ne.0.q0)THEN
!           WRITE(*,*)'coe_q(',h,k,')=',ncoe_q(h,k)
           WRITE(38,*)  h,k,ncoe_q(h,k)
           IF(abs(ncoe_q(h,k)).lt.mincoeq)THEN
              mincoeq=abs(ncoe_q(h,k))
           ENDIF
       ENDIF
     ENDDO
  ENDDO
  write(*,*)'maxcoeq=',maxcoeq,'mincoeq=',mincoeq

! initialization
  ev_a=0.q0
  ev_b=0.q0
! evaluate the poly coefficients a,b of Sylvester's matrix
  DO h = 1,21
     aa(h,1:21)=ncoe_p(h-1,0:20)
     ev_a(h,1:21)=ncoe_p(h-1,0:20)
     CALL rvfft_QP(ev_a(h,1:Nev),Nev,expo) 
!     write(*,*)'h,ev_a(h,1)=',h,ev_a(h,1)
  ENDDO
  DO h = 1,3
     bb(h,1:3)=ncoe_q(h-1,0:2)
     ev_b(h,1:3)=ncoe_q(h-1,0:2)
     CALL rvfft_QP(ev_b(h,1:Nev),Nev,expo) 
!     write(*,*)'h,ev_b(h,1)=',h,ev_b(h,1)
  ENDDO

! roots of unity
  ev_x=0.q0
  ev_x(2)=1.q0
  CALL rvfft_QP(ev_x(1:Nev),Nev,expo) 
  compl_evx(1) = DCMPLX(ev_x(1),0.q0)
  compl_evx(Nev/2+1) = DCMPLX(ev_x(Nev/2+1),0.q0)
  DO j = 1,Nev/2-1
     compl_evx(j+1) = DCMPLX(ev_x(j+1),ev_x(Nev-j+1))
     compl_evx(Nev/2+j+1) = DCMPLX(ev_x(Nev/2-j+1),-ev_x(Nev/2+1+j))
  ENDDO

! complexification of the result
  DO h = 1,21
    compl_eva(h,1) = QCMPLX(ev_a(h,1),0.q0)
    compl_eva(h,Nev/2+1) = QCMPLX(ev_a(h,Nev/2+1),0.q0)

     DO j = 1,Nev/2-1
        compl_eva(h,j+1) = QCMPLX(ev_a(h,j+1),ev_a(h,Nev-j+1))
        compl_eva(h,Nev/2+j+1) = QCMPLX(ev_a(h,Nev/2-j+1),-ev_a(h,Nev/2+1+j))
     ENDDO
  ENDDO
  DO h = 1,3
    compl_evb(h,1) = QCMPLX(ev_b(h,1),0.q0)
    compl_evb(h,Nev/2+1) = QCMPLX(ev_b(h,Nev/2+1),0.q0)
     DO j = 1,Nev/2-1
        compl_evb(h,j+1) = QCMPLX(ev_b(h,j+1),ev_b(h,Nev-j+1))
        compl_evb(h,Nev/2+j+1) = QCMPLX(ev_b(h,Nev/2-j+1),-ev_b(h,Nev/2+1+j))
     ENDDO
  ENDDO
  
! build Sylvester's matrixes
  DO j =1,Nev
!     write(*,*)'j,compl_eva(1,j)=',j,compl_eva(1,j)
     CALL compsylv48_QP(compl_eva(1:21,j),compl_evb(1:3,j),ev_Sylv_j)
     DO h = 1,22 
        DO k = 1,22 
           ev_Sylv(h,k,j) = ev_Sylv_j(h,k) 
        ENDDO
     ENDDO
! computation of determinants
     CALL cdetcomp22_QP(ev_Sylv_j,det_ev_Sylv_j)
     det_ev_Sylv(j) = det_ev_Sylv_j
!     write(*,*) j,'det_ev_Sylv_j=',det_ev_Sylv_j ! evaluated determinants
  ENDDO

! compute l^2 norm for resultant evaluations
  norm_res=0.q0
  DO j =1,Nev
     prea=QREAL(det_ev_Sylv(j))
     pimg=AIMAG(det_ev_Sylv(j))
     norm_res=norm_res + prea**2 + pimg**2
  ENDDO
  norm_res=sqrt(norm_res)  
!  write(*,*)'norm_res=',norm_res

! -------------------------------------------------------------------
! for the CHECK of RESULTANT EVALUATIONS
!  maxdeteval=maxval(abs(DBLE(det_ev_Sylv)))
!  DO j=1,Nev
!     write(*,*) j,det_ev_Sylv(j)/maxdeteval ! normalized evaluations of resultant
!  ENDDO
! -------------------------------------------------------------------
!  do h=1,22
!     do k=1,22
!        write(*,*) h,k,DBLE(ev_Sylv(k,h,1)) ! the 1st root of unity is real
!     enddo
!  enddo
!  write(*,*)'det(1)=',det_ev_Sylv(1)

  CALL code_input_QP(Nev,det_ev_Sylv,ev_detSylv) 

  ev_detSylv = ev_detSylv/norm_res

  CALL irvfft_QP(ev_detSylv,Nev,expo) 

! polynomial coefficients
  DO ncoe = 0,poldeg
     detS_coe(ncoe) = ev_detSylv(ncoe+1)
  ENDDO
! residual polynomial coefficients (should be zero)
  DO ncoe = poldeg+1,Nev-1
     detS_coe(ncoe) = ev_detSylv(ncoe+1)
  ENDDO

! -------------------------------------------------------------------
!  write(*,*)'CHECK the RESULTANT EVALUATION in the 64th ROOTS of UNITY'
!  DO j=1,Nev
!     deteval(j) = 0.q0
!     DO h=0,48
!        deteval(j) =  detS_coe(h)*compl_evx(j)**h + deteval(j)
!     ENDDO
!  ENDDO
!  maxdeteval=maxval(abs(DBLE(deteval)))
!  DO j=1,Nev
!     write(*,*) j,deteval(j)/maxdeteval
!  ENDDO
! -------------------------------------------------------------------

!! polynomial coefficients
  WRITE(*,*)'Resultant coefficients'
  DO ncoe = 0,poldeg
     write(*,*)ncoe,detS_coe(ncoe)
  ENDDO
!! residual polynomial coefficients (should be zero)
  write(*,*)'Residual coefficients (should be zero)'
  DO ncoe = poldeg+1,Nev-1
     write(*,*)ncoe,detS_coe(ncoe)
  ENDDO

! search for real roots of Res(rho2) := resultant(p,q,rho1)(rho2)
  CALL solve_poly_QP(poldeg,detS_coe(0:poldeg),roots,nroots,hzflag,multfl)
  IF(verb_2body.gt.10)THEN
     write(*,*)'nroots=',nroots
  ENDIF
  IF(nroots.eq.0)THEN
     nprel=nroots
     RETURN
  ENDIF
  DO ncoe = 1,nroots
     rho2sol(ncoe) = roots(ncoe)
  ENDDO

! compute for each rho2 the corresponding value of rho1
  CALL solvesys_QP(pmod,qmod,ncoe_p,ncoe_q,nroots, &
       & rho2sol(1:nroots),rho1pos(1:nroots),rho2pos(1:nroots),nsol)
  IF(verb_2body.gt.10)THEN
     WRITE(*,*)'nsol=',nsol
     WRITE(*,*)'      rho1             rho2'
     DO ncoe=1,nsol
        WRITE(*,111) ncoe,rho1pos(ncoe),rho2pos(ncoe)
     ENDDO
  ENDIF
111 FORMAT(i2,2x,2(f15.10,2x))
  
!  GOTO 313

  IF(verb_2body.gt.10) WRITE(*,*)'applying Newton-Raphson...'
! ***********************************************************************
! ======================= NEWTON - RAPHSON ========================
! improve the results 
  nsol_tmp=0
  DO ncoe=1,nsol
!     write(*,*)'ncoe',ncoe
     rho1store=rho1pos(ncoe)
     rho2store=rho2pos(ncoe)
     DO n=1,nitmax
        rho1=rho1pos(ncoe)
        rho2=rho2pos(ncoe)
        CALL orbprelim_newton_raphson_QP(pmod,qmod,ncoe_p,ncoe_q, &
             & rho1pos(ncoe),rho2pos(ncoe),pk,qk)
        drho1=rho1-rho1pos(ncoe)
        drho2=rho2-rho2pos(ncoe)

        minv=1.d-20*min(rho1,rho2)
        maxv=1.d0*min(rho1,rho2)

        IF(sqrt(drho1**2+drho2**2).lt.minv)THEN
           nsol_tmp=nsol_tmp+1
           rho1_tmp(nsol_tmp)=rho1pos(ncoe)
           rho2_tmp(nsol_tmp)=rho2pos(ncoe)
           WRITE(*,*)'good convergence!'
           WRITE(*,*)'ncoe,niter,diff:',ncoe,n,sqrt(drho1**2+drho2**2)
           WRITE(*,*)'pk,qk',pk,qk
           WRITE(*,*)'RHO1,RHO2 (dopo):',rho1pos(ncoe),rho2pos(ncoe)
           WRITE(*,*)''
           EXIT
        ELSEIF(sqrt(drho1**2+drho2**2).gt.maxv)THEN
           nsol_tmp=nsol_tmp+1
           rho1_tmp(nsol_tmp)=rho1store
           rho2_tmp(nsol_tmp)=rho2store
           WRITE(*,*)'correction too big! keep first step solution'
           WRITE(*,*)'ncoe,niter,diff:',ncoe,n,sqrt(drho1**2+drho2**2)
           WRITE(*,*)'RHO1,RHO2 (dopo):',rho1_tmp(nsol_tmp),rho2_tmp(nsol_tmp)
           EXIT
        ELSEIF(n.eq.nitmax)THEN
           WRITE(*,*)'convergence not attained but the solution is accepted'
           nsol_tmp=nsol_tmp+1
           rho1_tmp(nsol_tmp)=rho1pos(ncoe)
           rho2_tmp(nsol_tmp)=rho2pos(ncoe)
           WRITE(*,*)'ncoe,niter,diff:',ncoe,n,sqrt(drho1**2+drho2**2)
           WRITE(*,*)'pk,qk',pk,qk
           WRITE(*,*)'RHO1,RHO2 (dopo):',rho1pos(ncoe),rho2pos(ncoe)
           WRITE(*,*)''
        ENDIF
     ENDDO
  ENDDO
! rename solutions
  nsol=nsol_tmp
  IF(nsol.eq.0)THEN
     nprel=nsol
     RETURN
  ENDIF
  IF(verb_2body.gt.10)THEN
     WRITE(*,*)'nsol after Newton-Raphson:',nsol
     WRITE(*,*)'      rho1             rho2'
  ENDIF
  DO ncoe=1,nsol 
     rho1pos(ncoe)=rho1_tmp(ncoe)
     rho2pos(ncoe)=rho2_tmp(ncoe)
     IF(verb_2body.gt.10)THEN
        WRITE(*,111) ncoe,rho1pos(ncoe),rho2pos(ncoe)
     ENDIF
  ENDDO
!  WRITE(*,*)''
  IF(PRESENT(nsolstep))THEN
     nsolstep(1) = nsol
  ENDIF

313 CONTINUE

! check for positive rho1,rho2
  IF(verb_2body.gt.10)THEN
     WRITE(*,*) '(*) check for positive rho1,rho2 (>mindist AU)'
  ENDIF
  nsol_tmp=0
  DO ncoe=1,nsol 
     IF(rho1pos(ncoe).gt.mindist.AND.rho2pos(ncoe).gt.mindist)THEN
        nsol_tmp=nsol_tmp+1
        rho1_tmp(nsol_tmp)=rho1pos(ncoe)
        rho2_tmp(nsol_tmp)=rho2pos(ncoe)
     ELSE
        IF(verb_2body.gt.10)THEN
           WRITE(*,*)'very small or negative component: rho1,rho2',&
                & rho1pos(ncoe),rho2pos(ncoe)
        ENDIF
     ENDIF
  ENDDO
! rename solutions
  nsol=nsol_tmp
  IF(nsol.eq.0)THEN
     nprel=nsol
     RETURN
  ENDIF
  IF(verb_2body.gt.10)THEN
     WRITE(*,*)'nsol after removing very small or negative rho:',nsol
     WRITE(*,*)''
  ENDIF
  DO ncoe=1,nsol 
     rho1pos(ncoe)=rho1_tmp(ncoe)
     rho2pos(ncoe)=rho2_tmp(ncoe)
  ENDDO

! check for duplicates
  IF(verb_2body.gt.10)THEN
     WRITE(*,*) '(*) check for duplicate rho1,rho2'
  ENDIF
  nsol_tmp=0
  DO ncoe=1,nsol 
     duplfl=.FALSE.
     DO ncoe1=1,ncoe-1 
        norm2=sqrt((rho1pos(ncoe)-rho1pos(ncoe1))**2 + &
             & (rho2pos(ncoe)-rho2pos(ncoe1))**2)
        IF(norm2.lt.norm2x)THEN
           duplfl=.TRUE.
           EXIT
        ENDIF
     ENDDO
     IF(.NOT.duplfl)THEN
        nsol_tmp=nsol_tmp+1
        rho1_tmp(nsol_tmp)=rho1pos(ncoe)
        rho2_tmp(nsol_tmp)=rho2pos(ncoe)
     ELSE
        IF(verb_2body.gt.10)THEN
           WRITE(*,*)'duplicate solution'
        ENDIF
     ENDIF
  ENDDO
! rename solutions
  nsol=nsol_tmp
  IF(nsol.eq.0)THEN
     nprel=nsol
     RETURN
  ENDIF
  IF(PRESENT(nsolstep))THEN
     nsolstep(2) = nsol
  ENDIF
  IF(verb_2body.gt.10)THEN
     WRITE(*,*)'nsol after removing duplicates:',nsol
     WRITE(*,*)''
  ENDIF
  DO ncoe=1,nsol 
     rho1pos(ncoe)=rho1_tmp(ncoe)
     rho2pos(ncoe)=rho2_tmp(ncoe)
  ENDDO

! ***********************************************************************
!  GOTO 133 ! to avoid checks for spurious solutions
  IF(verb_2body.gt.10)THEN
     WRITE(*,*)'(*) spurious solutions: first check'
  ENDIF
! discard spurious solutions of first kind (S1)
  nsol_tmp=0
  DO ncoe=1,nsol
     CALL poly_eval_QP(fact2mod,coe_fact2,rho1pos(ncoe),rho2pos(ncoe),&
          & fact2eval)  
     IF(fact2eval.lt.1.d2*epsilon(1.q0))THEN
        nsol_tmp=nsol_tmp+1
        rho1_tmp(nsol_tmp)=rho1pos(ncoe)
        rho2_tmp(nsol_tmp)=rho2pos(ncoe)
     ELSEIF(fact2eval.gt.1.d2*epsilon(1.q0))THEN
        IF(verb_2body.gt.10)THEN
           WRITE(*,*)'spurious solution! ncoe,fact2eval=',ncoe,fact2eval
        ENDIF
     ELSE
        WRITE(*,*)'ERROR: cannot decide sign of fact2!!'
!        STOP
     ENDIF
  ENDDO
! rename solutions
  nsol=nsol_tmp
  IF(nsol.eq.0)THEN
     nprel=nsol
     RETURN
  ENDIF
  IF(verb_2body.gt.10)THEN
     WRITE(*,*)'nsol after first check for spurious:',nsol
     WRITE(*,*)''
  ENDIF
  DO ncoe=1,nsol 
     rho1pos(ncoe)=rho1_tmp(ncoe)
     rho2pos(ncoe)=rho2_tmp(ncoe)
  ENDDO

! ----------------------------------------------------------------------
!  GOTO 133 ! to avoid second check for spurious solutions
  IF(verb_2body.gt.10)THEN
     WRITE(*,*)'(*) spurious solutions: second check'
  ENDIF
  nsol_tmp=0
  DO ncoe=1,nsol
     CALL poly_eval_QP(dF1F2mod,coe_diffF1F2,rho1pos(ncoe),rho2pos(ncoe),&
          & dF1F2eval)
     CALL poly_eval_QP(2,coe_G1,rho1pos(ncoe),rho2pos(ncoe),&
          & G1eval)
     CALL poly_eval_QP(2,coe_G2,rho1pos(ncoe),rho2pos(ncoe),&
          & G2eval)
     IF(sign(1.q0,dF1F2eval).eq.sign(1.q0,sqrt(G2eval)-sqrt(G1eval)))THEN
        nsol_tmp=nsol_tmp+1
        rho1_tmp(nsol_tmp)=rho1pos(ncoe)
        rho2_tmp(nsol_tmp)=rho2pos(ncoe)
     ELSE
        IF(verb_2body.gt.10)THEN
           WRITE(*,*)'spurious solution! ncoe,dF1F2eval, &
             & 2k^2*(1/sqrtG2-1/sqrtG1)=',&
             & ncoe,dF1F2eval,2.q0*gmvarc*(1.q0/sqrt(G1eval)-1.q0/sqrt(G2eval))
        ENDIF
        CALL poly_eval_QP(4,coe_F1,rho1pos(ncoe),rho2pos(ncoe),&
          & F1eval)
        CALL poly_eval_QP(4,coe_F2,rho1pos(ncoe),rho2pos(ncoe),&
          & F2eval)
     ENDIF
  ENDDO
! rename solutions
  nsol=nsol_tmp
  IF(nsol.eq.0)THEN
     nprel=nsol
     RETURN
  ENDIF
  IF(verb_2body.gt.10)THEN
     write(*,*)'nsol after second check for spurious:',nsol
     WRITE(*,*)'      rho1             rho2'
     DO ncoe=1,nsol 
        WRITE(*,111) ncoe,rho1_tmp(ncoe),rho2_tmp(ncoe)
     ENDDO
  ENDIF
  DO ncoe=1,nsol 
     rho1pos(ncoe)=rho1_tmp(ncoe)
     rho2pos(ncoe)=rho2_tmp(ncoe)
  ENDDO
!  WRITE(*,*)''
133 CONTINUE
! ----------------------------------------------------------------------

! compute corresponding rho1dot, rho2dot solution
  DO ncoe=1,nsol
     rho1dot(ncoe)=0.q0; rho2dot(ncoe)=0.q0;
     DO h=0,2
        DO k=0,2
           IF(coe_rd1(h,k).ne.0.q0)THEN
              rho1dot(ncoe) = rho1dot(ncoe) + &
                   & coe_rd1(h,k)*rho1pos(ncoe)**h*rho2pos(ncoe)**k
           ENDIF
           IF(coe_rd2(h,k).ne.0.q0)THEN
              rho2dot(ncoe) = rho2dot(ncoe) + &
                   & coe_rd2(h,k)*rho1pos(ncoe)**h*rho2pos(ncoe)**k
           ENDIF
        ENDDO
     ENDDO
     write(*,*)'rho1dot(ncoe),rho2dot(ncoe)', rho1dot(ncoe),rho2dot(ncoe)
  ENDDO

! initialization
  elcar1=undefined_orbit_elem
  elcar2=undefined_orbit_elem

! Cartesian elements
  DO ncoe=1,nsol
     rho1pos(ncoe)=rho1pos(ncoe)/convdist
     rho1dot(ncoe)=rho1dot(ncoe)/convdist
     elcar1(ncoe)=undefined_orbit_elem
     elcar1(ncoe)%coo='CAR'
     IF(rhs.EQ.2)THEN
        elcar1(ncoe)%center=3
     END IF
     elcar1(ncoe)%t=t1-rho1pos(ncoe)/vlight ! aberration correction
     elcar1(ncoe)%coord(1:3)=pe(1:3,1)+rho1pos(ncoe)*rhohat(1:3,1)
     IF(rhs.NE.2)THEN
        elcar1(ncoe)%coord(1:3) = MATMUL(roteqec,elcar1(ncoe)%coord(1:3))
     END IF
     elcar1(ncoe)%coord(4:6)=pde(1:3,1) + rho1dot(ncoe)*rhohat(1:3,1) + &
          & rho1pos(ncoe)*(rhohat_a(1:3,1)*ad(1)+rhohat_d(1:3,1)*dd(1))
     IF(rhs.NE.2)THEN
        elcar1(ncoe)%coord(4:6) = MATMUL(roteqec,elcar1(ncoe)%coord(4:6))
     END IF
     elcar1(ncoe)%obscode=obscode1
!
     rho2pos(ncoe)=rho2pos(ncoe)/convdist
     rho2dot(ncoe)=rho2dot(ncoe)/convdist
     elcar2(ncoe)=undefined_orbit_elem
     elcar2(ncoe)%coo='CAR'
     IF(rhs.EQ.2)THEN
        elcar2(ncoe)%center=3
     ENDIF
     elcar2(ncoe)%t=t2-rho2pos(ncoe)/vlight ! aberration correction
     elcar2(ncoe)%coord(1:3)=pe(1:3,2) + rho2pos(ncoe)*rhohat(1:3,2)
     IF(rhs.NE.2)THEN
        elcar2(ncoe)%coord(1:3)=MATMUL(roteqec,elcar2(ncoe)%coord(1:3))
     END IF
     elcar2(ncoe)%coord(4:6)=pde(1:3,2)+rho2dot(ncoe)*rhohat(1:3,2)+&
          &rho2pos(ncoe)*(rhohat_a(1:3,2)*ad(2)+rhohat_d(1:3,2)*dd(2))
     IF(rhs.NE.2)THEN
        elcar2(ncoe)%coord(4:6) = MATMUL(roteqec,elcar2(ncoe)%coord(4:6))
     END IF
     elcar2(ncoe)%obscode=obscode2

! temp to see energy (change gmvar?)
!     WRITE(*,*)'Energy 1=', vsize(elcar1(ncoe)%coord(4:6))**2.d0/2 - &
!          & gmvar/vsize(elcar1(ncoe)%coord(1:3))
!     WRITE(*,*)'Energy 2=', vsize(elcar2(ncoe)%coord(4:6))**2.d0/2 - &
!          & gmvar/vsize(elcar2(ncoe)%coord(1:3))
  ENDDO
  nprel=nsol

  rrdot1=0.d0
  rrdot2=0.d0
  DO ncoe=1,nprel
     rrdot1(ncoe,1)=rho1pos(ncoe)
     rrdot1(ncoe,2)=rho1dot(ncoe)
     rrdot2(ncoe,1)=rho2pos(ncoe)
     rrdot2(ncoe,2)=rho2dot(ncoe)
  ENDDO
 ! ===============================================================
  
 CLOSE(3)
 CLOSE(2)
 CLOSE(1)
 
END SUBROUTINE orbprelim_2b

! ===============================================================
!Given two attributables, the two approximations for radial distance 
!and velocity, the observer position and velocity, computes the 
!partial derivatives of the vector field Phi.
!Written by Linda Dimare, September 2008.
SUBROUTINE derphi(att1,att2,rrd1,rrd2,qqd1,qqd2,phi_att,phi_r)
  USE attributable
!  USE fund_const
  IMPLICIT NONE  
  TYPE(attrib),INTENT(IN) :: att1,att2
  REAL(KIND=dkind), INTENT(IN), DIMENSION(2) :: rrd1, rrd2 !radial dist/vel
  REAL(KIND=dkind), INTENT(IN), DIMENSION(3,2) :: qqd1,qqd2 !obs pos/vel  
  REAL(KIND=dkind), INTENT(OUT), DIMENSION(4,4) :: phi_r
  REAL(KIND=dkind), INTENT(OUT), DIMENSION(4,8) :: phi_att
  ! end interface

  REAL(KIND=dkind), DIMENSION(2) :: rho, rhod ! radial dist.& vel.
  REAL(KIND=dkind), DIMENSION(2) :: alpha,delta,ad,dd ! attr.
  REAL(KIND=dkind), DIMENSION(3,2) :: q,qdot ! obs.position & vel.
  REAL(KIND=dkind), DIMENSION(4,4) :: phi_r_rdot    
  REAL(KIND=dkind), DIMENSION(3,2) :: rhohat_aa, rhohat_dd, rhohat_ad
  REAL(KIND=dkind), DIMENSION(3,4) :: D1_A1,E1_A1,F1_A1
  REAL(KIND=dkind), DIMENSION(3,4) :: D2_A2,E2_A2,F2_A2
  
  REAL(KIND=dkind), DIMENSION(3,2) :: q_rhata, q_rhatd
  REAL(KIND=dkind), DIMENSION(3,2) :: rhat_rhataa, rhat_rhatdd, rhat_rhatad
  REAL(KIND=dkind), DIMENSION(3,2) :: rhat_rhata, rhat_rhatd,rhata_rhatd
  REAL(KIND=dkind), DIMENSION(3,2) :: q_rhataa, q_rhatdd, q_rhatad
  REAL(KIND=dkind), DIMENSION(3,2) :: rhata_qdot, rhatd_qdot

  REAL(KIND=dkind), DIMENSION(4,2) :: c1_att, c2_att, c3_att, c5_att
  
  INTEGER :: i,j !loop parameters
 ! =========================================================== 

  alpha(1) = att1%angles(1)
  delta(1) = att1%angles(2)
  ad(1) = att1%angles(3)
  dd(1) = att1%angles(4)
  alpha(2) = att2%angles(1)
  delta(2) = att2%angles(2)
  ad(2) = att2%angles(3)
  dd(2) = att2%angles(4)

  rho(1) = rrd1(1)
  rhod(1) = rrd1(2)
  rho(2) = rrd2(1)
  rhod(2) = rrd2(2)

  q(1:3,1) = qqd1(1:3,1)
  qdot(1:3,1)=qqd1(1:3,2)
  q(1:3,2) = qqd2(1:3,1)
  qdot(1:3,2)=qqd2(1:3,2)

  phi_r_rdot=0.d0
  phi_att=0.d0  
  
!derivatives with respect to rho,rhodot
  DO j=1,2
     DO i=1,3
        phi_r_rdot(i,j)=(-1.d0)**(j+1)*(2.d0*Evec(i,j)*rho(j)+Fvec(i,j))
        phi_r_rdot(i,j+2)=(-1.d0)**(j+1)*Dvec(i,j)
     ENDDO
     phi_r_rdot(4,j)=(-1.d0)**(j+1)*(2.d0*c2(j)*rho(j)+c3(j)+ & 
          & gmvar*(rho(j)**2+c5(j)*rho(j)+c0(j))**(-1.5d0)*(2.d0*rho(j)+c5(j)))
     phi_r_rdot(4,j+2)=(-1.d0)**(j+1)*(2.d0*rhod(j)+c1(j)) 
  ENDDO
  
!derivatives with respect to attributables
  D1_A1=0.d0
  D2_A2=0.d0
  DO j=1,2 
     CALL cross_prod(q(1:3,j),rhohat_a(1:3,j),q_rhata(1:3,j))
     CALL cross_prod(q(1:3,j),rhohat_d(1:3,j),q_rhatd(1:3,j))
  ENDDO
 
  D1_A1(1:3,1)=q_rhata(1:3,1)
  D1_A1(1:3,2)=q_rhatd(1:3,1)
  D2_A2(1:3,1)=q_rhata(1:3,2)
  D2_A2(1:3,2)=q_rhatd(1:3,2) 

!second derivatives of rhohat
  DO j=1,2
     rhohat_aa(1,j)= -cosa(j)*cosd(j)
     rhohat_aa(2,j)= -sina(j)*cosd(j)
     rhohat_aa(3,j)= 0.d0

     rhohat_dd(1,j)= -cosa(j)*cosd(j)
     rhohat_dd(2,j)= -sina(j)*cosd(j)
     rhohat_dd(3,j)= -sind(j)

     rhohat_ad(1,j)= sina(j)*sind(j)
     rhohat_ad(2,j)= -cosa(j)*sind(j)
     rhohat_ad(3,j)= 0.d0
  ENDDO

  DO j=1,2 
     CALL cross_prod(rhohat(1:3,j),rhohat_aa(1:3,j),rhat_rhataa(1:3,j))
     CALL cross_prod(rhohat(1:3,j),rhohat_dd(1:3,j),rhat_rhatdd(1:3,j))
     CALL cross_prod(rhohat(1:3,j),rhohat_ad(1:3,j),rhat_rhatad(1:3,j))
     CALL cross_prod(rhohat_a(1:3,j),rhohat_d(1:3,j),rhata_rhatd(1:3,j))
     CALL cross_prod(rhohat(1:3,j),rhohat_a(1:3,j),rhat_rhata(1:3,j))
     CALL cross_prod(rhohat(1:3,j),rhohat_d(1:3,j),rhat_rhatd(1:3,j))
     CALL cross_prod(q(1:3,j),rhohat_aa(1:3,j),q_rhataa(1:3,j))
     CALL cross_prod(q(1:3,j),rhohat_dd(1:3,j),q_rhatdd(1:3,j))
     CALL cross_prod(q(1:3,j),rhohat_ad(1:3,j),q_rhatad(1:3,j))
     CALL cross_prod(rhohat_a(1:3,j),qdot(1:3,j),rhata_qdot(1:3,j))
     CALL cross_prod(rhohat_d(1:3,j),qdot(1:3,j),rhatd_qdot(1:3,j))
  ENDDO
 
  E1_A1(1:3,1)= ad(1)*rhat_rhataa(1:3,1)+dd(1)*(rhata_rhatd(1:3,1)+ &
       & rhat_rhatad(1:3,1))
  E1_A1(1:3,2)= ad(1)*(-rhata_rhatd(1:3,1)+rhat_rhatad(1:3,1))+ &
       & dd(1)*rhat_rhatdd(1:3,1)
  E1_A1(1:3,3)= rhat_rhata(1:3,1)
  E1_A1(1:3,4)= rhat_rhatd(1:3,1)

  E2_A2(1:3,1)= ad(2)*rhat_rhataa(1:3,2)+dd(2)*(rhata_rhatd(1:3,2)+ &
       & rhat_rhatad(1:3,2))
  E2_A2(1:3,2)= ad(2)*(-rhata_rhatd(1:3,2)+rhat_rhatad(1:3,2))+ &
       & dd(2)*rhat_rhatdd(1:3,2)
  E2_A2(1:3,3)= rhat_rhata(1:3,2)
  E2_A2(1:3,4)= rhat_rhatd(1:3,2)

  F1_A1(1:3,1)= ad(1)*q_rhataa(1:3,1)+dd(1)*q_rhatad(1:3,1)+rhata_qdot(1:3,1)
  F1_A1(1:3,2)= ad(1)*q_rhatad(1:3,1)+dd(1)*q_rhatdd(1:3,1)+rhatd_qdot(1:3,1)
  F1_A1(1:3,3)= q_rhata(1:3,1)
  F1_A1(1:3,4)= q_rhatd(1:3,1)

  F2_A2(1:3,1)= ad(2)*q_rhataa(1:3,2)+dd(2)*q_rhatad(1:3,2)+rhata_qdot(1:3,2)
  F2_A2(1:3,2)= ad(2)*q_rhatad(1:3,2)+dd(2)*q_rhatdd(1:3,2)+rhatd_qdot(1:3,2)
  F2_A2(1:3,3)= q_rhata(1:3,2)
  F2_A2(1:3,4)= q_rhatd(1:3,2)

  DO i=1,4
     phi_att(1:3,i)= rhod(1)*D1_A1(1:3,i)+rho(1)**2*E1_A1(1:3,i)+ &
          & rho(1)*F1_A1(1:3,i)
     phi_att(1:3,i+4)= -rhod(2)*D2_A2(1:3,i)-rho(2)**2*E2_A2(1:3,i)- &
          & rho(2)*F2_A2(1:3,i)
  ENDDO
 
  c1_att=0.d0
  c2_att=0.d0
  c3_att=0.d0
  c5_att=0.d0

  DO j=1,2
     c1_att(1,j)= 2.d0*DOT_PRODUCT(qdot(1:3,j),rhohat_a(1:3,j))
     c1_att(2,j)= 2.d0*DOT_PRODUCT(qdot(1:3,j),rhohat_d(1:3,j))
     c2_att(2,j)= -2.d0*ad(j)**2*cos(delta(j))*sin(delta(j))
     c2_att(3,j)= 2.d0*ad(j)*(cos(delta(j)))**2
     c2_att(4,j)= 2.d0*dd(j)
     c3_att(1,j)= 2.d0*(ad(j)*DOT_PRODUCT(qdot(1:3,j),rhohat_aa(1:3,j))+ &
          & dd(j)*DOT_PRODUCT(qdot(1:3,j),rhohat_ad(1:3,j)))
     c3_att(2,j)= 2.d0*(ad(j)*DOT_PRODUCT(qdot(1:3,j),rhohat_ad(1:3,j))+ &
          & dd(j)*DOT_PRODUCT(qdot(1:3,j),rhohat_dd(1:3,j)))
     c3_att(3,j)= 2.d0*DOT_PRODUCT(qdot(1:3,j),rhohat_a(1:3,j)) 
     c3_att(4,j)= 2.d0*DOT_PRODUCT(qdot(1:3,j),rhohat_d(1:3,j))
     c5_att(1,j)= 2.d0*DOT_PRODUCT(q(1:3,j),rhohat_a(1:3,j))
     c5_att(2,j)= 2.d0*DOT_PRODUCT(q(1:3,j),rhohat_d(1:3,j))
  ENDDO

  DO i=1,4
     phi_att(4,i)= rhod(1)*c1_att(i,1)+(rho(1)**2)*c2_att(i,1) + &
          & rho(1)*c3_att(i,1)+ &
          & gmvar*((rho(1)**2+rho(1)*c5(1)+c0(1))**(-1.5d0))*(rho(1)*c5_att(i,1))
     phi_att(4,i+4)= -(rhod(2)*c1_att(i,2)+(rho(2)**2)*c2_att(i,2)+ &
          & rho(2)*c3_att(i,2)+&
          & gmvar*((rho(2)**2+rho(2)*c5(2)+c0(2))**(-1.5d0))*(rho(2)*c5_att(i,2)))
  ENDDO

  phi_r=phi_r_rdot  
  phi_r(1:4,2)=phi_r_rdot(1:4,3)
  phi_r(1:4,3)=phi_r_rdot(1:4,2)

END SUBROUTINE derphi

END MODULE kepint

! ==============================================================
! subroutine matin_qp

! INPUT:    A(i,j)    -   Matrix of coefficients of the system          
!           N         -   Order of A(i,j)                               
!           L         -   Number of linear systems to be solved (the    
!                         right hand  sides are stored in A(i,N+j),     
!                         j=1,L)                                        
!           NL        -   First dimension of A(i,j)                     
!           INVOP     -   If =1 the inverse of A(i,j) is computed       
!                         explicitly,                                   
!                         if =0 only the linear systems are solved      
!                                                                       
! OUTPUT:   A(i,j)    -   If INVOP=1, contains the inverse of the input 
!                         if INVOP=0, contains the triangular           
!                         factorization of the input matrix.            
!                         In both cases A(i,N+j), j=1,L contain the     
!                         solutions of the input systems.               
!           DET       -   Determinant of the input matrix A(i,j)        
!           ISING     -   If =0 the algorithm was completed             
!                         successfully,                                 
!                         if =1 the input matrix was singular           
!        
subroutine matin_qp(a,det,n,l,nl,ising,invop) 
  USE fund_const, ONLY: qkind
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n,l,nl,invop
  INTEGER, INTENT(OUT) :: ising
  REAL(KIND=qkind), INTENT(INOUT) ::  a(nl,n+l)
  REAL(KIND=qkind), INTENT(OUT) :: det
  INTEGER, PARAMETER :: nmax=1000 
  REAL(KIND=qkind):: vet(nmax) ! workspace 
  INTEGER  p(nmax),j,i,k,nc,npiv,nsc
  REAL(KIND=qkind) amax,sc
if(n.gt.nmax)stop ' **** matin: n > nmax ****' 
   do  j=1,n 
      p(j)=j 
   enddo
   do 4 j=1,n-1 
      amax=0.d0 
      DO i=j,n 
         if(amax.ge.qabs(a(i,j))) CYCLE
         npiv=i 
         amax=qabs(a(i,j)) 
      ENDDO
! modification 2006: error flag must be set
      IF(amax.eq.0.d0)THEN
         ising=1
         RETURN 
      ENDIF
! end modification
      if(npiv.eq.j)goto 6 
      nsc=p(npiv) 
      p(npiv)=p(j) 
      p(j)=nsc 
      do  i=1,n+l 
         sc=a(npiv,i) 
         a(npiv,i)=a(j,i) 
         a(j,i)=sc
      ENDDO 
    6 DO i=j+1,n 
         a(i,j)=a(i,j)/a(j,j)
      ENDDO 
      do 4 i=j+1,n 
      do 4 k=j+1,n+l 
4  a(i,k)=a(i,k)-a(i,j)*a(j,k) 
   if(l.eq.0)goto 7 
   do 9 i=1,l 
   do 9 j=n,1,-1 
      if(j.eq.n)goto 9 
      do  k=j+1,n 
         a(j,n+i)=a(j,n+i)-a(j,k)*a(k,n+i)
      enddo 
9  a(j,n+i)=a(j,n+i)/a(j,j) 
7  det=1.d0 
   do j=1,n 
      det=det*a(j,j) 
   enddo
   ising=0 
   if(invop.eq.0)return 
   if(n.ne.1)goto 20 
   a(1,1)=1.d0/a(1,1) 
   return 
20 do 12 i=2,n 
   do 12 j=1,i-1 
      a(i,j)=-a(i,j) 
      if(j+2.gt.i)goto 12 
      do  k=j+1,i-1 
         a(i,j)=a(i,j)-a(i,k)*a(k,j)
      enddo 
12 continue 
   do 15 k=n,1,-1 
      do 14 i=1,n 
         vet(i)=0.d0 
         if(i.eq.k)vet(i)=1.d0 
         if(k.gt.i)vet(i)=a(k,i) 
         if(k.eq.n)goto 14 
         do  j=k+1,n 
            vet(i)=vet(i)-a(k,j)*a(j,i)
         enddo
   14 continue 
      sc=a(k,k) 
      do 15  i=1,n 
         a(k,i)=vet(i)/sc
15 continue
18 nc=0 
   do 16 j=1,n 
      if(j.eq.p(j))goto 16 
      nsc=p(j) 
      p(j)=p(nsc) 
      p(nsc)=nsc 
      do 17 i=1,n 
         sc=a(i,nsc) 
         a(i,nsc)=a(i,j) 
17       a(i,j)=sc 
         nc=1 
16 continue 
   if(nc.eq.1)goto 18 
   return 
END subroutine matin_qp           

! ===============================================================
SUBROUTINE cross_prod_QP(v1,v2,v3)
  USE fund_const, ONLY: qkind
  IMPLICIT NONE
  REAL(KIND=qkind),INTENT(IN),DIMENSION(3) :: v1,v2
  REAL(KIND=qkind),INTENT(OUT),DIMENSION(3) :: v3
  v3(1) = v1(2)*v2(3)-v1(3)*v2(2)
  v3(2) = v1(3)*v2(1)-v1(1)*v2(3)
  v3(3) = v1(1)*v2(2)-v1(2)*v2(1)
END SUBROUTINE cross_prod_QP

! ===============================================================
SUBROUTINE cross_prod(v1,v2,v3)
  USE fund_const, ONLY: dkind
  IMPLICIT NONE
  REAL(KIND=dkind),INTENT(IN),DIMENSION(3) :: v1,v2
  REAL(KIND=dkind),INTENT(OUT),DIMENSION(3) :: v3
  v3(1) = v1(2)*v2(3)-v1(3)*v2(2)
  v3(2) = v1(3)*v2(1)-v1(1)*v2(3)
  v3(3) = v1(1)*v2(2)-v1(2)*v2(1)
END SUBROUTINE cross_prod

! ===============================================================
SUBROUTINE matrix_mult(A,B,p,q,r,C)
  USE fund_const, ONLY: dkind
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: p,q,r
  REAL(KIND=dkind),INTENT(IN):: A(p,q),B(q,r)
  REAL(KIND=dkind),INTENT(OUT):: C(p,r)
! end interface
  INTEGER :: i,j,k  
  C=0.d0
  DO i=1,p
     DO j=1,r
        DO k=1,q
           C(i,j)=C(i,j)+A(i,k)*B(k,j)
        ENDDO
     ENDDO
  ENDDO
  
END SUBROUTINE matrix_mult


! ========================================================
! **************** N O T   U S E D  **********************

!! Earth elements at time t1
!  elem1=undefined_orbit_elem
!  CALL earth(t1,eqtmp1)
!  elem1%coo='EQU'
!  elem1%t=t1
!  elem1%coord(1:6)=eqtmp1(1:6)
!  call coo_cha(elem1,'CAR',ecpl1,fail_flag)
!!  write(*,*)'pos1 of the Earth center',ecpl1%coord(1:3)
!!  write(*,*)'vel1 of the Earth center',ecpl1%coord(4:6)
!  CALL  observer_position(t1,pos1,vel1,obscode1)
!!  write(*,*)'Observer position from Earth center',pos1
!!  write(*,*)'Observer velocity from Earth center',vel1
!  pe(1:3,1) = MATMUL(roteceq,ecpl1%coord(1:3)+pos1(1:3))
!  pde(1:3,1) = MATMUL(roteceq,ecpl1%coord(4:6)+vel1(1:3))
!  write(*,*)'Ecliptic observer position at t1:',ecpl1%coord(1:3)+pos1(1:3)

!! Earth elements at time t2
!  elem2=undefined_orbit_elem
!  CALL earth(t2,eqtmp2)
!  elem2%coo='EQU'
!  elem2%t=t2
!  elem2%coord(1:6)=eqtmp2(1:6)
!  call coo_cha(elem2,'CAR',ecpl2,fail_flag)
!!  write(*,*)'pos2 of the Earth center',ecpl2%coord(1:3)
!!  write(*,*)'vel2 of the Earth center',ecpl2%coord(4:6)
!  CALL  observer_position(t2,pos2,vel2,obscode2)
!!  write(*,*)'Observer velocity from Earth center',vel2
!  pe(1:3,2) = MATMUL(roteceq,ecpl2%coord(1:3)+pos2(1:3))
!  pde(1:3,2) = MATMUL(roteceq,ecpl2%coord(4:6)+vel2(1:3))

!  IF(PRESENT(qobs1).and.PRESENT(qpobs1).and. &
!       & PRESENT(qobs2).and.PRESENT(qpobs2))THEN
!     write(*,*)'interpolated pos. and .vel. of observer 1:'
!     write(*,*)qobs1,qpobs1
!     write(*,*)'interpolated pos. and .vel. of observer 2:'
!     write(*,*)qobs2,qpobs2
!  ENDIF

! ===============================================================
!SUBROUTINE orbdist(elem1,elem2,orb_dist)
!  USE orbit_elements
!  USE fund_const
!  IMPLICIT NONE
!  TYPE(orbit_elem),INTENT(IN) :: elem1,elem2 !angles in rad
!  REAL(KIND=dkind),INTENT(OUT) :: orb_dist
!!! end interface
!  REAL(KIND=dkind),DIMENSION(6) :: ek1,ek2
!  REAL(KIND=dkind) :: t1,t2,deltat,enne2
!  REAL(KIND=dkind) :: od2,pridif,diff,ell2_to_t1
!  REAL(KIND=dkind) :: w1,w2 ! weights
!  ek1(1:6) = elem1%coord(1:6)
!  ek2(1:6) = elem2%coord(1:6)
!  t1 =  elem1%t
!  t2 =  elem2%t
!  od2 = ((ek1(1)-ek2(1))/ek1(1))**2 + (ek1(2)-ek2(2))**2 
!  diff = pridif(ek1(3),ek2(3)) + pridif(ek1(4),ek2(4))
!  od2 = od2 + diff  
!  deltat=t2-t1
!  enne2 = ek2(1)**(-1.5d0)*sqrt(gmvar)
!  ell2_to_t1 =ek2(6)- enne2*deltat
!  w1 = 0.2d0
!  w2 = 0.1d0
!  diff = w1*pridif(ek1(5),ek2(5))**2 + w2*pridif(ek1(6),ell2_to_t1)**2
!  od2 = od2 + diff  
!  orb_dist=sqrt(od2)
!END SUBROUTINE orbdist


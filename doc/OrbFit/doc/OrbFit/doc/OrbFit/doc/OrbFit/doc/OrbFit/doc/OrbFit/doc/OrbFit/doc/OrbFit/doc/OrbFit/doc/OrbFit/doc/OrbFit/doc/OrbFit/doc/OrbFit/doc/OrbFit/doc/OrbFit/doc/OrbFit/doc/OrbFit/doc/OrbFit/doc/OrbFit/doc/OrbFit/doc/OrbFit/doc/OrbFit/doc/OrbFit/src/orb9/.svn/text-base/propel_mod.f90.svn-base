MODULE propel_mod
  USE name_rules, ONLY: name_len
  IMPLICIT NONE
  PRIVATE

! proper elements, synthetic version
TYPE proper_els
  CHARACTER*(name_len) :: name
  DOUBLE PRECISION               :: hmag  ! absol. mag
  DOUBLE PRECISION, DIMENSION(3) :: pr_el ! a(AU), e, sin I
  DOUBLE PRECISION, DIMENSION(3) :: pr_fr ! n, g, s in rad/y
  DOUBLE PRECISION               :: lce   ! max. Lyap exp. (1/y)
  DOUBLE PRECISION               :: t_int ! integration time span (y) 
END TYPE proper_els

TYPE close_couple
  CHARACTER*(name_len), DIMENSION(2) :: names
  INTEGER, DIMENSION(2)              :: addr
  DOUBLE PRECISION, DIMENSION(3)     :: delta_el ! da/a, de, dsinI
  DOUBLE PRECISION                   :: sigma
END TYPE close_couple


PUBLIC proper_els, close_couple

! public routines
PUBLIC input_propels, taxostep  

! shared data
INTEGER, PARAMETER :: nprox=300000 ! max size of array
TYPE(proper_els), DIMENSION(nprox) :: propel ! storage array
INTEGER    :: npro ! actual number
INTEGER, DIMENSION(nprox) :: isrti ! sort by sine of inclination
INTEGER, PARAMETER :: ncoupx=300000 ! max size of array
TYPE(close_couple), DIMENSION(ncoupx) :: clos 
INTEGER, DIMENSION(ncoupx) :: isrts ! sort by sigma

PUBLIC nprox, propel, npro, isrti, ncoupx, clos, isrts
 
CONTAINS

! input of synthetic proper element
 SUBROUTINE input_propels(iun,npr,err_line)
  INTEGER, INTENT(IN)  :: iun ! input unit
  INTEGER, INTENT(INOUT)  :: npr ! previously assigned in input
  INTEGER, INTENT(OUT) :: err_line ! error flag, if >0 is error location
! 
  CHARACTER*256 record
  INTEGER j
  CHARACTER*(name_len) nam
  DOUBLE PRECISION hmag, pr_el(3),pr_fr(3),lce,t_int
! skip header
  READ(iun,*)
  READ(iun,*)
! initialization
  npro=npr
  err_line=0
! main loop
  DO j=npr+1,nprox
    READ(iun,'(A)',END=3) record
    READ(record,*,ERR=2)nam, hmag,pr_el, pr_fr, lce, t_int
    npro=npro+1
    propel(npro)%name=nam
    propel(npro)%hmag=hmag
    propel(npro)%pr_el=pr_el
    propel(npro)%pr_fr=pr_fr
    propel(npro)%lce=lce
    propel(npro)%t_int=t_int
  ENDDO 
! too low nprox
  WRITE(*,*)' increase nprox, was ',nprox
  npr=npro 
  RETURN
2 WRITE(*,*)' read error in proper elements file, line ',j
  err_line=j
  npr=npro  ! stored anyway
  RETURN
! regular ending at eof
3 WRITE(*,*)' read ',npro,' proper elements'
  npr=npro
 END SUBROUTINE input_propels


  SUBROUTINE taxostep(ii,pro,nam0,hmag,d_max,iun_out,n1)
    INTEGER, INTENT(IN)          :: ii     ! address in array of pro
    DOUBLE PRECISION, INTENT(IN) :: pro(3), hmag  ! current proper elements, magnitude
    CHARACTER*(*), INTENT(IN)    :: nam0    ! orbit name
    INTEGER, INTENT(IN)          :: iun_out ! output unit
    DOUBLE PRECISION, INTENT(IN) :: d_max   ! control 
    INTEGER, INTENT(INOUT)         :: n1      ! number of filter 1 ids
! ==============end interface============================
    DOUBLE PRECISION, PARAMETER :: wa=1.d0, we=1.d0, wi=1.d0 !metrics
    INTEGER indal
    DOUBLE PRECISION d2_max,vmin, vmax
    d2_max=d_max*d_max
    vmax=pro(3)+d_max/wi
    vmin=pro(3)-d_max/wi
    CALL bin_search_rea(pro(3),npro,1,npro,propel(1:npro)%pr_el(3),isrti,indal)
    CALL scan_list(.true.)
    CALL scan_list(.false.)
  CONTAINS
! =========================================                             
!  SCAN_ALPHALIST
! ========================================= 
  SUBROUTINE scan_list(forward)
    LOGICAl, INTENT(IN) :: forward
! end interface
    INTEGER lmin, lmax, lstep, na, ni,i,nn1 ! loop indexes and control
    DOUBLE PRECISION di,de,da, sigma2,sigma
    CHARACTER*(name_len) nam
    INTEGER le, ll ! character manipulations
! find extremes of loop
    IF(forward)THEN
       lmin=indal+1
       lmax=npro
       lstep=1
    ELSE
       lmin=indal
       lmax=1
       lstep=-1
    ENDIF
! loop on list sorted by alpha
    DO 22 na=lmin, lmax, lstep ! na is the index of the list sorted by alpha
       ni=isrti(na) ! ni is the index of the list sorted by sinI
! remove duplication, keep with the first earlier in the original file
       IF(ni.le.ii) CYCLE
! test for distance in sinI
          di=propel(ni)%pr_el(3)
       IF(forward)THEN
          IF(di.gt.vmax) RETURN
       ELSE
          IF(di.lt.vmin) RETURN
       ENDIF
       di=di-pro(3)
       de=propel(ni)%pr_el(2)-pro(2)
       IF(ABS(de*we).gt.d_max) CYCLE
       da=(propel(ni)%pr_el(1)-pro(1))/pro(1)
       IF(ABS(da*wa).gt.d_max) CYCLE
       sigma2=(da*wa)**2+(de*we)**2+(di*wi)**2
       IF(sigma2.gt.d2_max) CYCLE
       sigma=sqrt(sigma2) 
! check for self identification
       IF(ii.eq.ni) CYCLE
! store for possible relationship
       IF(n1.eq.ncoupx)THEN
          WRITE(*,*)' increase ncoupx, was ',ncoupx
       ELSE
          n1=n1+1
          clos(n1)%names(1)=nam0
          clos(n1)%names(2)=propel(ni)%name
          clos(n1)%addr(1)=ii
          clos(n1)%addr(2)=ni
          clos(n1)%sigma=sigma
          clos(n1)%delta_el(1)=da
          clos(n1)%delta_el(2)=de
          clos(n1)%delta_el(3)=di
      ENDIF
22  ENDDO
  END SUBROUTINE scan_list
END SUBROUTINE taxostep



END MODULE propel_mod

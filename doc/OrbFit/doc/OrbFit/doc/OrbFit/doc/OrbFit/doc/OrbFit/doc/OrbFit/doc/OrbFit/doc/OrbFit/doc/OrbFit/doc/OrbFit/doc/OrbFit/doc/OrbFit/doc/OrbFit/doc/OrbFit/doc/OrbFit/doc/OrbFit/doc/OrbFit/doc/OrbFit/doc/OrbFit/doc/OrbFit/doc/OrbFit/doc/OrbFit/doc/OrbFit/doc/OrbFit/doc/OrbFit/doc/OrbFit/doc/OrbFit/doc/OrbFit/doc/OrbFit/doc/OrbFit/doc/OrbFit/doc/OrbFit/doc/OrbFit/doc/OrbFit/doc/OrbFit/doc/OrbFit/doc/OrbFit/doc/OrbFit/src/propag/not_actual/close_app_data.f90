MODULE close_app_data
USE orbit_elements
PUBLIC
! common data: former covariance.h
!      DOUBLE PRECISION gc_store(6,6) ! covariance matrix, to be used by strclan
      TYPE(orb_uncert) unc_store
      DOUBLE PRECISION coord_store(6)
      CHARACTER*3 coo_store
      LOGICAL covava ! availability flag
! common data: former cloapp.h
! Close approaches parameters
! Vers. 2.2.8, 5 Oct. 2001; A. Milani
! storage space for multiple minima/multiple target plane crossing
      INTEGER njcx,njc
      PARAMETER (njcx=50)
! planet approached previously
      INTEGER ipla0
! planetary coordinates 
      DOUBLE PRECISION xplaj(6,njcx)
! currently (last) approached planet, selected minimum 
      INTEGER iplam,jcsel
! LOCAL MOID at beginnign of encounter, angle of planetary position w.r. to MOID point
      DOUBLE PRECISION moid0,angmoid
! control on angular distance from MOID point
      INTEGER angx
      PARAMETER (angx=30)
! close approach time, relative position and velocity, with 
! partial derivatives with respect to cartesian coord, min. distance
      DOUBLE PRECISION tcla(njcx),xcla(21,njcx),vcla(21,njcx),rmin(njcx)

! common data; former tpcana.h
! reference system adapted to the MTP
      DOUBLE PRECISION vt3(3,3,njcx),v3(3,3,njcx)
! output from mtprot
      DOUBLE PRECISION tpc(3,njcx),dtpdet(6,3,njcx),sig(2,njcx),    &
     &         axes(2,2,njcx),tpr(njcx),svv(njcx),cxv(njcx),czv(njcx)
! trace of LOv on TP, time of close app.
      DOUBLE PRECISION wtp(2,njcx),wtpr(njcx),wtpal(njcx),wtpv(njcx)


END MODULE close_app_data

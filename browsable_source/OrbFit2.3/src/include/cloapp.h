* Close approaches parameters
c Vers. 2.2.8, 5 Oct. 2001; A. Milani
c storage space for multiple minima/multiple target plane crossing
      INTEGER njcx,njc
      PARAMETER (njcx=20)
c planet approached previously
      INTEGER ipla0
c planetary coordinates 
      DOUBLE PRECISION xplaj(6,njcx)
c currently (last) approached planet, selected minimum 
      INTEGER iplam,jcsel
c LOCAL MOID at beginnign of encounter, angle of planetary position w.r. to MOID point
      DOUBLE PRECISION moid0,angmoid
c control on angular distance from MOID point
      INTEGER angx
      PARAMETER (angx=3.d1)
c close approach time, relative position and velocity, with 
c partial derivatives with respect to cartesian coord, min. distance
      DOUBLE PRECISION tcla(njcx),xcla(21,njcx),vcla(21,njcx),rmin(njcx)
c common for al that
      COMMON/clos7/tcla,xcla,vcla,xplaj,rmin,moid0,angmoid,njc,
     +   iplam,ipla0,jcsel


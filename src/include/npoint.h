c common to handle data on multiple observations
      integer npoinx
      parameter(npoinx=4000)
c phase, distance to Earth, distance to Sun (to compute magnitude)
c elongation, galactic latitude, apparent motion
      double precision phav(npoinx),disv(npoinx),dsunv(npoinx)
      double precision elov(npoinx)
      double precision gallav(npoinx),adotv(npoinx),ddotv(npoinx)
      common/npoint/phav,disv,dsunv,elov,gallav,adotv,ddotv

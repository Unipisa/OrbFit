c common containing all physical information
c on the current asteroid needed to compute Yarkovsky acceleration
      DOUBLE PRECISION yarkp(9),beya(7),alya(7)
      DOUBLE PRECISION spya,sqya,etaya75,fmeaya,thfacya,radfluya
c logical falgs: availability of physical data, 
c  has yarkinit routine been called
      LOGICAL yarfil,yarini
      COMMON/yarko1/yarkp,alya,beya,yarfil,yarini
      COMMON/yarko2/spya,sqya,etaya75,fmeaya,thfacya,radfluya

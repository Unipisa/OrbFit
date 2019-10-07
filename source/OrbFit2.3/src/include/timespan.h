c extremes of time of jpleph
c units: MJD
      double precision tejpl1,tejpl2
      common/estjpl/tejpl1,tejpl2
c boundaries of timespan where ET-UT data are available
      double precision temut1,temut2
      common/estut1/temut1,temut2
c allow extrapolation beyond boundaries
      logical temute
      common/estut2/temute

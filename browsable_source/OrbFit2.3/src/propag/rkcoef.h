c ======================================
c  coefficients RKG for rkimp
      double precision a,b,c
      common/rkcoef/a(ismax,ismax),b(ismax),c(ismax)
c  coefficients interpolator for k array, for kintrp
      double precision a1
      common/rkint/a1(ismax,ismax)

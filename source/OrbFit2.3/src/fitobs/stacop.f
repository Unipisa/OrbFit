c =====================================================================
c STACOP
c =====================================================================
c copy state variables to current state
c ================INTERFACE=========================================
      subroutine stacop(icop,eq0,g0,c0,csino0,delno0,
     +         eqc,gc,cc,csinoc,delnoc)
      implicit none
      integer icop
      double precision eq0(6),g0(6,6),c0(6,6),csino0,delno0
      double precision eqc(6),gc(6,6),cc(6,6),csinoc,delnoc
c =================END INTERFACE====================================
c function: copy to/from
      if(icop.eq.1)then
         call vcopy(6,eq0,eqc)
         call mcopy(6,6,g0,gc)
         call mcopy(6,6,c0,cc)
         csinoc=csino0
         delnoc=delno0
      elseif(icop.eq.2)then
         call vcopy(6,eqc,eq0)
         call mcopy(6,6,gc,g0)
         call mcopy(6,6,cc,c0)
         csino0=csinoc
         delno0=delnoc
      else
         write(*,*)' stacop: copying not known, icop=',icop
         stop
      endif
      return
      end

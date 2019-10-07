MODULE offlov_checktp
  USE fund_const
  IMPLICIT NONE

PRIVATE
! common data  
  DOUBLE PRECISION, PUBLIC :: bsdmin,bsdmax  
! public routines
  PUBLIC riskchecktp, header_risk

CONTAINS
!                 newton_checktp
!  	          riskchecktp 
!                 big_vitp (NOT IN USE)
! =======================================================
! NEWTON_CHECK                                                          
! called by riskcheck                                                   
! target plane Newton's method to find minimum distance                 
! off the LOV
SUBROUTINE newton_checktp(iunnew,iunwarn,va_tracemin,b_e,distmin,va_tracenew,fold)
  USE obssto
  USE tp_trace
  USE multiple_sol
  USE least_squares  
  USE orbit_elements   
  USE propag_state 
  USE close_app, ONLY: kill_propag                        
  USE virtual_impactor
!===========================INPUT=======================================
  INTEGER, INTENT(IN)           :: iunnew,iunwarn  ! units for output
  TYPE(tp_point), INTENT(IN)    :: va_tracemin 
  DOUBLE PRECISION, INTENT(IN)  :: b_e             ! Earth radius increased 
                                                       ! by gravitational focusing  
!========================== OUTPUT======================================
  DOUBLE PRECISION, INTENT(OUT) :: distmin
  TYPE(tp_point), INTENT(OUT)   :: va_tracenew
  LOGICAL, INTENT(OUT)          :: fold 
! ======================END INTERFACE==================================                 
  TYPE(orbit_elem) :: el0          ! elements     
  TYPE(orb_uncert) :: unc0  ! normal and covariance matrices
  LOGICAL                           :: falsok 
!=============target plane before the current iteration==================
  DOUBLE PRECISION :: r,sina,cosa,str,wid,csi,zeta,dq 
  DOUBLE PRECISION, DIMENSION(2) :: dtpc
!==========correspondent ellipse in orbital elements space=============== 
  DOUBLE PRECISION, DIMENSION(6,2) :: dtpdt
  DOUBLE PRECISION, DIMENSION(2,2) :: b
  DOUBLE PRECISION, DIMENSION(4,2) :: ceicel
  DOUBLE PRECISION, DIMENSION(6,6) :: v 
!================corrections for next newton step========================
  DOUBLE PRECISION, DIMENSION(2)   :: del2
  DOUBLE PRECISION, DIMENSION(4)   :: del4
  DOUBLE PRECISION, DIMENSION(6)   :: dels,delem,elem 
  INTEGER :: jj, iunnew0,iunwarn0
  INTEGER, DIMENSION(6) :: icor
  INTEGER :: inew,inter,ncor,itsav,iun20m,jcselold 
  LOGICAL :: succ,error,bizarre 
! =======================================================================
  TYPE(orbit_elem) :: el,el1          ! corrected elements 
  TYPE(orb_uncert) :: uncnew,unc1
  DOUBLE PRECISION csinor,csinew,delnew,ecc ! for diff_cor, bizarre    
  DOUBLE PRECISION :: delnor                ! corr. norm
  DOUBLE PRECISION :: v_infi,v_infty        ! vel at infinity      
  DOUBLE PRECISION tafter,tbefore           ! time of the final elements,
                                            ! time too early for close app. 
                                            ! in the same shower 
!===================new TP data, rescaling in AU==============================
  DOUBLE PRECISION               :: xnew,rnew,rescov 
  DOUBLE PRECISION, DIMENSION(2) :: dtpf,dtpcau
  INTEGER, PARAMETER :: nfoldx=10           ! max no of folds to give up
  INTEGER  :: it, nfold                     ! iteration control,fold counter
  INTEGER, PARAMETER :: itnewmax=12 
  INTEGER  :: i,j                           ! loops
  LOGICAL  :: batch                         ! batch control
! ============================================================================
  va_tracenew=va_tracemin
  fold=.false. 
  stored_vi=.false.
  iunnew0=abs(iunnew)
  iunwarn0=abs(iunwarn)
! nominal rms
  csinor=csinom(imi0)
! if already a VI, nothing to do                                        
  IF(va_tracenew%b.lt.b_e*1.01d0)THEN 
! return with collision                                                 
     distmin=va_tracenew%b 
! get the initial conditions used in the last run                       
     CALL lovinterp(va_tracenew%rindex,deltasig,el,uncnew,falsok)
! store information to document the VI found 
     curr_vi%ele=el
     curr_vi%tp=va_tracenew
     curr_vi%unc=uncnew
     stored_vi=.true.
     RETURN 
  ENDIF
! count folds
  nfold=0
! get the initial conditions used in the last run                       
  CALL lovinterp(va_tracenew%rindex,deltasig,el0,unc0,falsok) 
  uncnew=unc0
  DO 1 it=1,itnewmax 
! check for cases in which we are already inside the Earth cross section
     r=va_tracenew%b 
! current radius, displacement required to get to b_e                   
     dtpc(1)=-va_tracenew%opik%coord(4)*(r-b_e)/r 
     dtpc(2)=-va_tracenew%opik%coord(5)*(r-b_e)/r 
     sina=sin(va_tracenew%alpha) 
     cosa=cos(va_tracenew%alpha) 
     str=va_tracenew%stretch 
     wid=va_tracenew%width
     dq=((-dtpc(1)*sina+dtpc(2)*cosa)/str)**2                       &
     &        +((dtpc(1)*cosa+dtpc(2)*sina)/wid)**2
     IF(iunwarn.lt.0)  WRITE(*,*)' dq ', dq                   
     WRITE(iunwarn0,*)' dq ', dq 
! compute ellipse in the elements space 
     dtpdt=TRANSPOSE(dtpde(1:2,1:6,jcsel))                                
! WARNING!!!! should use uncnew%g, uncnew%c ????
!     CALL slinel(dtpdt,unc0%g,unc0%c,ceicel,b,v)
! the one below should be better
     CALL slinel(dtpdt,uncnew%g,uncnew%c,ceicel,b,v)
! WARNING ??????????????????????????????????????
! dtpc needs to be rescaled back into AU!!! 
     dtpcau(1)=dtpc(1)*reau 
     dtpcau(2)=dtpc(2)*reau 
! find corresponding orbital elements at epoch                          
     del2=MATMUL(b,dtpcau) 
     del4=MATMUL(ceicel,del2) 
     dels(1:2)=del2 
     DO jj=1,4 
        dels(2+jj)=-del4(jj) 
     ENDDO
     delem=MATMUL(v,dels) 
     el=el0 
     el%coord=el0%coord+delem 
     IF(iunwarn.lt.0)THEN
        WRITE(*,200)delem 
200     FORMAT('change in els ',1P,6D12.4) 
        WRITE(*,201)el%coo,el%coord 
201     FORMAT(' new elements ',A3,1x,1P,6D12.4)
     ENDIF
     WRITE(iunwarn0,200)delem
     WRITE(iunwarn0,201)el%coo,el%coord
     IF(bizarre(el,ecc))THEN 
        IF(iunwarn.lt.0)WRITE(*,*)' fold found, bizarre orbit'
        IF(iunnew.lt.0)WRITE(*,*)' fold found, bizarre orbit' 
        WRITE(iunwarn0,*)' fold found, bizarre orbit' 
        WRITE(iunnew0,*)' fold found, bizarre orbit'
        fold=.true. 
        distmin=r 
        RETURN 
     ENDIF
! solve for all elements                                                
     inter=0 
     CALL whicor(inter,icor,ncor,inew) 
! zero iterations differential corrections, batch mode                  
     itsav=itmax 
     itmax=0 
     iun20m=-1 
     batch=.true.
! differential corrections       
     CALL diff_cor(m_m,obs_m,obsw_m,el,icor,iun20m,el,uncnew,csinew,delnew,succ)
     itmax=itsav 
! time up to which to search for close approaches                       
     v_infi=v_infty(el) 
     CALL aftclov(iplam,el0%t,va_tracemin%tcla,v_infi,tbefore,tafter) 
! make new covariance matrix available for target plane analysis        
     CALL cov_avai(uncnew,el0%coo,el0%coord) 
! reset storage of close encounters                                     
     njc=0 
     iplam=0 
! propagation to search for new target plane point                      
     CALL pro_ele(el,tafter,el1,uncnew,unc1)
! check that TP has been achieved, at the right time                    
     csi=va_tracenew%opik%coord(4) 
     zeta=va_tracenew%opik%coord(5)
! check presence of target plane                                        
     IF(njc.eq.0.or.iplam.ne.3)THEN 
        IF(iunnew.lt.0)WRITE(*,*)' fold found, no TP'
        WRITE(iunwarn0,*)' fold found, no TP' 
        WRITE(iunnew0,*)' fold found, no TP' 
        fold=.true. 
        distmin=r 
        RETURN 
     ELSEIF(tcla(njc).lt.tbefore.or.tcla(1).gt.tafter)THEN 
        IF(iunwarn.lt.0)WRITE(*,*)' fold found, no TP'
        WRITE(iunwarn0,*)' fold found, no TP'
        WRITE(iunnew0,*)' fold found, no TP'
        fold=.true. 
        distmin=r 
        RETURN 
     ENDIF
! get data on TP analysys of the new close approach                     
     CALL arrloadtp(va_tracenew,va_tracemin%rindex) 
     IF(.not.va_tracenew%tp_conv)THEN
        WRITE(iunwarn0,*)' elliptic encounter, no TP '
        WRITE(iunnew0,*)'  elliptic encounter, no TP '
        fold=.true. 
        distmin=r 
        RETURN 
     ENDIF
! store value of sigma corresponding to RMS of residuals                
     IF(csinew.gt.csinor)THEN 
        va_tracenew%sigma=sqrt((csinew**2-csinor**2)*m_m*2)                  &
     &           /rescov(6,2*m_m,csinor)                                  
     ELSE 
        va_tracenew%sigma=0.d0 
     ENDIF
     WRITE(iunwarn0,*)'csinor,csinew,sigq ',csinor,csinew,va_tracenew%sigma
! control of termination (either convergence or rebound)                
     rnew=va_tracenew%b 
     dtpf(1)=va_tracenew%opik%coord(4)-csi 
     dtpf(2)=va_tracenew%opik%coord(5)-zeta 
     WRITE(iunwarn0,100)it,r,rnew 
100  FORMAT(' iteration ',i3,' r, rnew ',f8.3,1x,f8.3) 
     WRITE(iunwarn0,*)'dtp_target ',dtpc 
     WRITE(iunwarn0,*)'dtp_found  ',dtpf 
     IF(dtpf(1)*dtpc(1)+dtpf(2)*dtpc(2).lt.0.d0)THEN 
! fold passed, no way to reach to b_e                                   
        nfold=nfold+1
        WRITE(iunwarn0,*)' fold number ',nfold,'  angle > 90' 
        WRITE(iunnew0,*)' fold number ',nfold,'  angle > 90'
        IF(nfold.gt.nfoldx)THEN
           fold=.true. 
           distmin=r 
           RETURN 
        ELSE
           distmin=rnew
        ENDIF
        kill_propag=.false.
     ELSEIF(va_tracenew%sigma.gt.10.d0)THEN 
        nfold=nfold+1
        WRITE(iunwarn0,*)' fold number ',nfold,' sigma=', va_tracenew%sigma
        WRITE(iunnew0,*)' fold number ',nfold,' sigma=', va_tracenew%sigma 
        IF(nfold.gt.nfoldx)THEN
           fold=.true. 
           distmin=r 
           RETURN
        ELSE
           distmin=rnew
        ENDIF
        kill_propag=.false.
     ELSEIF(rnew.lt.b_e*1.01d0)THEN 
! convergence to a VI
        IF(iunwarn.lt.0) WRITE(*,*)' VI found '
        WRITE(iunwarn0,*)' VI found ' 
        WRITE(iunnew0,*)' VI found ' 
        fold=.false. 
        distmin=rnew
        kill_propag=.false. 
! store information to document the VI found 
        curr_vi%ele=el
        curr_vi%tp=va_tracenew
        curr_vi%unc=uncnew
        stored_vi=.true.
        RETURN 
     ELSE 
! adopt new starting point                                              
        distmin=rnew 
        el0=el
     ENDIF
1 ENDDO
! too many iteration                                                    
  WRITE(iunnew0,*)' b_e not reached after ',itnewmax,' iterations ' 
  WRITE(iunwarn0,*)' b_e not reached after ',itnewmax,' iterations ' 
! we do not know if it is really a fold, but it must not be in the risk 
  fold=.true.
  kill_propag=.false. 
END SUBROUTINE newton_checktp

! ==========================================================================
! RISKCHECKTP compute risk and write risk file if necessary               
! ==========================================================================
SUBROUTINE riskchecktp(va_tracemin,t0,type,no_risk,   &
     &     iunnew,iunwarn,iunrisk,riskfile)
  USE multiple_sol, ONLY: lovmagn
  USE planet_masses 
  USE tp_trace 
  USE virtual_impactor                                          
  USE eval_risk
! ========================INPUT============================================
  TYPE(tp_point), INTENT(IN)      :: va_tracemin
  CHARACTER(LEN=7), INTENT(IN)    :: type 
  DOUBLE PRECISION, INTENT(IN)    :: t0
  INTEGER, INTENT(IN)             :: iunnew,iunwarn 
  CHARACTER(LEN=*), INTENT(IN),OPTIONAL    :: riskfile
  INTEGER, INTENT(IN), OPTIONAL   :: iunrisk
! =======================OUTPUT============================================
  INTEGER, INTENT(INOUT)            :: no_risk 
! ====================END INTERFACE========================================
  DOUBLE PRECISION  :: b_e,dcur,vsize,width,stretch,alpha,csi1,zeta1,U 
  DOUBLE PRECISION  :: p_imp,h,mass,v_imp,energy,e_tilde,fb 
  DOUBLE PRECISION  :: rel_prob,ps,prob 
  DOUBLE PRECISION  :: tcl,sigma,deltat,sigimp
  DOUBLE PRECISION  :: chi,fact,p_imp1 
  INTEGER           :: le, iunrisk0,iunrisk1,iunnew0,iunwarn0
  CHARACTER(LEN=14) :: calend             ! calendar date 
  TYPE(tp_point)    :: va_tracenew        ! for newton_check     
  DOUBLE PRECISION  :: distmin,vvvv,rrrr,mu
  LOGICAL           :: fold 
! ==========================================================================
  dcur=va_tracemin%b 
  width=va_tracemin%width 
! on TP, gravitational focusing                                         
  U=va_tracemin%opik%coord(1)
  mu=gmearth/reau**3
  b_e=sqrt(1.d0+(2.d0*mu)/U**2)
  iunwarn0=abs(iunwarn)
  iunnew0=abs(iunnew) 
! unidimensional case is handled as narrow strip                        
  IF(width.lt.1.d-4)width=1.d-4 
  IF(dcur-10*width.lt.b_e)THEN 
! interesting cases with possibility of low passage
     IF(iunwarn.lt.0)THEN                    
        WRITE(*,*)' risk case analysed' 
        WRITE(*,*)'b_e, U ',b_e,U
     ENDIF
     WRITE(iunwarn0,*)' risk case analysed' 
     WRITE(iunwarn0,*)'b_e, U ',b_e,U 
     alpha=va_tracemin%alpha 
     csi1=va_tracemin%opik%coord(4)*cos(alpha)+va_tracemin%opik%coord(5)*sin(alpha) 
     zeta1=-va_tracemin%opik%coord(4)*sin(alpha)+va_tracemin%opik%coord(5)*cos(alpha) 
     stretch=va_tracemin%stretch 
     sigma=va_tracemin%sigma
! probability of impact                                                 
     p_imp=prob(csi1,zeta1,sigma,width,stretch,b_e) 
     IF(iunwarn.lt.0)WRITE(*,*)'Impact Probability ',p_imp
     WRITE(iunwarn0,*)'Impact Probability ',p_imp
! write risk file if appropriate                                        
     IF(p_imp.gt.1e-11)THEN 


! palermo scale: need time span and magnitude, U                        
        tcl=va_tracemin%tcla 
        deltat=tcl-t0 
        rrrr=va_tracemin%rindex
        CALL lovmagn(rrrr,vvvv,h) ! vvvv is v_infty, replaced by current U 
        ps=palermo(U*reau,h,p_imp,deltat,                            &
     &           mass,v_imp,energy,e_tilde,fb,rel_prob,iunwarn)         
! calendar date                                                         
        CALL calendwri(tcl,calend)
        sigimp=MAX((abs(csi1)-b_e)/width,0.d0)
! before writing, check for spurious VI due to interrupted              
        WRITE(iunwarn0,*)' checking for risk of type ',type 
        CALL newton_checktp(iunnew,iunwarn,va_tracemin,b_e,distmin,va_tracenew,fold)
        IF(fold)THEN 
           RETURN 
        ELSEIF(va_tracemin%b.eq.distmin)THEN 
! nothing to change                                                     
           p_imp1=p_imp 
        ELSE 
! correct probability for value of chi**2                               
           chi=va_tracenew%sigma 
           fact=exp(-0.5d0*(chi**2-sigma**2-sigimp**2)) 
           p_imp1=p_imp*fact 
           ps=palermo(vvvv,h,p_imp1,deltat,                        &
     &              mass,v_imp,energy,e_tilde,fb,rel_prob,iunwarn)      
           IF(iunwarn.lt.0)WRITE(*,166)calend,p_imp1,fact
           WRITE(iunwarn0,166)calend,p_imp1,fact
166        FORMAT(a14,' probability recomputed is ',1p,d11.3,' fact= ', d11.3) 
           IF(p_imp1.lt.1.d-11)THEN 
              IF(iunnew.lt.0)WRITE(*,*)' spurious low prob=', p_imp1
              WRITE(iunnew0,*)' spurious low prob=', p_imp1 
              RETURN 
           ENDIF
           IF(p_imp1.gt.1.d0)THEN 
              IF(iunnew.lt.0)WRITE(*,*)' spurious high prob=', p_imp1 
              WRITE(iunnew0,*)' spurious high prob=', p_imp1 
              RETURN 
!           ELSEIF(ps.gt.0.d0.and.fact.gt.1.5.d0)THEN 
!              WRITE(iunnew0,*)' spurious high PS=', ps 
!              RETURN 
           ENDIF
        ENDIF
! store information on size of Earth impact cross section               
        IF(b_e.gt.bsdmax)bsdmax=b_e 
        IF(b_e.lt.bsdmin)bsdmin=b_e 
        no_risk=no_risk+1 
        IF(PRESENT(riskfile))THEN
           CALL rmsp(riskfile,le) 
           iunrisk0=44 
           iunrisk1=iunrisk0
           OPEN(UNIT=iunrisk0, FILE=riskfile(1:le),POSITION='APPEND')
        ELSEIF(PRESENT(iunrisk))THEN
           iunrisk0=abs(iunrisk)
           iunrisk1=iunrisk
        ELSE
           WRITE(*,*)'riskchecktp: missing output optional arguments'
           STOP
        ENDIF
        IF(iunrisk1.le.0) CALL header_risk(0)
        IF(width.ge.1.0d4)THEN 
           IF(iunrisk1.lt.0)WRITE(*,200)calend,tcl,sigma,sigimp,dcur,          &
     &              width,stretch,p_imp1,e_tilde,ps
           WRITE(iunrisk0,200)calend,tcl,sigma,sigimp,dcur,          &
     &              width,stretch,p_imp1,e_tilde,ps                     
200        FORMAT(a14,1x,f9.3,1x,f6.3,1x,f5.3,1x,f7.2,' +/- ',      &
     &              f8.2,1x,1p,e9.2,1x,e9.2,1x,e9.2,1x,0p,f6.2)         
        ELSE
           IF(iunrisk1.lt.0)WRITE(*,100)calend,tcl,sigma,sigimp,dcur,          &
     &              width,stretch,p_imp1,e_tilde,ps
           WRITE(iunrisk0,100)calend,tcl,sigma,sigimp,dcur,          &
     &              width,stretch,p_imp1,e_tilde,ps                     
100        FORMAT(a14,1x,f9.3,1x,f6.3,1x,f5.3,1x,f7.2,' +/- ',      &
     &              f8.3,1x,1p,e9.2,1x,e9.2,1x,e9.2,1x,0p,f6.2)         
        ENDIF
        IF(PRESENT(riskfile))CLOSE(iunrisk0) 
     ENDIF
  ENDIF

END SUBROUTINE riskchecktp

SUBROUTINE header_risk(iunrisk)
INTEGER, INTENT(IN) :: iunrisk ! output unit
  WRITE(iunrisk,200) '     date        MJD      sigma sigimp  ',    &
     &     ' dist +/-   width    stretch    p_RE    exp. en.    PS  '   
  WRITE(iunrisk,200) '   YYYY/MM                              ',    &
     &     ' (RE)       (RE)      RE/sig               MT           '   
  WRITE(iunrisk,200) '----------------------------------------',    &
     &     '--------------------------------------------------------'   
200 FORMAT(a40,a56) 
END SUBROUTINE header_risk
  
  SUBROUTINE big_vitp(tptrail, lre, tpmin, t0,nvai, nbigrisk,iunnew,iunwarn,riskfile)
    USE planet_masses, ONLY: gmearth
    USE tp_trace
    USE multiple_sol, ONLY: lovmagn
    USE eval_risk
    INTEGER, INTENT(IN) :: lre
    TYPE(tp_point), INTENT(IN)  :: tptrail(lre),tpmin
    DOUBLE PRECISION, INTENT(IN) :: t0 ! intial conditions time
    INTEGER, INTENT(IN) :: iunnew, iunwarn ! output units
    CHARACTER(LEN=*)  :: riskfile
    INTEGER, INTENT(OUT) :: nvai
    INTEGER, INTENT(OUT) :: nbigrisk
    LOGICAL, DIMENSION(lre) :: viflag
    INTEGER i,j,nvai_int,interv(2,10), iunrisk,le,ii2
    DOUBLE PRECISION v, b_e, be, mu, stretch, width, dvai, p_imp, gau, fb, rel_prob,sigma, dcur, rtp
    DOUBLE PRECISION rindex, vinf, sigimp, tcl, h, mass, deltat, ps, v_imp, energy, e_tilde,dsigma
    CHARACTER(LEN=14) :: calend             ! calendar date 
! count Virtual Asteroids Impacting (VAI)
    nvai=0
    nvai_int=0
    viflag=.false.
    mu=gmearth/reau**3
! on TP, gravitational focusing 
    v=tpmin%opik%coord(1) ! this is U if conversion to TP has taken place
    be=sqrt(1.d0+(2.d0*mu)/v**2) ! this formula is true if tp_coord are on TP
    DO i=1,lre
!       tp_flag_cur=tptrail(i)%tp_conv
!       IF(tp_flag_cur)THEN 
          b_e=be
       IF(tptrail(i)%b.le.b_e)THEN
          viflag(i)=.true.
          nvai=nvai+1
          IF(nvai.eq.0)THEN
             nvai_int=1
             interv(1,1)=i
             interv(1,2)=i
          ELSE
             IF(i.gt.1)THEN
                IF(.not.viflag(i-1))THEN
                   nvai_int=nvai_int+1
                   interv(nvai_int,1)=i
                   interv(nvai_int,2)=i
                ELSE
                   interv(nvai_int,2)=i
                ENDIF
             ELSE
                nvai_int=nvai_int+1
                interv(nvai_int,1)=i
                interv(nvai_int,2)=i
             ENDIF
         ENDIF
!          WRITE(*,*) tptrail(i)%b, b_e          
       ENDIF
    ENDDO
    WRITE(iunwarn,*)' risk case analysed' 
    WRITE(iunwarn,*)'b_e, v ',b_e,v 
    IF(nvai_int.gt.10)STOP ' *********big_vitp: nvai_int .gt.10 *************'
    IF(nvai.eq.0) RETURN
! there are risk cases
    WRITE(iunwarn,*)' nvai ', nvai
    WRITE(iunwarn,*)'intervals of vai ', nvai_int
    DO j=1,nvai_int
!      WRITE(iunwarn,*)' from ', interv(j,1), ' to ',interv(j,2)
      WRITE(iunwarn,*)' from ', tptrail(interv(j,1))%rindex, ' to ',tptrail(interv(j,2))%rindex
      IF(interv(j,2)-interv(j,1)+1.lt.10)THEN
         WRITE(iunwarn,*)' too few VAI ',interv(j,2)-interv(j,1)+1
         CYCLE
      ENDIF 
      width=0.d0
      p_imp=0.d0
      dcur=100.d0

      DO i=interv(j,1),interv(j,2)
        IF(i.eq.interv(j,1)) dsigma=abs(tptrail(i)%sigma-tptrail(i+1)%sigma)
        width=width+tptrail(i)%width
        gau=exp(-tptrail(i)%sigma**2/2)*(1.d0/sqrt(dpig))
        p_imp=p_imp+gau*dsigma ! for a tail it would crash
        rtp=sqrt(tptrail(i)%opik%coord(4)**2+tptrail(i)%opik%coord(5)**2)
        dcur=MIN(dcur,rtp)
      ENDDO
      WRITE(iunwarn,*)'Impact Probability ',p_imp 
      width=width/(interv(j,2)-interv(j,1)+1) 
      dvai=sqrt((tptrail(interv(j,2))%opik%coord(4)-tptrail(interv(j,1))%opik%coord(4))**2 + &
&      (tptrail(interv(j,2))%opik%coord(5)-tptrail(interv(j,1))%opik%coord(5))**2)
      stretch=dvai/abs(tptrail(interv(j,2))%sigma-tptrail(interv(j,1))%sigma)
!      WRITE(*,*)' dvai, stretch,wid, p_imp ',dvai,stretch, wid, p_imp
! write risk file if appropriate                                        
      IF(p_imp.gt.1e-11)THEN 
         nbigrisk=nbigrisk+1
! palermo scale: need time span and magnitude, U                        
         tcl=tpmin%tcla 
         deltat=tcl-t0 
         rindex=interv(j,1)
         CALL lovmagn(rindex,vinf,h) 
         ps=palermo(vinf,h,p_imp,deltat,                            &
     &           mass,v_imp,energy,e_tilde,fb,rel_prob,iunwarn)         
! calendar date                                                         
         CALL calendwri(tcl,calend)
         sigimp=0.d0
         WRITE(iunwarn,*)' Virtal Asteroids Impacting found, no check needed'
! store information on size of Earth impact cross section               
         IF(b_e.gt.bsdmax)bsdmax=b_e 
         IF(b_e.lt.bsdmin)bsdmin=b_e 
         CALL rmsp(riskfile,le) 
         iunrisk=44 
         OPEN(UNIT=iunrisk, FILE=riskfile(1:le),POSITION='APPEND')   
         ii2=(interv(j,1)+interv(j,2))/2
         sigma=tptrail(ii2)%sigma                                
         IF(width.ge.1.0d4)THEN 
            WRITE(iunrisk,200)calend,tcl,sigma,sigimp,dcur,          &
     &              width,stretch,p_imp,e_tilde,ps                     
200         FORMAT(a14,1x,f9.3,1x,f6.3,1x,f5.3,1x,f7.2,' +/- ',      &
     &              f8.2,1x,1p,e9.2,1x,e9.2,1x,e9.2,1x,0p,f6.2)         
         ELSE 
            WRITE(iunrisk,100)calend,tcl,sigma,sigimp,dcur,          &
     &              width,stretch,p_imp,e_tilde,ps                     
100         FORMAT(a14,1x,f9.3,1x,f6.3,1x,f5.3,1x,f7.2,' +/- ',      &
     &              f8.3,1x,1p,e9.2,1x,e9.2,1x,e9.2,1x,0p,f6.2)         
         ENDIF
      ENDIF
      CLOSE(iunrisk)
    ENDDO
  END SUBROUTINE big_vitp

!=======================================================================
END MODULE offlov_checktp







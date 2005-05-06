! ========MODULE target_plane==================  
! Last changes: 17/02/2003
! CONTAINS                                                              
! PUBLIC                                                                
!             fclan2   close approach analysis driver  
! ROUTINES                                                              
!             mtpsel2                                                   
!             newton_tp                                                 
!             tpslin                                                    
!             virimp                                                    
!  
! OUT of the module:
!             aftclov 
!             v_infty   vel at infinity (w.r. to Earth only) 
!                                                                       
! HEADERS       
!       target_plane.o: \
!	../include/nvarx.h90 \  for varwra in strclan3 
!	../suit/astrometric_observations.mod \
!	close_app.o \  module close_app_data
!	least_squares.o \  module
!	tp_trace.o         module  
!                                                                       
! ================================================

MODULE target_plane

USE output_control
USE orbit_elements
USE util_suit
IMPLICIT NONE

PRIVATE

! PUBLIC SUBROUTINES
PUBLIC fclan2 
CONTAINS
! =======================================                               
!  FCLAN2                                                               
! vers. 2.2.5, A. Milani, 5 april 2001                                  
! =======================================                               
! close approach analysis                                               
! ===============INTERFACE========================                      
SUBROUTINE fclan2(batchcl,t1,iunout,ok,siglim,              &
     &   el0,unc0,csinor,delnor,astna0,m,obs,obsw)
  USE astrometric_observations
  USE tp_trace
  USE propag_state 
! ================INPUT===========================                      
! control for batch mode                                                
  LOGICAL, INTENT(IN) :: batchcl 
! time to monitor up to, passed only in batch mode                      
  DOUBLE PRECISION, INTENT(IN) :: t1 
! max value of sigma change allowed for newton step                     
  DOUBLE PRECISION,INTENT(INOUT) :: siglim 
! ======observations==== 
  INTEGER,INTENT(IN) ::  m ! number of observations
! new data types
  TYPE(ast_obs),DIMENSION(m), INTENT(IN) :: obs
  TYPE(ast_wbsr),DIMENSION(m), INTENT(IN) :: obsw
 ! =======orbit==============  
  TYPE(orbit_elem), INTENT(INOUT) :: el0
  TYPE(orb_uncert), INTENT(INOUT) :: unc0
! asteroid names (18 CHARACTERs)                                        
  CHARACTER*18, INTENT(IN) :: astna0 
! norms of residuals, of last correction                                
  DOUBLE PRECISION, INTENT(INOUT) ::  csinor,delnor 
! output unit, verbosity                                                
  INTEGER,INTENT(IN) :: iunout 
! ===============HIDDEN INPUT========================                   
! arrays of close approach data and controls                            
! ===============OUTPUT=========================                        
! success flag                                                          
  LOGICAL ok 
! ==============HIDDEN OUTPUT===================                        
! covariance matrix passed through hidden common                        
! ==============END INTERFACE===================== 
  TYPE(orbit_elem) :: el1
  TYPE(orb_uncert) :: unc1
! target plane type, nonlinearity handling, flag for ephemerides, verbosity
  INTEGER iclan,inl,ikey,iverb 
! logical controls                                                      
  LOGICAL error 
! menu name                                                      
  CHARACTER*20 menunam 
! final elements after propagation                                      
! DOUBLE PRECISION tcl
! time of close approach passed to subroutines                          
  DOUBLE PRECISION tclo 
! initialisation                                                        
  INTEGER lflag 
! memory model                                                          
  SAVE 
  DATA lflag /0/ 
  IF(lflag.eq.0)THEN 
     ipla0=0 
     lflag=1 
  ENDIF
! ================================================================      
! t1 has been passed as dummy variable, nothing to do                   
! initialize close approach count                                       
  njc=0 
  iplam=0 
! check availability of JPL ephemerides                                 
  CALL chetim(el0%t,t1,ok) 
  IF(.not.ok)THEN 
     WRITE(*,*)' JPL ephemerides not available for t1=',t1 
     ok=.false. 
     RETURN 
  ENDIF
! availability of covariance                                            
  CALL cov_avai(unc0,el0%coo,el0%coord) 
! explorative propagation, incorporates linear target plane analysis    
  CALL pro_ele(el0,t1,el1,unc0,unc1) 
! select target plane                                                   
  IF(batchcl)THEN 
     iverb=1 
  ELSE 
     iverb=20 
  ENDIF
  CALL mtpsel2(batchcl,iverb,error) 
  IF(error)THEN 
     ok=.false. 
     RETURN 
  ENDIF
! ===================================================================== 
! selection of the procedure to use the target plane information        
! ===================================================================   
! restart from menu                                                     
58 CONTINUE 
! choose handling of nonlinearity                                       
  IF(batchcl)THEN 
! batch use is for Newton's method (from clonew.f)                      
     inl=4 
  ELSE 
     menunam='closnonl' 
     CALL menu(inl,menunam,4,'How to handle nonlinearity?=',        &
     &        'linear map=',                                            &
     &        'semilinear n-body=',                                     &
     &        'Newtons method on LOV, 1 step=',                         &
     &        'Newtons method on LOV, auto=')
! exit                                                                  
     IF(inl.eq.0)RETURN 
  ENDIF
! ===================================================================   
! operation mode depending upon inl                                     
  IF(inl.eq.1.or.inl.eq.2)THEN 
! ============target plane semilinear (and linear) analysis:            
     CALL cov_not_av ! semilinear method propagates without covariance
     CALL tpslin(batchcl,inl,iunout,astna0,el0,unc0) 
     CALL cov_avai(unc0,el0%coo,el0%coord) ! restore covariance availability
  ELSEIF(inl.eq.3.or.inl.eq.4)THEN 
! ============Newton's method to find virtual impactors:                
     IF(.not.batchcl)THEN 
        WRITE(*,*)' assign limit to change in sigma, 0=default' 
        READ(*,*)siglim 
        IF(siglim.eq.0.d0)siglim=3.d0 
     ENDIF
     tclo=tcla(jcsel) 
     CALL newton_tp(batchcl,inl,iunout,siglim,ok,                    &
  &        tclo,astna0,el0,unc0,csinor,delnor,m,obs,obsw)
     IF(batchcl)THEN 
        RETURN 
     ENDIF
  ENDIF
! ===========================================================           
! back to the menu                                                      
  GOTO 58 
END SUBROUTINE fclan2
! ===========================================================
!  MTPSEL2                                                              
! selects target plane (in case there are many)                         
! also handles output in interactive case                               
! in batch, it selects the lowest minimum distance                      
! ============================================================
SUBROUTINE mtpsel2(batchcl,verbose,error) 
  USE tp_trace
  USE planet_masses
! interface                                                             
! input: batch? write?                                                  
  LOGICAL, INTENT(IN) :: batchcl 
  INTEGER, INTENT(IN) :: verbose 
! output: something wrong                                               
  LOGICAL, INTENT(OUT) :: error 
! end interface                                                         
! hidden input:                                                         
! arrays of close approach data and controls                            
! basis adapted to MTP output from mtprot                               
! planet names                                                          
  CHARACTER*30 planam 
  INTEGER lpla,lench 
! index of target planes, indexes of coordinates on tp                  
  INTEGER jc,i,j,jj,k 
! minimum distance, angle on TP                                         
  DOUBLE PRECISION dist_min,tpth 
! are there close approaches?                                           
  IF(iplam.eq.0.or.njc.eq.0)THEN 
     WRITE(*,*)' no close approaches found', iplam, njc 
! error flag                                                            
     error=.true. 
     RETURN 
  ELSE 
     error=.false. 
     planam=ordnam(iplam) 
     lpla=lench(planam) 
     IF(.not.batchcl)THEN 
        WRITE(*,101)planam(1:lpla) 
101     FORMAT(' close approach to planet ',a) 
     ENDIF
  ENDIF
! selection of minimum                                                  
  IF(.not.batchcl.or. verbose.gt.19)THEN 
     WRITE(*,102)(jc,tcla(jc),rmin(jc),jc=1,njc) 
102  FORMAT(' time(MJD), min.dis.(AU), normal to MTP'/              &
     &        (i3,1x,f12.5,1x,f11.8/))                                  
  ENDIF
  IF(njc.eq.1)THEN 
     IF(batchcl)THEN 
        jcsel=1 
     ELSE 
198     WRITE(*,*)' select target plane 1=yes, 0=quit' 
        READ(*,*)jcsel 
        IF(jcsel.eq.0)THEN 
           WRITE(*,*)' no target plane analysis' 
           error=.true. 
           RETURN 
        ELSEIF(jcsel.lt.0.or.jcsel.gt.njc)THEN 
           WRITE(*,*) 'should be 0< ', jcsel,' < ',njc+1 
           GOTO 198 
        ENDIF
     ENDIF
  ELSE 
     IF(batchcl)THEN 
! selects the minimum; is this right???                                 
        dist_min=1.d2 
        DO jc=1,njc 
           IF(tpr(jc).lt.dist_min)THEN 
              dist_min=tpr(jc) 
              jcsel=jc 
           ENDIF
        ENDDO
     ELSE 
199     WRITE(*,*)' select target plane among the above,0=quit' 
        READ(*,*)jcsel 
        IF(jcsel.eq.0)THEN 
           WRITE(*,*)' no target plane analysis' 
           error=.true. 
           RETURN 
        ELSEIF(jcsel.lt.0.or.jcsel.gt.njc)THEN 
           WRITE(*,*) 'should be 0< ', jcsel,' < ',njc+1 
           GOTO 199 
        ENDIF
     ENDIF
  ENDIF
! output (only in interactive mode)                                     
  IF(verbose.lt.9) RETURN 
  tpth=atan2(tpc(2,jcsel),tpc(1,jcsel)) 
  WRITE(*,170) (tpc(j,jcsel),j=1,3),tpr(jcsel),tpth 
170 FORMAT('MTP coordinates ',3f12.8,' dist ',f12.8,' angle ',f8.4)
  WRITE(*,*)' partial derivatives' 
  DO k=1,3 
     WRITE(*,171) (dtpdet(jj,k,jcsel),jj=1,6) 
  ENDDO
171 FORMAT(6(f11.3,1x)) 
  WRITE(*,173)(sig(i,jcsel),i=1,2) 
173 FORMAT(' semiaxes of target ellipse: ',1p,2d10.3) 
  WRITE(*,172)((axes(i,j,jcsel),i=1,2),j=1,2) 
172 FORMAT(' directions ',2f10.6/12x,2f10.6) 
END SUBROUTINE mtpsel2
! ===========================================================================
! NEWTON_TP                                                               
! called by fclan2                                                      
! target plane Newton's method to find minimum distance                 
! on long axis of TP ellipse                                            
! ===========================================================================
SUBROUTINE newton_tp(batchcl,inl,iunout,siglim,ok,                 &
     &     tclo,astna0,el0,unc0,csinor,delnor,            &
     &     m,obs,obsw)
  USE astrometric_observations
  USE least_squares
  USE tp_trace
  USE propag_state 
! INPUT                                                                 
  LOGICAL, INTENT(IN) :: batchcl ! batch control
  INTEGER, INTENT(IN) :: inl ! inl=1 linear inl=2 semilinear 
  INTEGER, INTENT(IN) :: iunout ! log unit 
! observations
  INTEGER, INTENT(IN) ::  m ! number of observations
! new data types
  TYPE(ast_obs),DIMENSION(m), INTENT(IN) :: obs
  TYPE(ast_wbsr),DIMENSION(m), INTENT(IN) :: obsw
! max value of sigma change allowed for newton step                     
  DOUBLE PRECISION, INTENT(IN) :: siglim 
! time of close approach being explored                                 
  DOUBLE PRECISION, INTENT(IN) :: tclo 
! =========elements=================                                    
! asteroid names (18 CHARACTERs)                                        
  CHARACTER*18, INTENT(IN) :: astna0 
! =======orbit==============  
  TYPE(orbit_elem), INTENT(INOUT) :: el0
  TYPE(orb_uncert), INTENT(INOUT) :: unc0
! residuals norm (for nominal orbit), convergence of last diff. corr.   
  DOUBLE PRECISION, INTENT(INOUT) :: csinor,delnor 
! =========observations=================                                
! HIDDEN INPUT                                                          
! max distance (from nominal point on target plane) at which we can     
! reasonably apply a linear approx                                      
  DOUBLE PRECISION, PARAMETER  :: linlim=0.02d0
! arrays of close approach data                                         
! basis adapted to MTP output from mtprot                               
! OUTPUT                                                                
! success?                                                              
  LOGICAL, INTENT(OUT) ::  ok 
! ==================END INTERFACE=====================                  
! ====== new starting point=================                            
! correspondent ellipse in orbital elements space                       
  DOUBLE PRECISION b(2,2),ceicel(4,2),v(6,6) 
! variables for Newton's method                                         
  DOUBLE PRECISION s1,tpc1(2),dtpc(2),s1max 
  LOGICAL s1_change 
  TYPE(orbit_elem) :: elc,el1
  TYPE(orb_uncert) :: uncc,unc1
  DOUBLE PRECISION del2(2),del4(4),dels(6),delem(6),elem(6) ,elemp(6)
! new expected min.distance, normal and for shortened step,change in dis
  DOUBLE PRECISION dminew,dminew2,deldis 
! ================for difcor=================                           
  DOUBLE PRECISION csinew,delnew,cnew(6,6),gnew(6,6) 
  INTEGER icor(6),inew,inter,ncor,itsav,iun20m,jcselold 
  LOGICAL succ,error 
! =============== for new propagation===============                    
! vel at infinity                                                       
  DOUBLE PRECISION v_inf,v_infty 
! time of the final elements, new elements, covariance                  
! time too early for close app. in the same shower                      
  DOUBLE PRECISION tafter,eq1(6),c1(6,6),g1(6,6),tbefore 
! newton iterations                                                     
  INTEGER itnewma,it,it_too_long, itlongmax,iverb 
  DOUBLE PRECISION newcont, disnow 
!  moid routine                                                         
  DOUBLE PRECISION moid,dnp,dnm 
  INTEGER iconv 
! ephem of VI                                                           
  INTEGER ikey 
  CHARACTER*18 astnavi
  CHARACTER*60 file
! =================================                                     
! loop indexes                                                          
  INTEGER jj 
! =============================================                         
! Newton's method setup                                                 
  IF(.not.batchcl)WRITE(*,*)' Newton method' 
  ok=.false. 
  IF(inl.eq.3)THEN 
     itnewma=1 
  ELSEIf(inl.eq.4)THEN 
     itnewma=10 
  ENDIF
  itlongmax=5 
  newcont=1.d-6 
  it_too_long=0 
! iteration loop                                                        
  DO 1 it=1,itnewma 
!         IF(.not.batchcl)                                              
     WRITE(*,*)' ===== iteration ',it,' ========' 
     disnow=sqrt(tpc(1,jcsel)**2+tpc(2,jcsel)**2) 
! find closest point on MTP according to linear approx                  
     s1=-tpc(1,jcsel)*axes(1,2,jcsel)-tpc(2,jcsel)*axes(2,2,jcsel) 
! new target point on MTP                                               
     DO jj=1,2 
        dtpc(jj)=s1*axes(jj,2,jcsel) 
        tpc1(jj)=tpc(jj,jcsel)+dtpc(jj) 
     ENDDO
     dminew=sqrt(tpc1(1)**2+tpc1(2)**2) 
!         IF(.not.batchcl)                                              
     WRITE(*,120) s1,s1/sig(2,jcsel),dminew 
120  FORMAT(' displ.=',1p,d9.2,' deltasig=',d9.2,' dmin=',d9.2) 
     WRITE(iunout,*) s1,s1/sig(2,jcsel),dminew 
! control on the size of displacement s1 (in AU)                        
     s1_change=.false. 
     s1max=linlim 
     IF(abs(s1).gt.s1max)THEN 
!            IF(.not.batchcl)                                           
        WRITE(*,121)s1,s1max 
121     FORMAT(' Newton step too long s1=',1p,d9.2,' max=',d9.2) 
        s1=s1*s1max/abs(s1) 
        s1_change=.true. 
     ENDIF
! control on the size of change in the sigma space                      
     IF(abs(s1/sig(2,jcsel)).gt.siglim)THEN 
        WRITE(*,122)s1/sig(2,jcsel),siglim 
122     FORMAT(' too large deltasig=',1p,d9.2,' max=',d10.2) 
!           ok=.false.                                                  
!           RETURN                                                      
! experiment!!!                                                         
        it_too_long=it_too_long+1 
        s1=s1*siglim/abs(s1/sig(2,jcsel)) 
        s1_change=.true. 
     ENDIF
! new target point on MTP                                               
     IF(s1_change)THEN 
        DO jj=1,2 
           dtpc(jj)=s1*axes(jj,2,jcsel) 
           tpc1(jj)=tpc(jj,jcsel)+dtpc(jj) 
        ENDDO
        dminew2=sqrt(tpc1(1)**2+tpc1(2)**2) 
!           IF(.not.batchcl)                                            
        WRITE(*,120) s1,s1/sig(2,jcsel),dminew2 
        WRITE(iunout,120) s1,s1/sig(2,jcsel),dminew 
     ENDIF
! compute ellipse in the elements space                                 
     CALL slinel(dtpdet(1,1,jcsel),unc0%g,unc0%c,ceicel,b,v) 
! find corresponding orbital elements at epoch                          
     del2=MATMUL(b,dtpc) 
     del4=MATMUL(ceicel,del2) 
     dels(1:2)=del2 
     dels(3:6)=-del4(1:4) 
     delem=MATMUL(v,dels)
     elc=el0 
     elc%coord=el0%coord+delem 
! solve for all elements                                                
     inter=0 
     CALL whicor(inter,icor,ncor,inew) 
! zero iterations differential corrections, batch mode                  
     itsav=itmax 
     itmax=0 
     iun20m=-iunout 
!        batch=.true.                                                   
     CALL diff_cor(m,obs,obsw,elc,icor,iun20m,elc,uncc,csinew,delnew,succ)
     itmax=itsav 
! time up to which to search for close approaches                       
     v_inf=v_infty(elc) 
     CALL aftclov(iplam,el0%t,tclo,v_inf,tbefore,tafter) 
! make new covariance matrix available for target plane analysis        
     CALL cov_avai(uncc,el0%coo,el0%coord) 
! reset storage of close encounters                                     
     njc=0 
     iplam=0 
! propagation to search for new target plane point                      
     CALL pro_ele(elc,tafter,el1,uncc,unc1) 
! select target plane: the same local minimum!!!????                    
     jcselold=jcsel 
!    IF(batchcl)THEN                                                 
!        iverb=1                                                      
!    ELSE                                                            
     iverb=20 
!    ENDIF                                                           
     CALL mtpsel2(.true.,iverb,error) 
! different causes of failure                                           
     IF(error)THEN 
        ok=.false. 
        WRITE(*,*)' No target plane, Newton failed iter ',it 
        RETURN 
     ELSEIF(jcselold.ne.jcsel)THEN 
        WRITE(*,*)' change of jcsel, was ',jcselold,' now ',jcsel 
! should we stop????                                                    
!       WRITE(*,*)' Newton failed at iteration ',it                  
!       ok=.false.                                                   
!       RETURN                                                       
     ELSEIF(tcla(jcsel).gt.max(tbefore,tafter).or.                   &
     &          tcla(jcsel).lt.min(tbefore,tafter))THEN                 
        WRITE(*,130)tcla(jcsel), tbefore, tafter 
130     FORMAT(' time_cla=',f11.2,' out of range, must be in '       &
     &          ,2f11.2)                                                
        WRITE(*,*)' Newton failed iter ',it 
        ok=.false. 
        RETURN 
     ENDIF
! if there is an acceptable close approach, adopt the newton orbit      
! as nominal                                                            
     el0=elc 
     unc0=uncc 
     csinor=csinew 
     delnor=delnew 
! assess where we stand after this iteration                            
     deldis=tpr(jcsel)-disnow 
     disnow=tpr(jcsel) 
     write(*,*)' Newton controls=',dminew-disnow,deldis 
     write(iunout,*)' Newton controls=',dminew-disnow,deldis 
     IF(abs(dminew-disnow).lt.newcont.or.abs(deldis).lt.newcont)THEN 
        ok=.true. 
        IF(batchcl)THEN 
! Output the data for resret.pl:                                        
!($time,$mindist,$tpc[0],$tpc[1],$semiwidth,$stretch,$rms,$nconv)=      
!         write (iunout,145) tcla(1),rmin(1),tpc(1),tpc(2),           
!     +           sig(1),sig(2), csinor, dminew-tpr                     
!         write (*,145) tcla(1),rmin(1),tpc(1),tpc(2),               
!     +            sig(1),sig(2),csinor,dminew-tpr                      
! assuming reference system on MTP has first axis along the             
! direction of the last close approach found                            
           write (iunout,145) tcla(jcsel),disnow,disnow,0.d0,         &
     &             sig(1,jcsel),sig(2,jcsel), csinor, dminew-disnow     
           write (*,145) tcla(jcsel),disnow,disnow,0.d0,             &
     &             sig(1,jcsel),sig(2,jcsel),csinor,dminew-disnow       
  145         format(f14.7,7(1x,g18.10)) 
        ELSE 
! convergence already achieved                                          
           write(*,*)' newton converged' 
           astnavi='virimp'
           CALL  write_elems(el0,astnavi,'ML',file,iunout,unc0)                     
           CALL nomoid(el0%t,el0,moid,dnp,dnm) 
           write(iunout,198)moid,dnp,dnm 
198        format('!MOID ',f8.5/'!NODES ',f8.5,1x,f8.5) 
           write(*,*)' ephemeris of virtual impactor? 1=yes 0=no' 
           read(*,*)ikey 
           IF(ikey.eq.1)THEN 
              CALL virimp(tpc(1,jcsel),dtpdet(1,1,jcsel),           &
   &          axes(1,1,jcsel),sig(1,jcsel),ceicel,b,v,         &
     &                 el0,unc0,iunout)                    
           ENDIF
        ENDIF
        GOTO 2 
     ELSEIF(it_too_long.ge.itlongmax)THEN 
        GOTO 3 
     ENDIF
! end iteration loop                                                    
1 ENDDO
  IF(.not.batchcl.and.inl.eq.4)THEN 
     WRITE(*,*)' Done ',itmax,' Newton iterations' 
  ENDIF
!     controversial choice: if Newton has not converged, but not gone   
!     out of target plane, then take where it stops as output           
  ok=.true. 
  RETURN 
3 IF(.not.batchcl.and.inl.eq.4)THEN 
     WRITE(*,*)' Done ',it,' Newton iterations ', it_too_long,      &
     &        ' too long'                                               
  ENDIF
!     controversial choice: if Newton has stopped because of too many   
!     short steps, but not gone out of target plane, then take where it 
!     stops as output                                                   
  ok=.true. 
  RETURN 
2 IF(.not.batchcl)WRITE(*,*)' Newton convergence at iter. ',it 
  RETURN 
END SUBROUTINE newton_tp
! ========================================================================
! TPSLIN                                                                
! called by fclan2                                                      
! target plane semilinear (and linear) analysis:                        
! ========================================================================
SUBROUTINE tpslin(batchcl,inl,iunout,astna0,el0,unc0) 
  USE tp_trace
  USE propag_state 
! INPUT                                                                 
! batch control                                                         
  LOGICAL, INTENT(IN) ::  batchcl 
! inl=1 linear inl=2 semilinear                                         
  INTEGER, INTENT(IN) :: inl 
! log unit                                                              
  INTEGER, INTENT(IN) :: iunout 
! asteroid names (18 CHARACTERs)                                        
  CHARACTER*18, INTENT(IN) :: astna0 
! nominal initial conditions: epoch, elements  
  TYPE(orbit_elem), INTENT(INOUT) :: el0
! normal and covariance matrices                                        
  TYPE(orb_uncert), INTENT(INOUT) :: unc0
! HIDDEN INPUT/OUTPUT                                                          
! arrays of close approach data                                         
! basis adapted to MTP output from mtprot                               
! END INTERFACE                                                         
! propagations: v_inf, close approach times                
  DOUBLE PRECISION v_inf,v_infty,tbefore,tafter,tcla0
  TYPE(orbit_elem) :: el1, elv
! ====== multiple target plane intersections============                
! string of virtual asteroids: elements, target plane points            
  INTEGER, PARAMETER :: npox=4000 
  DOUBLE PRECISION elm(6,npox),xcl(npox),ycl(npox) 
! confidence boundary, line of max variation                            
  INTEGER ibv,npo,npo1,npoc,nc 
  DOUBLE PRECISION sigma,maxsig,minsig 
! correspondent ellipse in orbital elements space                       
  DOUBLE PRECISION b(2,2),ceicel(4,2),v(6,6) 
  DOUBLE PRECISION xc(3),vc(3),dxl,dyl 
! index of planet (for comparison)                                      
  INTEGER iplamold 
! output units                                                          
  INTEGER iun7,iun8,iun9 
! device flag,labels for plots                                          
  INTEGER idev 
  CHARACTER*60 ylabel,xlabel,title 
! loop indexes                                                          
  INTEGER i,n 
!==========================================                             
! input specification of set of points                                  
  CALL asscbd(iunout,npox,npo,sigma,ibv) 
! If ibv=0 then use automatic selection method                          
  IF(ibv.eq.0)THEN 
     maxsig=max(sig(1,jcsel),sig(2,jcsel)) 
     minsig=min(sig(1,jcsel),sig(2,jcsel)) 
     if(maxsig/minsig.le.200.d0)then 
        ibv=1 
     else 
        ibv=2 
     endif
  endif
! compute ellipse in the elements space                                 
  CALL slinel(dtpdet(1,1,jcsel),unc0,ceicel,b,v) 
! compute line of orbital elements                                      
  CALL linobs(ibv,npo,el0,axes(1,1,jcsel),sig(1,jcsel),    &
     &        b,v,sigma,ceicel,elm,npo1)                                
! open output files for graphics                                        
  IF(inl.eq.1)THEN 
! linear analysis is by definition on a fixed TP                        
     tcla0=tcla(jcsel) 
     CALL filopn(iun7,'mtp.fla','unknown') 
  ELSEIF(inl.eq.2)THEN 
! semilinear analysys is on a variable TP, the close approach manifold  
     CALL filopn(iun8,'clo.fla','unknown') 
! also resonant return analysys could be done by looking at the elements
     CALL filopn(iun9,'a.fla','unknown') 
! need to find tolerances on target plane time                          
     v_inf=v_infty(el0) 
     tcla0=tcla(jcsel) 
     CALL aftclov(iplam,el0%t,tcla0,v_inf,tbefore,tafter) 
  ELSE 
     WRITE(*,*)'tpslin: option inl=',inl,' unknown' 
     RETURN 
  ENDIF
! ===========================================================           
! main loop on the number of output points                              
  nc=0 
  iplamold=iplam 
  DO 7 n=1,npo1 
     IF(inl.eq.1)THEN 
! linear map from ellipse                                               
        dxl=DOT_PRODUCT(elm(1:6,n),dtpdet(1:6,1,jcsel)) 
        dyl=DOT_PRODUCT(elm(1:6,n),dtpdet(1:6,2,jcsel)) 
        nc=nc+1 
        xcl(nc)=dxl+tpc(1,jcsel) 
        ycl(nc)=dyl+tpc(2,jcsel) 
! file output: target plane linear map                                  
        WRITE(iun7,107)xcl(nc),ycl(nc) 
107     FORMAT(7e20.12) 
        elm(1:6,n)=el0%coord+elm(1:6,n) 
     ELSEIF(inl.eq.2)THEN 
! full n-body propagation from ellipse                                  
        elm(1:6,n)=el0%coord+elm(1:6,n) 
        njc=0 
        iplam=0 
        elv=el0
        elv%coord=elm(1:6,n)
        CALL pro_ele(elv,tafter,el1) 
        WRITE(iun9,109)el1%coord 
109     FORMAT(6f20.15) 
! check for existence of data                                           
! PROBLEM: jcsel is always the same???                                  
        IF(njc.eq.0.or.iplam.ne.iplamold)THEN 
           WRITE(*,*)' no close approach to planet ',iplamold 
        ELSEIF(njc.lt.jcsel)THEN 
           WRITE(*,*)n,' multiple close approach',njc,jcsel 
        ELSEIF(tcla(jcsel).gt.max(tbefore,tafter).or.                &
     &          tcla(jcsel).lt.min(tbefore,tafter))THEN                 
           WRITE(*,130)n,tcla(jcsel), tbefore, tafter 
130        FORMAT(i4,' time_cla=',f11.2,' out of range, must be in ' &
     &             ,2f11.2)                                             
        ELSE 
           nc=nc+1 
! reference system with third axis normal to ecliptic, first=velocity   
! is available from tpcana.h
           xc=MATMUL(vt3(1:3,1:3,jcsel),xcla(1:3,jcsel)) 
           vc=MATMUL(vt3(1:3,1:3,jcsel),vcla(1:3,jcsel)) 
           xcl(nc)=xc(2) 
           ycl(nc)=xc(3) 
! file output: close approach manifold                                  
           WRITE(*,131)xc(3),xc(2),tcla(jcsel),n,nc 
131        FORMAT(f7.4,1x,f7.4,1x,f11.4,1x,i4,1x,i4) 
           WRITE(iun8,107)xc,vc,tcla(jcsel) 
        ENDIF
     ENDIF
7 ENDDO
  IF(inl.eq.1)THEN 
     CALL filclo(iun7,' ') 
  ELSEIF(inl.eq.2)THEN 
     CALL filclo(iun8,' ') 
     CALL filclo(iun9,' ') 
  ENDIF
  iplam=iplamold 
! ======================================================================
! graphics:                                                             
! close approach manifold number of points                              
  IF(ibv.eq.1.and.nc.gt.0)THEN 
     npoc=nc+1 
     xcl(npoc)=xcl(1) 
     ycl(npoc)=ycl(1) 
  ELSE 
     npoc=nc 
  ENDIF
  IF(npoc.eq.0)THEN 
     WRITE(*,*)' nothing to plot!' 
     RETURN 
  ENDIF
! labels                                                                
  IF(inl.eq.1)THEN 
     xlabel=' target plane xi (AU), linear approximation ' 
  ELSE 
     xlabel=' target plane xi (AU), semilinear approximation ' 
  ENDIF
  ylabel=' target plane zeta (AU) (proj, ecl. normal)' 
! selection of graphics device                                          
  WRITE(title,200)astna0,tcla0 
200 FORMAT(a18,' encounter at MJD ',f9.2) 
2 CALL getdev(idev) 
  IF(idev.eq.0)RETURN 
! If we are making .ps files then the MTP will be overwritten!          
  if(idev.eq.5.or.idev.eq.6)then 
     write(*,*)'The PostScript file ''giffv.ps'' is '//             &
     &        'about to be overwritten. If you wish to save it'//       &
     &        ' you should rename it before proceeding.'                
     pause 
  endif
! plot, mark Earth                                                      
  CALL plotob(xcl,ycl,0.d0,0.d0,npoc,xlabel,ylabel,title,idev,1) 
  GOTO 2 
END SUBROUTINE tpslin
! ======================================================================
!  VIRtual IMPactor ephemerides                                         
! ======================================================================
SUBROUTINE virimp(tpc,dtpdet,axes,sig,ceicel,b,v,el0,unc0,iunout)
  USE fund_const 
  USE pred_obs                                            
! ===============input=================                                 
! nominal initial conditions: epoch, elements  
  TYPE(orbit_elem), INTENT(IN) :: el0
! normal and covariance matrices                                        
  TYPE(orb_uncert), INTENT(IN) :: unc0
! output unit                                                           
  INTEGER, INTENT(IN) :: iunout 
! target plane coordinates: cartesian                                   
  DOUBLE PRECISION, INTENT(IN) :: tpc(2) 
! computation of target ellipse                                         
  DOUBLE PRECISION, INTENT(IN) :: axes(2,2),sig(2),dtpdet(6,2) 
! correspondent ellipse in orbital lements space                        
  DOUBLE PRECISION, INTENT(IN) :: b(2,2),ceicel(4,2),v(6,6) 
! ==============end interface=================                          
! option flag                                                           
  INTEGER ivir 
! sigma level across the LOV                                            
  DOUBLE PRECISION sigimp 
! safe distance from the Earth                                          
  DOUBLE PRECISION dsafe 
! rectangle, and its copy in the elements space                         
  DOUBLE PRECISION tpstr(2),tpwea(2) 
  TYPE(orbit_elem) els
  DOUBLE PRECISION elems(6,4),del2(2),del4(4),dels(6) 
  DOUBLE PRECISION delems(6),delemw(6),selems(6),selemw(6) 
! observation data (proper motion, elongation, distance) header removed
! times                                                                 
  DOUBLE PRECISION tmjd,tut,sec1,sec2 
  INTEGER mjd1,mjd2 
! call to seleph                                                        
  DOUBLE PRECISION  tut1,tdt1,tut2,tdt2,dt 
  INTEGER idsta 
  CHARACTER*3 scale 
! for preobs                                                            
  DOUBLE PRECISION alpha,delta,appmag 
! station code, observation type                                        
  INTEGER ids,iob1 
! integer indexes, lengths, functions, units                            
  INTEGER j,jj,ln,iuneph,iunvir 
! for outobc: covariance of the observations                            
  DOUBLE PRECISION gamad(2,2),axesky(2,2),sigsky(2) 
  DOUBLE PRECISION :: adot,ddot ! proper motion
  DOUBLE PRECISION :: pha,dis,dsun,elo,gallat ! phase, distance to Earth, distance to Sun
! ephemeris options                                                     
  CHARACTER*80 fields 
  DOUBLE PRECISION mass 
  CHARACTER*60 file 
! file name                                                             
  CHARACTER*80 titnam 
! confidence boundary, line of max variation (alpha, delta, app. magnitu
  INTEGER npo, ibv, npo1, inl 
  INTEGER, PARAMETER :: npox=4000 
  DOUBLE PRECISION sigma, al(npox),de(npox),appmagv(npox) 
  CHARACTER*1 type
! line of elements                                                      
  DOUBLE PRECISION elm(6,npox) 
! menu                                                       
  CHARACTER*20 menunam 
! =====================================================                 
! chase is open for the virtual impactor; what to do?                   
10 menunam='null' 
  CALL menu(ivir,menunam,3,'how to catch it?=',                     &
     &     'exploratory ephemerides=',                                  &
     &     'select observation time=',                                  &
     &     'output impact orbital elements=')
  IF(ivir.eq.0)RETURN 
  IF(ivir.eq.1)THEN 
! ======= GENERATE EPHEMERIS =========                                  
! select time interval, step                                            
     CALL seleph(tut1,tdt1,tut2,tdt2,dt,idsta) 
     CALL filnam('.','virimp','eph',file,ln) 
     CALL filopn(iuneph,file(1:ln),'unknown') 
     fields='cal,mjd,coord,mag,elong,glat,r,delta,appmot,skyerr' 
     scale='UTC' 
     IF(nint(abs(tdt2-tdt1)/dt).gt.500)THEN 
        write(*,*)'Too many ephemeris points: ',                    &
   &           nint(abs(tdt2-tdt1)/dt)                                
        write(*,*)'Select a time interval and time span to ',       &
     &        'ensure that there are fewer than 500 points.'            
        goto 10 
     ELSE 
        CALL ephemc(iuneph,el0,unc0,.true.,tdt1,tdt2,        &
     &           dt,idsta,scale,fields)                  
     ENDIF
     CALL filclo(iuneph,' ') 
     WRITE(*,*)' look at the ephemerides in file ./virimp.eph ' 
  ELSEIF(ivir.eq.3)THEN 
     CALL filopn(iunvir,'virimp.eq0','unknown') 
! output header                                                         
     CALL wromlh (iunvir,'ECLM','J2000') 
     CALL write_elems(el0,'virimp   ','ML',file,iunvir,unc0)
!    CALL wromlr (iunvir,'virimp',eqc,'EQU',tc,gc,.true.,           &
!     &        cc,.true.,hmag,gmag,0.d0)                                 
     CALL filclo(iunvir,' ') 
  ELSEIF(ivir.eq.2)THEN 
! ============PREDICT OBSERVATION==================                     
! ======two dimensional preimage===============                         
! compute rectangle on the MTP  enclosing all the impact points;        
! note that the target plane point tpc is used as origin                
     sigimp=5.d0 
     dsafe=5.d0*4.2e-5 
     DO j=1,2 
        tpstr(j)=axes(j,1)*sig(1)*sigimp 
        tpwea(j)=axes(j,2)*dsafe 
     ENDDO
! find corresponding orbital elements at epoch                          
     del2=MATMUL(b,tpstr) ! CALL mulmav(b,2,2,tpstr,2,del2) 
     del4=MATMUL(ceicel,del2) ! CALL mulmav(ceicel,4,2,del2,2,del4) 
     dels(1:2)=del2 ! CALL vcopy(2,del2,dels) 
     dels(3:6)=-del4(1:4) 
     delems=MATMUL(v,dels) ! CALL mulmav(v,6,6,dels,6,delems) 
! linear map from ellipse                                               
!        dxl=prscag(6,delems,dtpdet(1,1))                               
!        dyl=prscag(6,delems,dtpdet(1,2))                               
!        WRITE(*,*)dxl,dyl                                              
! find corresponding orbital elements at epoch                          
     del2=MATMUL(b,tpwea) ! CALL mulmav(b,2,2,tpwea,2,del2) 
     del4=MATMUL(ceicel,del2) ! CALL mulmav(ceicel,4,2,del2,2,del4) 
     dels(1:2)=del2 ! CALL vcopy(2,del2,dels) 
     DO jj=1,4 
        dels(2+jj)=-del4(jj) 
     ENDDO
     delemw=MATMUL(v,dels) ! CALL mulmav(v,6,6,dels,6,delemw) 
! linear map from ellipse                                               
!        dxl=prscag(6,delemw,dtpdet(1,1))                               
!        dyl=prscag(6,delemw,dtpdet(1,2))                               
!        WRITE(*,*)dxl,dyl                                              
! compute 4 corners                                                     
     elems(1:6,1)=el0%coord+delems
     elems(1:6,2)=el0%coord+delemw
     selemw=-delemw 
     selems=-delems 
     elems(1:6,1)=elems(1:6,1)+selemw
     elems(1:6,3)=el0%coord+selems
     elems(1:6,4)=elems(1:6,3)+delemw
     elems(1:6,3)=elems(1:6,3)+selemw
! give time                                                             
     WRITE(*,*)' give time for prediction (MJD)' 
     READ(*,*)tmjd 
! universal time of the required observation                            
     mjd1=FLOOR(tmjd) 
     sec1=(tmjd-float(mjd1))*86400.d0 
     CALL cnvtim(mjd1,sec1,'TDT',mjd2,sec2,'UTC') 
     tut=sec2/86400.d0+float(mjd2) 
     ids=500 
     iob1=1001 
     type='O'
! predict                           
     inl=1                                    
     CALL predic_obs(el0,ids,tmjd,type,      &
     &        alpha,delta,appmag,inl,                                    &
     &        UNCERT=unc0,GAMAD=gamad,SIG=sigsky,AXES=axesky,             &
     &        ADOT0=adot,DDOT0=ddot,DIS0=dis,PHA0=pha,DSUN0=dsun,       &
     &        ELO0=elo,GALLAT0=gallat)
     CALL outobc(iunout,type,ids,tut,alpha,delta,appmag,adot,ddot,   &
     &        elo,dis,2,gamad,sigsky,axesky)                            
     DO jj=1,4 
        els=el0
        els%coord=elems(1:6,jj)
        CALL  predic_obs(els,ids,tmjd,type,al(jj),de(jj),appmag,inl)
        al(jj)=al(jj)-alpha 
        de(jj)=de(jj)-delta 
!          write(*,*)jj,al(jj)*degrad,de(jj)*degrad                    
     ENDDO
! ======four dimensional preimage===============                        
     menunam='prednonl' 
     CALL menu(inl,menunam,3,'How to handle nonlinearity?=',        &
     &        'linear map=',                                            &
     &        '2-body nonlinearity=',                                   &
     &        'full n-body nonlinearity=')
     IF(inl.eq.0)RETURN 
! input specification of set of points                                  
     CALL asscbd(iunout,npox,npo,sigma,ibv) 
     CALL preob4(el0,ids,tmjd,unc0,                      &
     &        v,sigma,npo,ibv,inl,al(5),de(5),appmagv,elm,           &
     &        alpha,delta,appmag,gamad,sigsky,axesky,npo1)              
     al(5+npo1)=al(5) 
     de(5+npo1)=de(5) 
     DO j=1,npo1+5 
        write(*,*)'diff. observation ', j,al(j)*degrad,de(j)*degrad 
     ENDDO
     titnam='virtual impactor' 
     CALL plocbd(titnam,alpha,delta,5.d0,tut,al,de,5+npo1,iob1) 
  ENDIF 
  GOTO 10 
END SUBROUTINE virimp                                          
!                                        
END MODULE target_plane

! ================================================================      
! AFTCLOV                                                               
! select time interval to get after the close approach                  
! taking into account the relative velocity w.r. to Earth               
      SUBROUTINE aftclov(iplam,t0,tcla,v_inf,tbefore,tafter) 
      USE planet_masses
      IMPLICIT NONE 
! input: planet number, time of initial conditions,                     
! of known close approach, velocity                                     
      INTEGER iplam 
      DOUBLE PRECISION t0,tcla,v_inf 
! output: time "after", time "before" (depending upon sense of propagati
      DOUBLE PRECISION tafter,tbefore 
! end interface                                                         
! time interval to exit from TP disk                                    
      DOUBLE PRECISION delt_tp 
! target plane disk radius from planet_masses
! ===================================================================   
! warning: really done only for Earth                                   
      IF(ordnam(iplam).ne.'EARTH') THEN 
         WRITE(*,*)' aftclov: not to be used for planet ',ordnam(iplam) 
         delt_tp=180.d0 
! time interval to exit from TP disk
      ELSEIF(v_inf.gt.0.d0)THEN                                 
         delt_tp=2*dmin(iplam)/v_inf
      ELSE
         delt_tp=365.25d0
      ENDIF 
! forced to avoid infinte intervals for v-inf=0
      delt_tp=MIN(delt_tp,365.25d0)  
! forced to avoid short intervals for fast encounters                   
      IF(delt_tp.lt.50.d0)delt_tp=50.d0 
! time interval to be clear out of the TP disk is given as deltat       
      IF(t0.lt.tcla)THEN 
! future close approaches                                               
         tafter=tcla+delt_tp 
         tbefore=tcla-delt_tp 
      ELSE 
! past close approaches                                                 
         tafter=tcla-delt_tp 
         tbefore=tcla+delt_tp 
      ENDIF 
      RETURN 
      END SUBROUTINE aftclov 
! ===================================================
! velocity at infinity with repect to Earth
! (circular approx) in AU/day
DOUBLE PRECISION FUNCTION v_infty(el0)
  USE fund_const
  USE orbit_elements 
  IMPLICIT NONE 
! input elements
  TYPE(orbit_elem), INTENt(IN) :: el0
! END INTERFACE
  TYPE(orbit_elem) eleq
  DOUBLE PRECISION eq0(6) 
  DOUBLE PRECISION cosi,v2
  CHARACTER*3 coo
  INTEGER fail_flag
  coo=el0%coo
  CALL coo_cha(el0,'EQU',eleq,fail_flag)
  eq0=eleq%coord
  IF(fail_flag.eq.0)THEN 
! v_infinity computation in equinoctal
     cosi=(1.d0-eq0(4)**2-eq0(5)**2)/                                  &
     &     (1.d0+eq0(4)**2+eq0(5)**2)                                   
     v2=3.d0-1.d0/eq0(1) -                                             &
     &     2.d0*sqrt(eq0(1)*(1.d0-eq0(2)**2-eq0(3)**2))*cosi            
  ELSEIF(fail_flag.eq.5)THEN
! v_infinity computation in cometary
     cosi=cos(eq0(3))
     v2=3.d0-(1.d0-eq0(2))/eq0(1) - 2.d0*sqrt(eq0(1)*(1.d0+eq0(2)))*cosi
  ELSE
     WRITE(*,*)' v_infty: coordinate conversion failed ', el0,fail_flag  
  ENDIF
! handle non hyperbolic (w.r. to Earth) case
  IF(v2.gt.0.d0)THEN 
     v_infty=sqrt(v2) 
  ELSE 
     v_infty=0.d0 
  ENDIF
! normalization in AU/day by Gauss constant
  v_infty=v_infty*gk 
END FUNCTION v_infty

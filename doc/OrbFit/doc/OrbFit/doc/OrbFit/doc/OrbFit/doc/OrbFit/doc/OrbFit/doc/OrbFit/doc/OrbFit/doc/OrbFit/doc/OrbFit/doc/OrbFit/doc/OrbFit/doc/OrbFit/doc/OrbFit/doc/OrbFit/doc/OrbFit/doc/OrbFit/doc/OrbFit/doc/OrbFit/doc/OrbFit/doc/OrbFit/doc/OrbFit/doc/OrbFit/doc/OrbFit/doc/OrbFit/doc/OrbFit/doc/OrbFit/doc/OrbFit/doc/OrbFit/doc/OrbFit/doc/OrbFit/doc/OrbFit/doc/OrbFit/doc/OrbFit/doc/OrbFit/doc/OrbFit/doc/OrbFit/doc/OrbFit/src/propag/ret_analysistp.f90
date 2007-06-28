!==================MODULE ret_analysistp======================            
! Created by Giacomo Tommei, 03 January 2005                           
!===========================================================
! CONTAINS
! LIST OF PUBLIC ENTITIES:
!
!  SUBROUTINES 
!             ret_minctp
!
! PRIVATE ROUTINES: 
!                 achillestp
!                 falslog4tp
!  	          findminctp
!===========================================================
! DEP.
!ret_analysistp.o: \
!	../suit/fund_const.mod \
!	../suit/planet_masses.mod \
!	close_app.o \
!	tp_trace.o 
!===========================================================

MODULE ret_analysistp
  USE fund_const
  USE tp_trace

  IMPLICIT NONE
  PRIVATE

! LIST OF PUBLIC ENTITIES
! SUBROUTINEs
  PUBLIC :: ret_minctp, ret_min, arrcut
  PUBLIC :: wrireptp, header_rep, wriwarn, header_new, wriouttp

  DOUBLE PRECISION, PUBLIC :: beta_factor

! controls for methods
  LOGICAl, PUBLIC :: segment_method                                          
CONTAINS 

!============================================================================
! RET_MIN (public subroutine) FITOBS
!============================================================================
  SUBROUTINE ret_min(lre,vas_tr,t0,dnewton,no_risk)  
    USE offlov_checktp 
    USE output_control    
    USE virtual_impactor
!=======================INPUT================================================
    INTEGER, INTENT(IN) :: lre  ! no close approach records in filament  
    TYPE(tp_point), DIMENSION(lre),INTENT(IN) :: vas_tr 
    DOUBLE PRECISION, INTENT(IN) :: t0           ! epoch 
    DOUBLE PRECISION, INTENT(IN)   :: dnewton    ! control for falsi/newton 
!=======================OUTPUT===============================================
    INTEGER, INTENT(OUT) :: no_risk   ! number of risk cases found
    INTEGER iunnew,iunwarn,iunrisk    ! output units (redirected)
!=================== END INTERFACE===========================================
    DOUBLE PRECISION, DIMENSION(lre) :: dminpos
    INTEGER :: j 
    DOUBLE PRECISION, DIMENSION(lre) :: dr2ds 
    DOUBLE PRECISION :: dminpo
    DOUBLE PRECISION :: x1,x2 
    CHARACTER(LEN=7) :: type          ! output from falslog               
    TYPE(tp_point) :: va_tracemin      ! output from findminctp   
    DOUBLE PRECISION :: tmin,smin 
    LOGICAL :: fals_conv,fals_notp 
    INTEGER :: nvai 
    INTEGER :: niter,iunwarn0,iunnew0                                                         
    DOUBLE PRECISION :: siglim,dalpha1,dalpha2,pridif,dfalsi
    INTEGER :: i1,i2 
    DOUBLE PRECISION, PARAMETER :: del_fal=2.3d-3    ! convergence control 
                  ! this is the limit change in distance in RE (=1e-7 AU)  
!=============================================================================
    dminpos=vas_tr%minposs
    no_risk=0    
    num_vi=0
    iunnew= -iun_log ! standard output + logfile 
    iunwarn= -iun_log ! standard output + logfile 
    iunrisk= -iun_log ! standard output + logfile 
    iunnew0=iun_log
    iunwarn0=iun_log
! derivative of distance squared with respect to sigma                  
    dr2ds(1:lre)=vas_tr(1:lre)%dd2_ds 
! scan filament                                                         
    i1=vas_tr(1)%rindex 
    i2=vas_tr(lre)%rindex 
    DO 1 j=1,lre+1
       IF(j.eq.1)THEN 
! handle heads (including singletons)                                   
          IF(dr2ds(1).gt.0.d0)THEN 
             IF(lre.eq.1)THEN 
                type='SINGLET' 
             ELSE 
                type='ASCHEAD' 
             ENDIF 
! increasing head                                                       
             IF(dminpos(1).lt.dnewton)THEN 
                siglim=deltasig/2.d0 
                CALL achillestp(vas_tr(1),siglim,type,        &
     &                iunwarn,iunnew,va_tracemin,niter,fals_conv,fals_notp)
                IF(fals_notp)THEN 
                   WRITE(iunnew0,*)' achillestp no TP ',type 
                ELSE 
                   CALL riskchecktp(va_tracemin,t0,type,no_risk,iunnew,iunwarn,iunrisk)
                   CALL store_vi()
                ENDIF
             ENDIF
          ENDIF 
       ELSEIF(j.eq.lre+1)THEN 
! handle tails (including singletons)                                   
          IF(dr2ds(lre).lt.0.d0)THEN 
             IF(lre.eq.1)THEN 
                type='SINGLET' 
             ELSE 
                type='DESTAIL' 
             ENDIF 
! decreasing tail                                                       
             IF(dminpos(lre).lt.dnewton)THEN 
                siglim=deltasig/2.d0 
                CALL achillestp(vas_tr(lre),siglim,type,      &
     &                iunwarn,iunnew,va_tracemin,niter,fals_conv,fals_notp)
! if failure write failure record, else possibly risk record            
                IF(fals_notp)THEN 
                   WRITE(iunnew0,*)' achillestp no TP ',type 
                ELSE 
                   CALL riskchecktp(va_tracemin,t0,type,no_risk,iunnew,iunwarn,iunrisk)
                   CALL store_vi()
                ENDIF
             ENDIF
          ENDIF 
       ELSE 
          x2=vas_tr(j)%rindex 
          x1=vas_tr(j-1)%rindex 
          CALL falslog4tp(vas_tr(j-1),iunwarn,type) 
! check for interesting cases                                           
          dminpo=min(dminpos(j),dminpos(j-1)) 
          IF(dminpo.lt.dnewton)THEN 
!  if falslog proposes to try regula falsi, and the MOID is low, do it  
             IF(type.eq.'SIMPMIN'.or.type.eq.'INTEMIN')THEN 
                CALL findminctp(vas_tr(j-1),x1,x2,type,                  &
       &                iunwarn,iunnew,va_tracemin,                  &
     &                fals_conv,niter,fals_notp)  
! if failure write failure record, else possibly risk record            
                IF(fals_notp)THEN 
                   WRITE(iunnew0,*)' falsi no TP ' 
                   siglim=deltasig/2.d0 
                   CALL achillestp(vas_tr(j-1),siglim,     &
     &                      type,iunwarn,iunnew,va_tracemin,               &
     &                      niter,fals_conv,fals_notp)                  
                   IF(fals_notp.or..not.fals_conv)THEN 
                      WRITE(iunnew0,*)' achillestp no TP/conv ',type 
                   ELSE 
                      CALL riskchecktp(va_tracemin,t0,type,no_risk,iunnew,iunwarn,iunrisk)
                      CALL store_vi()
                   ENDIF
                   siglim=deltasig/2.d0 
                   CALL achillestp(vas_tr(j),siglim,type,iunwarn,iunnew,va_tracemin,    &
     &                      niter,fals_conv,fals_notp)                  
                   IF(fals_notp.or..not.fals_conv)THEN 
                      WRITE(iunnew0,*)' achillestp no TP/conv ',type 
                   ELSE 
                      CALL riskchecktp(va_tracemin,t0,type,no_risk,iunnew,iunwarn,iunrisk) 
                      CALL store_vi()
                   ENDIF
                ELSE 
                      CALL riskchecktp(va_tracemin,t0,type,no_risk,iunnew,iunwarn,iunrisk)
                      CALL store_vi()
                ENDIF
! but for the interrupted minima, there is always the suspicion that the
! composite, that is two minima and one maximum                         
                IF(type.eq.'INTEMIN')THEN 
! change in direction of the LOV                                        
                   dalpha1=pridif(vas_tr(j-1)%alpha_lov,va_tracemin%alpha_lov) 
                   dalpha2=pridif(vas_tr(j)%alpha_lov,va_tracemin%alpha_lov)
! store the minimum distance output from falsi                          
                   dfalsi=va_tracemin%b 
                   WRITE(iunwarn0,*)'check for almost interrupted',     &
     &                   dalpha1*degrad,dalpha2*degrad                  
                   IF(cos(dalpha2).gt.0.70d0)THEN 
                      siglim=deltasig*abs(vas_tr(j-1)%rindex-va_tracemin%rindex)/2.d0 
                      CALL achillestp(vas_tr(j-1),siglim,     &
     &                      type,iunwarn,iunnew,va_tracemin,               &
     &                      niter,fals_conv,fals_notp)                  
                      IF(fals_notp.or..not.fals_conv)THEN 
                         WRITE(iunnew0,*)' achillestp no TP/conv ',type 
                      ELSE 
                         IF(abs(dfalsi-va_tracemin%b).gt.del_fal*dfalsi)   &
        &                   CALL riskchecktp(va_tracemin,t0,type,no_risk,iunnew,iunwarn,iunrisk) 
                         CALL store_vi()
                      ENDIF
                   ELSEIF(cos(dalpha1).gt.0.7d0)THEN 
!                       x1=vas_tr(j)%rindex 
                      siglim=deltasig*abs(vas_tr(j)%rindex-va_tracemin%rindex)/2.d0 
                      CALL achillestp(vas_tr(j),siglim,        &
     &                      type,iunwarn,iunnew,va_tracemin,               &
     &                      niter,fals_conv,fals_notp)                  
                      IF(fals_notp.or..not.fals_conv)THEN 
                         WRITE(iunnew0,*)' achillestp no TP/conv ', type 
                      ELSE 
                         IF(abs(dfalsi-va_tracemin%b).gt.del_fal*dfalsi)   &
        &                   CALL riskchecktp(va_tracemin,t0,type,no_risk,iunnew,iunwarn,iunrisk) 
                         CALL store_vi()
                      ENDIF
                   ENDIF
                                                                        
                ENDIF
             ELSEIF(type.eq.'INTEUNK')THEN 
                IF(dr2ds(j-1).lt.0.d0)THEN 
                   siglim=deltasig/2.d0 
                   CALL achillestp(vas_tr(j-1),siglim,type,   &
     &                iunwarn,iunnew,va_tracemin,niter,fals_conv,fals_notp)
                ELSEIF(dr2ds(j).gt.0.d0)THEN 
                   siglim=deltasig/2.d0 
                   CALL achillestp(vas_tr(j),siglim,type,     &
     &                iunwarn,iunnew,va_tracemin,niter,fals_conv,fals_notp)
                ELSE 
                   WRITE(iunnew0,*)' case not understood ', type 
                ENDIF
                IF(fals_notp)THEN 
                   WRITE(iunnew0,*)' achillestp no TP ',type 
                ELSE 
                   CALL riskchecktp(va_tracemin,t0,type,no_risk,iunnew,iunwarn,iunrisk)
                   CALL store_vi()
                ENDIF
             ELSEIF(type.eq.'ENTADES')THEN 
! look for first minimum                                                
                siglim=deltasig/2.d0 
                CALL achillestp(vas_tr(j-1),siglim,type,      &
     &                iunwarn,iunnew,va_tracemin,niter,fals_conv,fals_notp)
                IF(fals_notp)THEN 
                   WRITE(iunnew0,*)' achillestp no TP ',type 
                ELSE 
                   CALL riskchecktp(va_tracemin,t0,type,no_risk,iunnew,iunwarn,iunrisk)
                   CALL store_vi()
                ENDIF
! missing second half of the entagled case                              
                WRITE(iunnew0,*)type, ' second minimum not handled' 
             ELSEIF(type.eq.'ENTAASC')THEN 
! look for last minimum                                                 
                siglim=deltasig/2.d0 
                CALL achillestp(vas_tr(j),siglim,type,        &
     &                iunwarn,iunnew,va_tracemin,niter,fals_conv,fals_notp)
                IF(fals_notp)THEN 
                   WRITE(iunnew0,*)' achillestp no TP ',type 
                ELSE 
                   CALL riskchecktp(va_tracemin,t0,type,no_risk,iunnew,iunwarn,iunrisk)               
                   CALL store_vi()
                ENDIF
! missing second half of the entagled case                              
                WRITE(iunnew0,*)type, ' second minimum not handled' 
             ELSE 
                IF(type.eq.'NOWORRY'.or.type.eq.'INTEMAX'.or.          &
     &                   type.eq.'SIMPMAX')THEN                         
!  ok not to be handled                                                 
                ELSE 
                   WRITE(iunnew0,*)' case not handled, type=',type,     &
   &                   vas_tr(j-1)%tcla, vas_tr(j-1)%rindex    
                   WRITE(iunwarn0,*)' case not handled, type=',type 
                   CALL wriouttp(vas_tr(j-1),iunwarn) 
                   CALL wriouttp(vas_tr(j),iunwarn) 
                ENDIF
             ENDIF
          ENDIF
       ENDIF
                                                                        
1   ENDDO 

  END SUBROUTINE ret_min              
!============================================================================
! RET_MINCTP (public subroutine) RESRET2TP
!============================================================================
  SUBROUTINE ret_minctp(lre,vas_tr,t0,iunwarn,iunnew,riskfile,dnewton,no_risk)
    USE offlov_checktp     
!=======================INPUT================================================
    INTEGER, INTENT(IN) :: lre                   ! number of close approach
                                                   ! records in filament  
    TYPE(tp_point), DIMENSION(lre), INTENT(IN) :: vas_tr 
    DOUBLE PRECISION, INTENT(IN) :: t0           ! epoch 
    DOUBLE PRECISION, INTENT(IN)   :: dnewton     ! control for falsi/newton 
!=====================OUTPUT UNITS===========================================
    INTEGER            :: iunwarn,iunnew 
    CHARACTER(LEN=100) :: riskfile  
!=======================OUTPUT===============================================
    INTEGER, INTENT(OUT) :: no_risk 
!=================== END INTERFACE===========================================
    DOUBLE PRECISION, DIMENSION(lre) :: dminpos
    INTEGER :: j 
    DOUBLE PRECISION, DIMENSION(lre) :: dr2ds 
    DOUBLE PRECISION :: dminpo
    DOUBLE PRECISION :: x1,x2 
    CHARACTER(LEN=7) :: typ          ! output from falslog               
    TYPE(tp_point) :: va_tracemin      ! output from findminctp   
    DOUBLE PRECISION :: tmin,smin 
    LOGICAL :: fals_conv,fals_notp 
    INTEGER :: nbigrisk, nsegrisk
    INTEGER :: nvai 
    INTEGER :: niter                                                         
    DOUBLE PRECISION :: siglim,dalpha1,dalpha2,pridif,dfalsi
    INTEGER :: i1,i2 
    DOUBLE PRECISION, PARAMETER :: del_fal=2.3d-3    ! convergence control 
                      ! this is the limit change in distance in RE (=1e-7 AU)  
!============================================================================
    dminpos=vas_tr%minposs
    no_risk=0
! derivative of distance squared with respect to sigma
    dr2ds(1:lre)=vas_tr(1:lre)%dd2_ds 
! scan filament                                                         
    i1=vas_tr(1)%rindex 
    i2=vas_tr(lre)%rindex 
    DO 1 j=1,lre+1
       IF(j.eq.1)THEN 
! handle heads (including singletons)                                   
          IF(dr2ds(1).gt.0.d0)THEN 
             IF(lre.eq.1)THEN 
                typ='SINGLET' 
             ELSE 
                typ='ASCHEAD'
             ENDIF
! increasing head                                                       
             IF(dminpos(1).lt.dnewton)THEN 
                siglim=deltasig/2.d0 
                CALL achillestp(vas_tr(1),siglim,typ,        &
     &                iunwarn,iunnew,va_tracemin,niter,fals_conv,fals_notp)
                IF(fals_notp)THEN 
                   WRITE(iunnew,*)' achillestp no TP ',typ 
                ELSE 
                   CALL riskchecktp(va_tracemin,t0,typ,no_risk,    &
     &                   iunnew,iunwarn,RISKFILE=riskfile)
                ENDIF
             ENDIF
          ENDIF
       ELSEIF(j.eq.lre+1)THEN 
! handle tails (including singletons)                                   
          IF(dr2ds(lre).lt.0.d0)THEN 
             IF(lre.eq.1)THEN 
                typ='SINGLET' 
             ELSE 
                typ='DESTAIL' 
             ENDIF 
! decreasing tail                                                       
             IF(dminpos(lre).lt.dnewton)THEN 
                siglim=deltasig/2.d0 
                CALL achillestp(vas_tr(lre),siglim,typ,      &
     &                iunwarn,iunnew,va_tracemin,niter,fals_conv,fals_notp)
! if failure write failure record, else possibly risk record            
                IF(fals_notp)THEN 
                   WRITE(iunnew,*)' achillestp no TP ',typ 
                ELSE 
                   CALL riskchecktp(va_tracemin,t0,typ,no_risk,    &
     &                   iunnew,iunwarn,RISKFILE=riskfile)
                ENDIF
             ENDIF
          ENDIF
       ELSE 
          x2=vas_tr(j)%rindex 
          x1=vas_tr(j-1)%rindex 
          CALL falslog4tp(vas_tr(j-1),iunwarn,typ) 
! check for interesting cases                                           
          dminpo=min(dminpos(j),dminpos(j-1)) 
          IF(dminpo.lt.dnewton)THEN 
!  if falslog proposes to try regula falsi, and the MOID is low, do it  
             IF(typ.eq.'SIMPMIN'.or.typ.eq.'INTEMIN')THEN 
                CALL findminctp(vas_tr(j-1),x1,x2,typ,                  &
       &                iunwarn,iunnew,va_tracemin,                  &
     &                fals_conv,niter,fals_notp)                        
! if failure write failure record, else possibly risk record            
                IF(fals_notp)THEN 
                   WRITE(iunnew,*)' falsi no TP ' 
                   siglim=deltasig/2.d0 
                   CALL achillestp(vas_tr(j-1),siglim,     &
     &                      typ,iunwarn,iunnew,va_tracemin,               &
     &                      niter,fals_conv,fals_notp)                  
                   IF(fals_notp.or..not.fals_conv)THEN 
                      WRITE(iunnew,*)' achillestp no TP/conv ',typ 
                   ELSE 
                      CALL riskchecktp(va_tracemin,t0,typ,      &
     &                         no_risk,iunnew,iunwarn,RISKFILE=riskfile)
                   ENDIF
                   siglim=deltasig/2.d0 
                   CALL achillestp(vas_tr(j),siglim,       &
     &                      typ,iunwarn,iunnew,va_tracemin,               &
     &                      niter,fals_conv,fals_notp)                  
                   IF(fals_notp.or..not.fals_conv)THEN 
                      WRITE(iunnew,*)' achillestp no TP/conv ',typ 
                   ELSE 
                      CALL riskchecktp(va_tracemin,t0,typ,      &
     &                         no_risk,iunnew,iunwarn,RISKFILE=riskfile)
                   ENDIF
                ELSE 
                   nbigrisk=0
!                   CALL big_vitp(vas_tr,lre,va_tracemin,t0,nvai,nbigrisk,iunnew,iunwarn,RISKFILE=riskfile)
!                   IF(nbigrisk.gt.0)THEN
!                      no_risk=no_risk+nbigrisk
!                   ELSE
!                      IF(segment_method)THEN
!                         IF(nvai.eq.0)CALL risk_segmenttp(vas_tr(j-1),va_tracemin,t0,   &   !&                      nsegrisk,iunnew,iunwarn,riskfile)
!                         no_risk=no_risk+nsegrisk
!                      ENDIF
                   CALL riskchecktp(va_tracemin,t0,typ,no_risk,    &
     &                      iunnew,iunwarn,RISKFILE=riskfile)
!                   ENDIF                       
                ENDIF
! but for the interrupted minima, there is always the suspicion that the
! composite, that is two minima and one maximum                         
                IF(typ.eq.'INTEMIN')THEN 
! change in direction of the LOV                                        
                   dalpha1=pridif(vas_tr(j-1)%alpha_lov,va_tracemin%alpha_lov) 
                   dalpha2=pridif(vas_tr(j)%alpha_lov,va_tracemin%alpha_lov)
! store the minimum distance output from falsi                          
                   dfalsi=va_tracemin%b 
                   WRITE(iunwarn,*)'check for almost interrupted',     &
     &                   dalpha1*degrad,dalpha2*degrad                  
                   IF(cos(dalpha2).gt.0.70d0)THEN 
                      siglim=deltasig*abs(vas_tr(j-1)%rindex-va_tracemin%rindex)/2.d0 
                      CALL achillestp(vas_tr(j-1),siglim,     &
     &                      typ,iunwarn,iunnew,va_tracemin,               &
     &                      niter,fals_conv,fals_notp)                  
                      IF(fals_notp.or..not.fals_conv)THEN 
                         WRITE(iunnew,*)' achillestp no TP/conv ',typ 
                      ELSE 
                         IF(abs(dfalsi-va_tracemin%b).gt.del_fal*dfalsi)   &
        &                         CALL riskchecktp(va_tracemin,t0,typ, &
     &                         no_risk,iunnew,iunwarn,RISKFILE=riskfile)
                      ENDIF
                   ELSEIF(cos(dalpha1).gt.0.7d0)THEN 
!                       x1=vas_tr(j)%rindex 
                      siglim=deltasig*abs(vas_tr(j)%rindex-va_tracemin%rindex)/2.d0 
                      CALL achillestp(vas_tr(j),siglim,        &
     &                      typ,iunwarn,iunnew,va_tracemin,               &
     &                      niter,fals_conv,fals_notp)                  
                      IF(fals_notp.or..not.fals_conv)THEN 
                         WRITE(iunnew,*)' achillestp no TP/conv ', typ 
                      ELSE 
                         IF(abs(dfalsi-va_tracemin%b).gt.del_fal*dfalsi)   &
        &                         CALL riskchecktp(va_tracemin,t0,typ, &
     &                         no_risk,iunnew,iunwarn,RISKFILE=riskfile)
                      ENDIF
                   ENDIF                                                     
                ENDIF
             ELSEIF(typ.eq.'INTEUNK')THEN 
                IF(dr2ds(j-1).lt.0.d0)THEN 
                   siglim=deltasig/2.d0 
                   CALL achillestp(vas_tr(j-1),siglim,typ,   &
     &                iunwarn,iunnew,va_tracemin,niter,fals_conv,fals_notp)
                ELSEIF(dr2ds(j).gt.0.d0)THEN 
                   siglim=deltasig/2.d0 
                   CALL achillestp(vas_tr(j),siglim,typ,     &
     &                iunwarn,iunnew,va_tracemin,niter,fals_conv,fals_notp)
                ELSE 
                   WRITE(iunnew,*)' case not understood ', typ 
                ENDIF
                IF(fals_notp)THEN 
                   WRITE(iunnew,*)' achillestp no TP ',typ 
                ELSE 
                   CALL riskchecktp(va_tracemin,t0,typ,            &
     &                   no_risk,iunnew,iunwarn,RISKFILE=riskfile)
                ENDIF
             ELSEIF(typ.eq.'ENTADES')THEN 
! look for first minimum                                                
                siglim=deltasig/2.d0 
                CALL achillestp(vas_tr(j-1),siglim,typ,      &
     &                iunwarn,iunnew,va_tracemin,niter,fals_conv,fals_notp)
                IF(fals_notp)THEN 
                   WRITE(iunnew,*)' achillestp no TP ',typ 
                ELSE 
                   CALL riskchecktp(va_tracemin,t0,typ,            &
     &                   no_risk,iunnew,iunwarn,RISKFILE=riskfile)
                ENDIF
! missing second half of the entagled case                              
                WRITE(iunnew,*)typ, ' second minimum not handled' 
             ELSEIF(typ.eq.'ENTAASC')THEN 
! look for last minimum                                                 
                siglim=deltasig/2.d0 
                CALL achillestp(vas_tr(j),siglim,typ,        &
     &                iunwarn,iunnew,va_tracemin,niter,fals_conv,fals_notp)
                IF(fals_notp)THEN 
                   WRITE(iunnew,*)' achillestp no TP ',typ 
                ELSE 
                   CALL riskchecktp(va_tracemin,t0,typ,            &
     &                   no_risk,iunnew,iunwarn,RISKFILE=riskfile)
                ENDIF
! missing second half of the entagled case                              
                WRITE(iunnew,*)typ, ' second minimum not handled' 
             ELSE 
                IF(typ.eq.'NOWORRY'.or.typ.eq.'INTEMAX'.or.          &
     &                   typ.eq.'SIMPMAX')THEN                         
!  ok not to be handled                                                 
                ELSE 
                   WRITE(iunnew,*)' case not handled, typ=',typ,     &
   &                   vas_tr(j-1)%tcla, vas_tr(j-1)%rindex    
                   WRITE(iunwarn,*)' case not handled, typ=',typ 
                   CALL wriouttp(vas_tr(j-1),iunwarn) 
                   CALL wriouttp(vas_tr(j),iunwarn) 
                ENDIF
             ENDIF
          ENDIF
       ENDIF
1   ENDDO 

  END SUBROUTINE ret_minctp
!============================================================================
  SUBROUTINE falslog4tp(vas_tracec,iunwar,type) 
!===================== INPUT=========================================
    TYPE(tp_point), DIMENSION(2), INTENT(IN) :: vas_tracec    
                                                 ! 2 close approach records
!====================OUTPUT UNIT=====================================
    INTEGER :: iunwar 
!===================== OUTPUT========================================
!  ENCOUNTER TYPE                                                        
!  SIMPMIN= simple minimum (not interrupted, presumably single minimum) 
!  SIMPMAX= simple maximum (not interrupted, presumably single maximum) 
!  NOWORRY= increasing or decreasing (not interrupted) NO               
!  ENTAASC= entangled, multiple stationary points, first stationary poin
!  ENTADES= entangled, multiple stationary points, first stationary poin
!  INTEMIN= simple interrupted; at least one minimum in between, possibl
!  INTEMAX= simple maximum (interrupted, presumably single maximum) NO  
!  INTEUNK= interrupted, but not expected  DENSIFY                      
!  UNKNOWN= LOV turning (change of direction between 60 and 90 deg) DENS
    CHARACTER(LEN=7) :: type 
!================== END INTERFACE=======================================
    DOUBLE PRECISION :: dalpha,dbeta1,dbeta2,num,den,pridif 
    DOUBLE PRECISION :: betalim1,betalim2,gk,factor,sin20
!======================================================================
! change in direction of the LOV 
    IF(((vas_tracec(2)%opik%coord(4)-vas_tracec(1)%opik%coord(4)))==0.d0.and.  &
     &  (vas_tracec(2)%opik%coord(5)-vas_tracec(1)%opik%coord(5))== 0.d0 )THEN
       type='NOWORRY'
    ELSE                                
       dalpha=pridif(vas_tracec(1)%alpha_lov,vas_tracec(2)%alpha_lov) 
! angle between tangent to LOV and vector between two TP points         
       num=(vas_tracec(2)%opik%coord(4)-vas_tracec(1)%opik%coord(4))*         &
     &    (-vas_tracec(1)%stretch_lov*sin(vas_tracec(1)%alpha_lov))+        &
     &     (vas_tracec(2)%opik%coord(5)-vas_tracec(1)%opik%coord(5))         &
     &    *(vas_tracec(1)%stretch_lov*cos(vas_tracec(1)%alpha_lov))         
       den=vas_tracec(1)%stretch_lov*           &
     &     sqrt((vas_tracec(2)%opik%coord(4)-vas_tracec(1)%opik%coord(4))**2+                   &
     &              (vas_tracec(2)%opik%coord(5)-vas_tracec(1)%opik%coord(5))**2)  
       dbeta1=acos(num/den) 
       num=(vas_tracec(2)%opik%coord(4)-vas_tracec(1)%opik%coord(4))      &
     &    *(-vas_tracec(2)%stretch_lov*sin(vas_tracec(2)%alpha_lov))+        &
     &     (vas_tracec(2)%opik%coord(5)-vas_tracec(1)%opik%coord(5))       &
     &    *(vas_tracec(2)%stretch_lov*cos(vas_tracec(2)%alpha_lov))         
       den=vas_tracec(2)%stretch_lov*           &
     &     sqrt((vas_tracec(2)%opik%coord(4)-vas_tracec(1)%opik%coord(4))**2+                   &
     &     (vas_tracec(2)%opik%coord(5)-vas_tracec(1)%opik%coord(5))**2)       
       dbeta2=acos(num/den) 
! selection of type                                                     
       factor=beta_factor
       gk=0.01720209895d0 
       sin20=0.34202014d0
       betalim1=MAX(factor*gk*vas_tracec(1)%b/vas_tracec(1)%opik%coord(1),sin20) 
       betalim2=MAX(factor*gk*vas_tracec(2)%b/vas_tracec(2)%opik%coord(1),sin20)
       IF(cos(dalpha).gt.0.d0)THEN 
          IF(abs(sin(dbeta1)).lt.betalim1.or.abs(sin(dbeta2)).lt.betalim2)THEN
             IF(vas_tracec(1)%dd2_ds.lt.0.d0.and.vas_tracec(2)%dd2_ds.gt.0.0d0)THEN 
! simple minimum                                                        
                type='SIMPMIN' 
             ELSEIF(vas_tracec(1)%dd2_ds.gt.0.d0.and.vas_tracec(2)%dd2_ds.lt.0.0d0)THEN 
                type='SIMPMAX' 
             ELSE 
! could be increasing or decreasing, we do not care                     
                type='NOWORRY' 
             ENDIF
          ELSE 
! this is a complicated case, to be further discriminated               
             IF(vas_tracec(2)%dd2_ds.gt.0.d0)THEN 
                type='ENTAASC' 
                IF(vas_tracec(1)%dd2_ds.lt.0.d0)THEN 
                   WRITE(iunwar,*)' ENTAASC but first decreasing' 
                ENDIF
             ELSE 
                type='ENTADES' 
                IF(vas_tracec(1)%dd2_ds.gt.0.d0)THEN 
                   WRITE(iunwar,*)' ENTADES but second increasing' 
                ENDIF
             ENDIF
          ENDIF
       ELSEIF(cos(dalpha).lt.0.d0)THEN 
          IF(vas_tracec(1)%dd2_ds.lt.0.d0.and.vas_tracec(2)%dd2_ds.gt.0.0d0)THEN 
! there is at least one minimum of the interrupted type in between      
             type='INTEMIN' 
          ELSEIF(vas_tracec(1)%dd2_ds.gt.0.d0.and.vas_tracec(2)%dd2_ds.lt.0.0d0)THEN 
! there is at least one maximum of the almost interrupted type          
             type='INTEMAX' 
          ELSE 
! interrputed, but distance derivative does not change???               
! could have one maximum and one minimum                                
             type='INTEUNK' 
          ENDIF
       ELSE 
! case abolished                                                        
! dalpha is between 90 and 90 degrees, can this happen?                 
          type='UNKNOWN' 
       ENDIF
    ENDIF
  END SUBROUTINE falslog4tp
!===================================================================================
  SUBROUTINE achillestp(va_trace,siglim,type,                  &
     &     iunwar0,iunnew0,va_tracemin,niter,fals_conv,fals_notp) 
    USE output_control            
!========================== INPUT================================
    TYPE(tp_point), INTENT(IN) :: va_trace
    DOUBLE PRECISION, INTENT(IN) :: siglim 
    CHARACTER(LEN=7), INTENT(INOUT) :: type 
    INTEGER, INTENT(IN) :: iunwar0,iunnew0 
!========================== OUTPUT===============================     
    TYPE(tp_point), INTENT(OUT) :: va_tracemin                
    LOGICAL, INTENT(OUT) :: fals_notp,fals_conv 
    INTEGER, INTENT(OUT) :: niter 
!========================END INTERFACE===========================
    CHARACTER(LEN=4) :: method 
    DOUBLE PRECISION :: strlov,alphalov,csi,zeta,tc 
    INTEGER :: j,it,nn, iunwar,iunnew 
    INTEGER, PARAMETER :: nit=8 
    DOUBLE PRECISION, DIMENSION(nit+1) :: xx
    TYPE(tp_point), DIMENSION(nit+1) :: ar 
    DOUBLE PRECISION, DIMENSION(nit+1) ::f,g,ddx,xlim,d 
    TYPE(tp_point), DIMENSION(2) :: va2_trace       ! for call to findminctp 
    DOUBLE PRECISION :: dddx,cont_dist,cont_prob,fm,dds,           &
     &     ddw,cont_dist1,cont_dist2               ! convergence control      
    DOUBLE PRECISION, PARAMETER :: eps_fal=1.d-10    ! this is the limit probability 
                                                       ! (with uniform prob.density)     
    DOUBLE PRECISION, PARAMETER :: del_fal=2.3d-3    ! this is the limit change 
                                                       ! in distance in RE (=1e-7 AU)
!=================================================================
    WRITE(*,*)' achillestp entry ', type,va_trace%tcla,va_trace%rindex 
    iunwar=abs(iunwar0)
    iunnew=abs(iunnew0)
    WRITE(iunwar,*)' beginning Achillestp: ' 
    WRITE(iunwar,*)' x,t,f',va_trace%rindex,va_trace%tcla,va_trace%dd2_ds 
    ar(1)=va_trace
    xlim(1)=siglim/deltasig 
    xx(1)=va_trace%rindex 
! iterations                                                            
    nn=0 
    DO 1 it=1,nit 
       nn=nn+1 
! computation of f,g; uses the stretching as given in input, the LOV one
       strlov=ar(nn)%stretch_lov 
       alphalov=ar(nn)%alpha_lov 
       csi=ar(nn)%opik%coord(4) 
       zeta=ar(nn)%opik%coord(5) 
       tc=ar(nn)%tcla 
       d(nn)=ar(nn)%b
       f(nn)=ar(nn)%dd2_ds 
       g(nn)=2.d0*strlov**2 
! newton method                                                         
       dds=-f(nn)/g(nn) 
       ddx(nn)=dds/deltasig 
       IF(abs(ddx(nn)).gt.xlim(nn))THEN 
          ddx(nn)=ddx(nn)*abs(xlim(nn)/ddx(nn)) 
       ENDIF
       xx(nn+1)=xx(nn)+ddx(nn) 
       WRITE(iunwar,*)nn,ddx(nn),xlim(nn) 
! find and propagate to close approach the interpolated orbit
       CALL lovclosapptp(xx(nn+1),tc,tc,iunwar,fals_notp,ar(nn+1))
       IF(fals_notp)THEN 
          WRITE(iunwar,*)' achillestp: no TP',ddx(nn),xx(nn),it 
          WRITE(*,*)' achillestp: no TP',ddx(nn),xx(nn),it 
          nn=nn-1 
          method='ANTP' 
          niter=it 
          CALL  wriwarn(ar(nn+1)%tcla,ar(nn+1)%rindex,ar(nn+1)%b, &
     &           ar(nn+1)%dd2_ds,                                       &
     &           ar(nn+1)%stretch,ar(nn+1)%width,ar(nn+1)%moid,       &
     &           cont_dist,cont_prob,niter,type,method,iunwar)          
          fals_conv=.false. 
          xlim(nn+1)=min(xlim(nn+1)/2.d0,abs(ddx(nn+1)/2.d0)) 
       ELSE 
! convergence criteria                                                  
! intervals for control                                                 
          fm=abs(ar(nn+1)%dd2_ds) 
          d(nn+1)=ar(nn+1)%b 
          dddx=abs(xx(nn+1)-xx(nn)) 
          cont_dist1=fm*deltasig*dddx/(d(nn+1)+d(nn)) 
          cont_dist2=abs(d(nn+1)-d(nn)) 
          IF(it.eq.1)THEN 
             cont_dist=cont_dist1 
          ELSE 
             cont_dist=min(cont_dist1,cont_dist2) 
          ENDIF
          cont_prob=dddx*deltasig/6 
! write warn record                                                     
          method='Ac' 
          niter=nn 
          CALL  wriwarn(ar(nn+1)%tcla,ar(nn+1)%rindex,ar(nn+1)%b, &
               &           ar(nn+1)%dd2_ds,                                       &
     &           ar(nn+1)%stretch,ar(nn+1)%width,ar(nn+1)%moid,       &
     &           cont_dist,cont_prob,niter,type,method,iunwar0)            
          ddw=d(nn+1) 
          IF(cont_dist.lt.del_fal.or.cont_dist.lt.del_fal*ddw         &
     &           .or.cont_prob.lt.eps_fal.or.ddw.lt.1.d0)THEN           
! convergence achieved                                                  
             WRITE(iunwar,*)' minimum found, iter= ',it,nn 
             WRITE(*,*)' minimum found, iter= ',it,nn 
! write new record                                                      
             CALL  wriwarn(ar(nn+1)%tcla,ar(nn+1)%rindex,ar(nn+1)%b, &
     &           ar(nn+1)%dd2_ds,                                       &
     &           ar(nn+1)%stretch,ar(nn+1)%width,ar(nn+1)%moid,       &
     &           cont_dist,cont_prob,niter,type,method,iunnew0)      
! copy output close approach array                                      
             va_tracemin=ar(nn+1) 
             fals_conv=.true. 
             RETURN 
          ELSE 
! see if possible to switch to falsi                                    
             IF(xx(nn).lt.xx(nn+1))THEN 
                IF(ar(nn)%dd2_ds.lt.0.d0.and.ar(nn+1)%dd2_ds.gt.0.d0)THEN 
                   WRITE(iunwar,*)' achillestp: minimum identified' 
! copy arrays in the order nn, nn+1]
                   va2_trace(1)=ar(nn)
                   va2_trace(2)=ar(nn+1)                                     
                   CALL falslog4tp(va2_trace,iunwar,type) 
! new type                                                              
                   IF(type.eq.'SIMPMIN')THEN 
                      type='ASIMPMI' 
                   ELSEIF(type.eq.'INTEMIN')THEN 
                      type='AINTEMI' 
! use falsi only for interrupted                                        
                      WRITE(iunwar,*)' achillestp: switch to falsi' 
                      CALL findminctp(va2_trace,xx(nn),xx(nn+1),type,        &
     &                       iunwar,iunnew,va_tracemin,                      &
     &                       fals_conv,niter,fals_notp)                 
                      RETURN 
                   ENDIF
                ENDIF
             ELSE 
                IF(ar(nn+1)%dd2_ds.lt.0.d0.and.ar(nn)%dd2_ds.gt.0.d0)THEN 
                   WRITE(iunwar,*)' achillestp: minimum identified' 
! copy arrays in the order nn+1,nn                                      
                   va2_trace(1)=ar(nn+1)
                   va2_trace(2)=ar(nn)
                   CALL falslog4tp(va2_trace,iunwar,type) 
! new type                                                              
                   IF(type.eq.'SIMPMIN')THEN 
                      type='ASIMPMI' 
                   ELSEIF(type.eq.'INTEMIN')THEN 
                      type='AINTEMI' 
! use falsi only for interrupted                                        
                      WRITE(iunwar,*)' achillestp: switch to falsi' 
                      CALL findminctp(va2_trace,xx(nn+1),xx(nn),type,        &
              &                       iunwar,iunnew,va_tracemin,             &
              &                       fals_conv,niter,fals_notp)                 
                      RETURN 
                   ENDIF
                ENDIF
             ENDIF
          ENDIF
          xlim(nn+1)=xlim(nn)/2.d0 
       ENDIF
1   ENDDO
! case not convergent                                                   
    IF(fals_notp)THEN 
       method='ANTP' 
    ELSE 
       method='AcNC' 
    ENDIF
    niter=nn 
    CALL  wriwarn(ar(nn+1)%tcla,xx(nn+1),ar(nn+1)%b,ar(nn+1)%dd2_ds,         &
     &     ar(nn+1)%stretch,ar(nn+1)%width,ar(nn+1)%moid,                 &
     &     cont_dist,cont_prob,niter,type,method,iunnew0)                
    CALL  wriwarn(ar(nn+1)%tcla,xx(nn+1),ar(nn+1)%b,ar(nn+1)%dd2_ds,         &
     &     ar(nn+1)%stretch,ar(nn+1)%width,ar(nn+1)%moid,                 &
     &     cont_dist,cont_prob,niter,type,method,iunwar)                
! this case not convergent, but the output can be used                  
    fals_conv=.false.
    va_tracemin=ar(nn+1) 
    RETURN 
  END SUBROUTINE achillestp
!===============================================================
! ==============================================================        
! FINDMINCTP  Regula falsi, to be used only when the                      
!           two previous points have opposite sign                      
! ============================================================          
!                                                                       
  SUBROUTINE findminctp(arrc,a0,b0,type,iunwar0,iunnew0,arrmin,         &
     &            fals_conv,niter,fals_notp) 
    USE planet_masses
!========================== INPUT===============================
    DOUBLE PRECISION, INTENT(IN)             :: a0,b0 
    TYPE(tp_point), DIMENSION(2), INTENT(IN) :: arrc   ! the last 
                                               !two used in the return  
    CHARACTER(LEN=7), INTENT(IN)             :: type 
    INTEGER, INTENT(IN)                      :: iunwar0,iunnew0 
!======================= INPUT/OUTPUT===========================
!    TYPE(tp_point), DIMENSION(2) :: arrf       ! array contains the last two 
                                                 ! used in this routine   
!========================== OUTPUT==============================
    TYPE(tp_point), INTENT(OUT)  :: arrmin
    INTEGER, INTENT(OUT)         :: niter 
    LOGICAL, INTENT(OUT)         :: fals_conv,fals_notp 
!======================= END INTERFACE==========================
    INTEGER :: it,nit 
!    INTEGER, PARAMETER :: itmax=200,itma1=itmax+1,itma2=itmax+2 
    INTEGER, PARAMETER :: itmax=30,itma1=itmax+1,itma2=itmax+2 
    DOUBLE PRECISION, DIMENSION(itma1) :: a,b,ta,tb
    DOUBLE PRECISION, DIMENSION(itma1) :: da,db,fa,fb 
    DOUBLE PRECISION, DIMENSION(itma1) :: fw,tw,dw
    DOUBLE PRECISION, DIMENSION(itma2) :: w
    DOUBLE PRECISION, DIMENSION(itma1) :: wstr,wwid,wmoidgf 
    DOUBLE PRECISION :: ffa,ffb 
    CHARACTER(LEN=7) :: type2 
    DOUBLE PRECISION :: ddx,dold,cont_dist,cont_prob,fm  &
    &     ,ddw,cont_dist1,cont_dist2       ! convergence control       
    DOUBLE PRECISION, PARAMETER :: eps_fal=1.d-10 ! limit probability 
                                      ! (with uniform prob.density)
    DOUBLE PRECISION, PARAMETER :: del_fal=2.3d-3 ! limit change in 
                                          ! distance in RE (=1e-7 AU) 
    DOUBLE PRECISION :: tafter 
    LOGICAl          :: falsok
    CHARACTER(LEN=4) :: method 
    INTEGER          :: i, npset, iunwar, iunnew
    DOUBLE PRECISION  :: b_e,vsize,v,mu
!=========================================================================
    iunwar=abs(iunwar0)
    iunnew=abs(iunnew0) 
    WRITE(iunwar,*)'findminctp ',type,arrc(1)%tcla,a0,b0 
    IF(a0.eq.b0)THEN
       WRITE(iunwar,*) ' supect multiple minimum '
       fals_conv=.false.
       niter=0       
       fals_notp=.false.
       arrmin=arrc(1)
       RETURN
    ENDIF
!    arrf(1)=arrc(1)
!    arrf(2)=arrc(2)
! consult arrays for initial two points                                 
    a(1)=a0 
    b(1)=b0 
    ta(1)=arrc(1)%tcla 
    tb(1)=arrc(2)%tcla 
    da(1)=arrc(1)%b 
    db(1)=arrc(2)%b 
    fa(1)=arrc(1)%dd2_ds 
    fb(1)=arrc(2)%dd2_ds 
! store a as w(1)                                                       
      w(1)=a0 
      fw(1)=fa(1) 
      dw(1)=da(1) 
      tw(1)=ta(1) 
      wstr(1)=arrc(1)%stretch 
      wwid(1)=arrc(1)%width 
      wmoidgf(1)=arrc(1)%moid 
! beginning record                                                      
      WRITE(iunwar,120)ta(1),a(1),da(1),fa(1),tb(1),b(1),db(1),fb(1) 
  120 FORMAT('beginning falsi:'/'t1,x1,d1,f1 ',f10.2,1x,f11.5,1p,2d11.3,&
     &  0p/'t2,x2,d2,f2 ',f10.2,1x,f11.5,1p,2d11.3)                     
! control on sign                                                       
      IF(fa(1)*fb(1).gt.0.d0)THEN 
         WRITE(iunwar,*)' this case not suitable for regula falsi' 
      ENDIF 
! first point to be tested                                              
      w(2)=(fb(1)*a(1)-fa(1)*b(1))/(fb(1)-fa(1)) 
! iteartion loop                                                        
      DO 1 it=1,itmax 
         npset=0
! find and propagate to close approach the imterpolated orbit
 11      CALL lovclosapptp(w(it+1),ta(it),tb(it),iunwar,          &
     &        fals_notp,arrmin)                                      
! target plane found?                                                   
         IF(fals_notp)THEN 
            fals_conv=.false. 
            IF(it.eq.1)THEN 
               cont_dist=-1.d0 
               cont_prob=-1.d0 
            ENDIF 
            niter=it 
            write(iunwar,*)' falsi failed: no target plane', w(it+1) 
            method='FNTP' 
            CALL wriwarn(tw(it),w(it),dw(it),fw(it),wstr(it),wwid(it),  &
     &        wmoidgf(it),cont_dist,cont_prob,niter,type,method,iunwar) 
            CALL wriwarn(tw(it),w(it),dw(it),fw(it),wstr(it),wwid(it),  &
     &        wmoidgf(it),cont_dist,cont_prob,niter,type,method,iunnew) 
            RETURN 
         ELSEIF(.not.arrmin%tp_conv)THEN
! set a larger value for npoint
            npset=npset+1
            IF(npset.gt.0)THEN
               niter=it 
               write(iunwar,*)' falsi failed: singularity', w(it+1) 
               write(iunnew,*)' falsi failed: singularity', w(it+1) 
               method='FSIN' 
               tw(it+1)=arrmin%tcla
               dw(it+1)=arrmin%b
               wstr(it+1)=arrmin%stretch 
               wwid(it+1)=arrmin%width
               fw(it+1)=arrmin%dd2_ds
               wmoidgf(it+1)=arrmin%moid
               CALL wriwarn(tw(it+1),w(it+1),dw(it+1),fw(it+1),wstr(it+1),wwid(it+1),  &
     &             wmoidgf(it+1),cont_dist,cont_prob,niter,type,method,iunwar) 
               CALL wriwarn(tw(it+1),w(it+1),dw(it+1),fw(it+1),wstr(it+1),wwid(it+1),  &
     &                wmoidgf(it+1),cont_dist,cont_prob,niter,type,method,iunnew)
               arrmin%rindex=w(it+1)
               RETURN
            ELSE
               write(iunwar,*)' findminctp: npoint increased by ',2**npset
               write(iunnew,*)' findminctp: npoint increased by ',2**npset
 !             CALL npoint_set(2**npset) ! try increasing npoint by a factor 2
               GOTO 11
            ENDIF
         ENDIF
! data on the new point                                                 
         tw(it+1)=arrmin%tcla
         dw(it+1)=arrmin%b
         wstr(it+1)=arrmin%stretch 
         wwid(it+1)=arrmin%width
         fw(it+1)=arrmin%dd2_ds
         wmoidgf(it+1)=arrmin%moid 
! falsi logic                                                           
         IF(fa(it)*fw(it+1).le.0.d0.and.fb(it)*fw(it+1).ge.0.d0)THEN 
! standard regula falsi                                                 
            b(it+1)=w(it+1) 
            fb(it+1)=fw(it+1) 
            db(it+1)=dw(it+1) 
            tb(it+1)=tw(it+1) 
!...........................
!            CALL arrloadtp(arrf(2),b(it+1)) 
            a(it+1)=a(it) 
            fa(it+1)=fa(it) 
            da(it+1)=da(it) 
            ta(it+1)=ta(it) 
! improved regula falsi                                                 
            IF(fw(it)*fw(it+1).gt.0.d0.and.it.gt.1)THEN 
               IF(it.gt.1.and.fw(it-1)*fw(it).gt.0.d0)THEN 
! superaccelerator                                                      
                  ffa=fa(it)/4.d0 
               ELSE 
! accelerator                                                           
                  ffa=fa(it)/2.d0 
               ENDIF 
            ELSE 
! ordinary regula falsi                                                 
               ffa=fa(it) 
            ENDIF 
            ffb=fb(it) 
         ELSEIF(fb(it)*fw(it+1).le.0.d0.and.fa(it)*fw(it+1).ge.0.d0)THEN 
! standard regula falsi                                                 
            a(it+1)=w(it+1) 
            fa(it+1)=fw(it+1) 
            da(it+1)=dw(it+1) 
            ta(it+1)=tw(it+1) 
!            CALL arrloadtp(arrf(1),a(it+1)) 
            b(it+1)=b(it) 
            fb(it+1)=fb(it) 
            db(it+1)=db(it) 
            tb(it+1)=tb(it) 
! improved regula falsi                                                 
            IF(fw(it)*fw(it+1).gt.0.d0.and.it.gt.1)THEN 
               IF(it.gt.1.and.fw(it-1)*fw(it).gt.0.d0)THEN 
! superaccelerator                                                      
                  ffb=fb(it)/4.d0 
               ELSE 
! accelerator                                                           
                  ffb=fb(it)/2.d0 
               ENDIF 
            ELSE 
! ordinary regula falsi                                                 
               ffb=fb(it) 
            ENDIF 
            ffa=fa(it) 
         ELSE 
! these are not regula falsi, but the secant method                     
            WRITE(iunwar,*)' findminctp should not be used with same sign' 
            WRITE(*,*) it, a(it),fa(it),b(it),fb(it),w(it+1),fw(it+1) 
         ENDIF 
! compute next step                                                     
         w(it+2)=(ffb*a(it+1)-ffa*b(it+1))/(ffb-ffa) 
! convergence criteria                                                  
! intervals for control                                                 
         fm=abs(fw(it+1)) 
         ddx=abs(b(it+1)-a(it+1)) 
         cont_dist1=fm*deltasig*ddx/(da(it+1)+db(it+1)) 
         cont_dist2=abs(dw(it+1)-dw(it)) 
         IF(it.eq.1)THEN 
            cont_dist=cont_dist1 
         ELSE 
            cont_dist=min(cont_dist1,cont_dist2) 
         ENDIF 
         cont_prob=ddx*deltasig/6 
         niter=it 
         method='Fa' 
         CALL wriwarn(tw(it+1),w(it+1),dw(it+1),fw(it+1),wstr(it+1),    &
     &        wwid(it+1),wmoidgf(it+1),cont_dist,cont_prob,niter,       &
     &        type,method,iunwar)                                       
         ddw=dw(it+1) 
! check type           
         CALL falslog4tp(arrc,iunwar,type2) 
         IF(type2.ne.type)THEN 
            IF(type.eq.'AINTEMI'.and.type2.eq.'INTEMIN')THEN 
            ELSEIF(type.eq.'ASIMPMI'.and.type2.eq.'SIMPMIN')THEN 
            ELSE 
               WRITE(iunwar,*)' type change, from ',type,' to ',type2 
            ENDIF 
         ENDIF 
! convergence control
!         tp_flag_cur=arrmin%tp_conv
!         IF(tp_flag_cur)THEN 
! on TP, gravitational focusing 
         v=arrmin%opik%coord(1) 
         mu=gmearth/reau**3
         b_e=sqrt(1.d0+(2.d0*mu)/v**2)
         IF(cont_dist.lt.del_fal  & ! .or.cont_dist.lt.del_fal*ddw      &
     &           .or.cont_prob.lt.eps_fal.or.ddw.lt.b_e)THEN           
            WRITE(iunwar,*)' minimum found, iter= ',niter 
            WRITE(*,*)' minimum found, iter= ',niter 
            fals_conv=.true. 
! write record in newton file                                           
            CALL wriwarn(tw(it+1),w(it+1),dw(it+1),fw(it+1),wstr(it+1), &
     &           wwid(it+1),wmoidgf(it+1),cont_dist,cont_prob,niter,    &
     &           type,method,iunnew)                                    
            arrmin%rindex=w(it+1) 
            RETURN 
         ENDIF 
1     ENDDO
! non convergent case                                                   
      niter=itmax 
      WRITE(iunwar,*)' too many iter, ',cont_dist,cont_prob,itmax 
      IF(dw(it).lt.1.d0)THEN 
         fals_conv=.true. 
      ELSE 
         fals_conv=.false. 
      ENDIF 
! write record in newton file                                           
      nit=niter+1 
      method='FaNC' 
      CALL wriwarn(tw(nit),w(nit),dw(nit),fw(nit),wstr(nit),wwid(nit),  &
     &   wmoidgf(nit),cont_dist,cont_prob,niter,type,method,iunwar)     
      CALL wriwarn(tw(nit),w(nit),dw(nit),fw(nit),wstr(nit),wwid(nit),  &
     &   wmoidgf(nit),cont_dist,cont_prob,nit,type,method,iunnew)     
      arrmin%rindex=w(it+1) 
  END SUBROUTINE findminctp
! ==========================================================================
! =========================================================================
SUBROUTINE lovclosapptp(rindex,t1,t2,iunwar,fals_notp,va_tracemin) 
  USE close_app, ONLY: kill_propag
  USE orbit_elements
  USE propag_state
  USE multiple_sol, ONLY: lovinterp
! ========================= INPUT========================================
  DOUBLE PRECISION, INTENT(IN) :: rindex 
  DOUBLE PRECISION, INTENT(IN) ::  t1,t2
  INTEGER, INTENT(IN) :: iunwar 
! ==========================OUTPUT==========================================
  LOGICAL, INTENT(OUT) :: fals_notp 
  TYPE(tp_point), INTENT(OUT) :: va_tracemin        
! =======================END INTERFACE======================================
  LOGICAL :: falsok 
  DOUBLE PRECISION :: tafter,tbefore           ! epoch time
  TYPE(orbit_elem) :: el0,el
  TYPE(orb_uncert) :: unc0,unc
  INTEGER :: ipla 
  DOUBLE PRECISION :: v_infty, v_inf0             ! velocity at the boundary 
                                                      ! of Earth's sphere 
                                                      ! of influence
  LOGICAL :: batch
! =========================================================================== 

  fals_notp=.false. 
! find orbit with interpolated value of sigma                           
  CALL lovinterp(rindex,deltasig,el0,unc0,falsok) 
  IF(.not.falsok) THEN 
     WRITE(iunwar,*)' lovinterp fails', rindex 
     fals_notp=.true. 
     RETURN 
  ENDIF
! select propagation time                                               
  ipla=3 
  v_inf0=v_infty(el0) 
  CALL aftclo2v(ipla,el0%t,t1,t2,v_inf0,tbefore,tafter) 
! availability of covariance for TP analysis online                     
  CALL cov_avai(unc0,el0%coo,el0%coord) 
! reset close approach storage                                          
  njc=0 
  iplam=0 
! propagate                                                             
  batch=.true. 
  CALL pro_ele(el0,tafter,el,unc0,unc) 
!  WRITE(*,*)'lovclosapptp: t1,t2,tbefore,tafter,tclas: ',t1,t2,tbefore,tafter,tcla(1:njc)
  IF(njc.eq.0)THEN 
     WRITE(iunwar,*)' lovclosapp failed; no close approach' 
     falsok=.false. 
  ELSEIF(iplam.ne.3)THEN 
     WRITE(iunwar,*)' lovclosapp failed; close app. planet ', iplam
     falsok=.false. 
     kill_propag=.false.
  ELSEIF(tp_store(njc)%tcla.gt.max(tbefore,tafter).or.tp_store(njc)%tcla.lt.min(tbefore,tafter))THEN   
     IF(kill_propag)THEN
        WRITE(iunwar,*)' lovclosapp failed: previous collision at ',tcla(njc)
     ELSE
        WRITE(iunwar,*)' lovclosapp failed; tcla too far ',tcla(njc) 
     ENDIF
     falsok=.false.
     kill_propag=.false.
  ELSE 
     IF(kill_propag)THEN
        WRITE(iunwar,*)' lovclosapp stopped, collision found '
     ENDIF
     kill_propag=.false.
     falsok=.true.
  ENDIF
  IF(falsok)THEN 
     CALL arrloadtp(va_tracemin,rindex)
     IF(.not.va_tracemin%tp_conv)THEN
        WRITE(iunwar,*)' lovclosapp failed: close approach elliptic ',va_tracemin%tcla
        fals_notp=.true.
     ENDIF  
  ELSE 
     fals_notp=.true. 
  ENDIF
CONTAINS
! ================================================================      
! SUBROUTINE aftclo2v                                                  
! select time interval to get after the close approach                  
! taking into account the relative velocity w.r. to Earth
! ================================================================        
  SUBROUTINE aftclo2v(iplam,t0,t1,t2,v_inf0,tbefore,tafter) 
    USE planet_masses
! ======================INPUT===================================== 
    INTEGER, INTENT(IN) :: iplam             ! planet number
    DOUBLE PRECISION, INTENT(IN) :: t0       ! time of initial 
                                             ! conditions
    DOUBLE PRECISION, INTENT(IN) :: t1,t2    ! time of two known 
                                             ! close approaches
    DOUBLE PRECISION, INTENT(IN) :: v_inf0   ! velocity 
! =====================OUTPUT=====================================
    DOUBLE PRECISION, INTENT(OUT) :: tafter    ! time "after"
    DOUBLE PRECISION, INTENT(OUT) :: tbefore   ! time "before"
                                !  depending upon sense of propagation  
! ===================END INTERFACE=================================      
    DOUBLE PRECISION :: delt_tp       ! time interval to 
                                      ! exit from TP disk    
! ===================================================================   
! warning: really done only for Earth                                   
    IF(ordnam(iplam).ne.'EARTH') THEN 
       WRITE(*,*)' aftclo2v: not to be used for planet ',ordnam(iplam) 
    ENDIF
! time interval to exit from TP disk   
    IF(v_inf0.gt.0.d0)THEN                                 
       delt_tp=4*dmin(iplam)/v_inf0
    ELSE
       delt_tp=365.25d0
    ENDIF
! forced to avoid infinte intervals for v-inf=0
    delt_tp=MIN(delt_tp,365.25d0)
! forced to avoid short intervals for fast encounters                   
    IF(delt_tp.lt.50.d0)delt_tp=50.d0 
! time interval to be clear out of the TP disk is given as deltat       
    IF(t0.lt.t1)THEN 
! future close approaches                                               
       tafter=max(t1,t2)+delt_tp 
       tbefore=min(t1,t2)-delt_tp 
    ELSE 
! past close approaches                                                 
       tafter=min(t1,t2)-delt_tp 
       tbefore=max(t1,t2)+delt_tp 
    ENDIF
  END SUBROUTINE aftclo2v
! =========================================================================
                                                                        


END SUBROUTINE lovclosapptp



! ===========================================================
! SMALL I/O ROUTINES
! ===========================================================

SUBROUTINE wrireptp(vatr,imul1,imul2,iunrep)
! =================================================================
  INTEGER, INTENT(IN)  :: iunrep ! output unit
  DOUBLE PRECISION, INTENT(IN) :: imul1,imul2 ! return interval 
  TYPE(tp_point), INTENT(IN) :: vatr ! tp point with minimum distance
  CHARACTER(LEN=14) :: calend                ! calendar date
  INTEGER iun
! =================================================================
  CALL calendwri(vatr%tcla,calend) 
! write .rep file
  IF(iunrep.ge.0)THEN
     WRITE(iunrep,100)calend,vatr%tcla,vatr%rindex,imul1,imul2,vatr%b,  &
     &     vatr%opik%coord(4),vatr%opik%coord(5),vatr%stretch,  &
     &     vatr%width,vatr%alpha*degrad,vatr%minposs,      &
     &     vatr%stretch_lov,vatr%alpha_lov*degrad,vatr%opik%coord(1)
  ELSEIF(iunrep.lt.0)THEN
     iun=abs(iunrep)
     WRITE(iun,100)calend,vatr%tcla,vatr%rindex,imul1,imul2,vatr%b,  &
     &     vatr%opik%coord(4),vatr%opik%coord(5),vatr%stretch,  &
     &     vatr%width,vatr%alpha*degrad,vatr%minposs,      &
     &     vatr%stretch_lov,vatr%alpha_lov*degrad,vatr%opik%coord(1)
     WRITE(*,100)calend,vatr%tcla,vatr%rindex,imul1,imul2,vatr%b,  &
     &     vatr%opik%coord(4),vatr%opik%coord(5),vatr%stretch,  &
     &     vatr%width,vatr%alpha*degrad,vatr%minposs,      &
     &     vatr%stretch_lov,vatr%alpha_lov*degrad,vatr%opik%coord(1)
  ENDIF
100 FORMAT(a14,1x,f12.5,1x,f6.1,1x,f6.1,1x,f6.1,1x,3(1x,f12.6),1p,          &
     &     2(1x,e10.3),                                                 &
     &     0p,1x,f8.3,1x,f12.6,1x,1p,e11.4,0p,1x,f10.5,                 &
     &     1p,e12.4)               
END SUBROUTINE wrireptp

SUBROUTINE header_rep(iunrep)
  INTEGER, INTENT(IN) :: iunrep ! output unit
  WRITE(iunrep,199) 
199 FORMAT('  CALENDAR DATE ',1x,'  MJD DATE  ',1x,'IMIN'             &
     &     ,1x,'IMU1',1x,'IMU2',1x,' MIN_FOUND ',1x                     &
     &           ,1x,'  CSI (RE) ',1x,' ZETA (RE) ',1x,'STRETCHING', 1x,&
     &  '   WIDTH  ',1x,' ANG_ELL',1x,'  MOID     ',                    &
     &        1x,'STRETC_LOV',1x,' ANG_LOV')     
END SUBROUTINE header_rep

SUBROUTINE wriouttp(vatr,iunout) 
! =========================================================================
  INTEGER :: iunout
  TYPE(tp_point), INTENT(IN) :: vatr  
! ====================================================================
! write .out file                                                       
  WRITE(iunout,100)vatr%rindex,vatr%tcla,vatr%opik%coord(4), &
     &     vatr%opik%coord(5),                 &
     &     vatr%stretch,vatr%width,vatr%alpha*degrad,       &
     &     vatr%minposs,                  &
     &     vatr%stretch_lov,vatr%alpha_lov*degrad,vatr%opik%coord(1)
100 FORMAT(f6.1,1x,f12.5,2(1x,f12.6),1p,2(1x,e10.3),                    &
     &     0p,1x,f8.3,1x,f12.6,1x,1p,e11.4,0p,1x,f8.3,                 &
     &     1p,e12.4)                         
END SUBROUTINE wriouttp
         
!========================================================================= 
SUBROUTINE wriwarn(tw,w,dw,fw,wstr,wwid,wmoidgf,                  &
     &     cont_dist,cont_prob,niter,type,method,iunwar)                
! ===========================INPUT=========================================
  DOUBLE PRECISION, INTENT(IN) :: tw,w,dw,fw,wstr,wwid,wmoidgf ! data 
  ! on current close approach 
  DOUBLE PRECISION, INTENT(IN) ::  cont_dist,cont_prob  ! convergence control 
  INTEGER, INTENT(IN) :: niter               ! number of iterations so far
  CHARACTER(LEN=7), INTENT(IN) :: type       ! string comments 
  CHARACTER(LEN=4), INTENT(IN) :: method 
  INTEGER,INTENT(IN) :: iunwar               ! output unit
! ==========================================================================
  CHARACTER(LEN=14) :: calend                ! calendar date 
  INTEGER iun
! =========================================================================
  CALL calendwri(tw,calend) 
  IF(iunwar.gt.0)THEN
     WRITE(iunwar,119)calend,tw,w,dw,fw,wstr                           &
     &     ,wwid,wmoidgf,cont_dist,cont_prob,niter,type,method
  ELSEIF(iunwar.eq.0)THEN
     CALL header_new(iunwar)
     WRITE(*,119)calend,tw,w,dw,fw,wstr                           &
     &     ,wwid,wmoidgf,cont_dist,cont_prob,niter,type,method
  ELSEIF(iunwar.lt.0)THEN
     iun=abs(iunwar)
     CALL header_new(0)
     WRITE(*,119)calend,tw,w,dw,fw,wstr                           &
     &     ,wwid,wmoidgf,cont_dist,cont_prob,niter,type,method
     WRITE(iun,119)calend,tw,w,dw,fw,wstr                           &
     &     ,wwid,wmoidgf,cont_dist,cont_prob,niter,type,method
  ENDIF
119 FORMAT(a14,1x,f12.5,1x,f10.5,1x,f7.2,1x,1p,d9.2,1x,d9.2,1x,d9.2,  &
     &           1x,d9.2,1x,d9.2,1x,d9.2,1x,i3,1x,a7,1x,a4)             
END SUBROUTINE wriwarn

SUBROUTINE header_new(iunnew)
INTEGER, INTENT(IN) :: iunnew ! output unit
  WRITE(iunnew,198) 
198 FORMAT('  CALENDAR DATE ',1x,'      MJD  ',1x,'  INDEX  ',1x,     &
     &       'LOV_DIST ',1x,'DR2/DSIG',1x,'  STRETCH.', 1x,'  WIDTH ',  &
     &       '   MOIDGF',1x,'  DIST.CON',1x,' PROB.CON',1x,' IT',1x,    &
     & 'TYPE',1x,'METH')  
END SUBROUTINE header_new

!===========================================================================
! SUBROUTINE arrcut 
!===========================================================================
SUBROUTINE arrcut(vas_trace,ire,lre,nox,iunout,vas_traceloc) 
  INTEGER, INTENT(IN) :: iunout   ! unit for out record
  INTEGER, INTENT(IN) :: lre,ire,nox  ! address and length of filament, dimension
  TYPE(tp_point), DIMENSION(nox), INTENT(IN) :: vas_trace 
                                      ! global close approach records,   
! -------------------LOCAL ARRAYS (for a single trail)--------------------
  TYPE(tp_point), DIMENSION(lre), INTENT(OUT) :: vas_traceloc     
  INTEGER j,k 
!-------------------------------------------------------------------------
! copy and traslate
  DO j=1,lre
     vas_traceloc(j)=vas_trace(j+ire-1) 
     CALL wriouttp(vas_traceloc(j),iunout) 
  ENDDO
END SUBROUTINE arrcut
! ===========================================================================
END MODULE ret_analysistp





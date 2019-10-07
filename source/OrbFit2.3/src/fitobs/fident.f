c ========================================
c FIDENT compute identification norm
c =======================================
      SUBROUTINE fident(igue,t0,tp,cov0,covp,eq0,eqp,g0,c0,gp,cp,
     +    to0,top,iun20,eqid,tid,ff)
      IMPLICIT NONE
      LOGICAL cov0,covp,ff,fail
      INTEGER iun20,igue
      DOUBLE PRECISION eq0(6),eqp(6),t0,tp
      DOUBLE PRECISION c0(6,6),cp(6,6),g0(6,6),gp(6,6)
c similarity norms 
      DOUBLE PRECISION to0,top,eqid(6),tid,eqf(6)
      DOUBLE PRECISION d2,da2,dista,dq,dqalt
c determinants and eigenvalues
      DOUBLE PRECISION detc2,detc5,detc6,eigen5(5),eigen6(6)
c ========================================
c check requirements
      IF(.not.cov0)THEN
         WRITE(*,*)' covariance matrix for arc 1 not ready'
         ff=.false.
         RETURN
      ENDIF
      IF(.not.covp)THEN
         WRITE(*,*)' covariance matrix for arc 2 not ready'
         ff=.false.
         RETURN
      ENDIF
      IF(t0.ne.tp)THEN
         WRITE(*,*)' to=',t0,' tp=',tp
         WRITE(*,*)' elements and covariances must be available'
         WRITE(*,*)' for the same time, use propagation first'
         RETURN
      ENDIF
c selection of algorithm
      IF(igue.eq.5)THEN
c identification by 5x5 matrix
         CALL idno5(eq0,eqp,g0,c0,gp,cp,
     +    d2,da2,dq,dqalt,dista,detc2,detc5,eqf,fail,eigen5)
      ELSEIF(igue.eq.6)THEN
c identification by 6x6 matrix
         CALL idno6(eq0,eqp,g0,c0,gp,cp,
     +    d2,da2,dq,dqalt,dista,detc2,detc6,eqf,fail,eigen6)
      ELSE
         WRITE(*,*)' fident: igue=',igue,' not understood'
         RETURN
      ENDIF
c output
      IF(fail)THEN
         CALL tee(iun20,' FAILED IDENTIFICATION ALGORITHM=')
         ff=.false.
      ELSE
         WRITE(*,*)' INITIAL COND. AVAILABLE, NOW USE DIFF. CORR.'
         ff=.true.
      ENDIF
c penalties
      WRITE(*,*)' d2, d2alt,  dq,   dqalt,    dista'
      WRITE(*,194)d2,da2,dq,dqalt,dista  
 194  FORMAT(4(f13.4,1x),f10.6)
      WRITE(iun20,*)' d2, d2alt,  dq,   dqalt,    dista'
      WRITE(iun20,194)d2,da2,dq,dqalt,dista  
c proposed elements
      WRITE(iun20,208)t0
 208  FORMAT(' ident. elem (a,h,k,p,q,lam), epoch=',f8.1)
      WRITE(*,208)t0
      WRITE(*,104)eqf
      WRITE(iun20,104)eqf 
 104  FORMAT(6f13.7)
c store as proposed identification
      IF(ff)THEN
         CALL vcopy(6,eqf,eqid)
         tid=(t0+tp)/2.d0
      ENDIF
      RETURN
      END







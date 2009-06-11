SUBROUTINE determinants(pe,pde,rho,rhohat,rhohat_a,rhohat_d,rhodot,&
     & detM, detJ1, detK1, detJ2, detK2, detJ3, detK3)
  IMPLICIT NONE
  DOUBLE PRECISION,DIMENSION(3,2),INTENT(IN) :: pe,pde
  DOUBLE PRECISION,DIMENSION(2),INTENT(IN) :: rho
  DOUBLE PRECISION,DIMENSION(3,2),INTENT(IN) :: rhohat,rhohat_a,rhohat_d
  DOUBLE PRECISION,DIMENSION(2),INTENT(IN) :: rhodot
  DOUBLE PRECISION,INTENT(OUT) :: detM,detJ1,detK1,detJ2
  DOUBLE PRECISION,INTENT(OUT) :: detK2,detJ3,detK3
! -----------------------------------------------------
! auxiliary
  DOUBLE PRECISION :: q11,q12,q13,q21,q22,q23
  DOUBLE PRECISION :: qd11,qd12,qd13,qd21,qd22,qd23
  DOUBLE PRECISION :: r11,r12,r13,r21,r22,r23
  DOUBLE PRECISION :: rhoh11,rhoh12,rhoh13,rhoh21,rhoh22,rhoh23
  DOUBLE PRECISION :: rhoha11,rhoha12,rhoha13,rhoha21,rhoha22,rhoha23
  DOUBLE PRECISION :: rhohd11,rhohd12,rhohd13,rhohd21,rhohd22,rhohd23
  DOUBLE PRECISION :: rho1,rho2,rhodot1,rhodot2
  DOUBLE PRECISION :: s1,s2,s3,s4,s5,s6,s7,s8

  q11=pe(1,1)
  q12=pe(2,1)
  q13=pe(3,1)
  q21=pe(1,2)
  q22=pe(2,2)
  q23=pe(3,2)

  qd11=pde(1,1)
  qd12=pde(2,1)
  qd13=pde(3,1)
  qd21=pde(1,2)
  qd22=pde(2,2)
  qd23=pde(3,2)

  rhoh11=rhohat(1,1)
  rhoh12=rhohat(2,1)
  rhoh13=rhohat(3,1)
  rhoh21=rhohat(1,2)
  rhoh22=rhohat(2,2)
  rhoh23=rhohat(3,2)

  rho1=rho(1)
  rho2=rho(2)

  r11= q11 + rho1*rhoh11
  r12= q12 + rho1*rhoh12
  r13= q13 + rho1*rhoh13
  r21= q21 + rho2*rhoh21
  r22= q22 + rho2*rhoh22
  r23= q23 + rho2*rhoh23

  rhoha11=rhohat_a(1,1)
  rhoha12=rhohat_a(2,1)
  rhoha13=rhohat_a(3,1)
  rhoha21=rhohat_a(1,2)
  rhoha22=rhohat_a(2,2)
  rhoha23=rhohat_a(3,2)

  rhohd11=rhohat_d(1,1)
  rhohd12=rhohat_d(2,1)
  rhohd13=rhohat_d(3,1)
  rhohd21=rhohat_d(1,2)
  rhohd22=rhohat_d(2,2)
  rhohd23=rhohat_d(3,2)

  rhodot1=rhodot(1)
  rhodot2=rhodot(2)

  detM = -rho1**2*rho2*(-rhoha13*r11*rhohd12+rhoha11*r13*rhohd12+r11&
     &*rhoha12*rhohd13-r13*rhoha12*rhohd11+r12*rhoha13*rhohd11-r12*rhoha&
     &11*rhohd13)*(r11*r22*rhoha23-r11*r23*rhoha22-r13*r22*rhoha21+r13*r&
     &21*rhoha22+r12*r23*rhoha21-r12*r21*rhoha23)                       
                                                                        
      s1 = -rho1 
      s3 = rho2 
      s8 = -rhodot1*q11*rhoh13*r22*rhoha23*r12*rhohd11+rhodot1*q12*rhoh1&
     &3*r13*rhohd11*r22*rhoha21+rhodot2*q23*rhoh21*r23*rhoha22*r12*rhohd&
     &11-q23*qd22*r11*rhohd13*r22*rhoha21-rhodot1*q12*rhoh13*r13*rhohd11&
     &*r21*rhoha22+rhodot1*q13*rhoh11*r12*rhohd13*r21*rhoha22+q13*qd12*r&
     &23*rhoha21*r12*rhohd11+rhodot1*q12*rhoh11*r22*rhoha23*r11*rhohd13+&
     &rhodot1*q12*rhoh11*r23*rhoha22*r13*rhohd11+rhodot2*q22*rhoh21*r23*&
     &rhoha22*r11*rhohd13-rhodot1*q11*rhoh12*r12*rhohd13*r23*rhoha21-rho&
     &dot2*q22*rhoh21*r23*rhoha22*r13*rhohd11                           
      s7 = s8+rhodot1*q13*rhoh12*r13*rhohd11*r21*rhoha22-q12*qd13*r11*rh&
     &ohd13*r22*rhoha21+rhodot2*q23*rhoh21*r22*rhoha23*r11*rhohd12+rhodo&
     &t1*q12*rhoh13*r23*rhoha21*r11*rhohd12-q11*qd12*r22*rhoha23*r11*rho&
     &hd13+q11*qd12*r22*rhoha23*r13*rhohd11-q13*qd11*r13*rhohd12*r21*rho&
     &ha22-q11*qd12*r12*rhohd13*r23*rhoha21+q22*qd21*r12*rhohd13*r21*rho&
     &ha23-q23*qd21*r23*rhoha22*r11*rhohd12-rhodot2*q23*rhoh21*r13*rhohd&
     &12*r22*rhoha21-rhodot1*q13*rhoh12*r23*rhoha21*r11*rhohd12         
      s8 = rhodot2*q23*rhoh21*r13*rhohd12*r21*rhoha22-q21*qd23*r23*rhoha&
     &22*r12*rhohd11+q22*qd21*r22*rhoha23*r13*rhohd11+rhodot2*q21*rhoh22&
     &*r13*rhohd12*r21*rhoha23+s7+rhodot1*q13*rhoh12*r11*rhohd13*r22*rho&
     &ha21+rhodot1*q12*rhoh11*r12*rhohd13*r23*rhoha21-rhodot1*q12*rhoh11&
     &*r12*rhohd13*r21*rhoha23-rhodot2*q23*rhoh22*r23*rhoha21*r12*rhohd1&
     &1-rhodot2*q21*rhoh23*r22*rhoha23*r11*rhohd12+rhodot2*q21*rhoh22*r2&
     &3*rhoha22*r13*rhohd11+rhodot2*q21*rhoh23*r12*rhohd13*r21*rhoha22  
      s6 = s8-q22*qd21*r22*rhoha23*r11*rhohd13-rhodot2*q23*rhoh21*r12*rh&
     &ohd13*r21*rhoha22+rhodot1*q13*rhoh12*r21*rhoha23*r11*rhohd12-rhodo&
     &t1*q11*rhoh12*r13*rhohd12*r21*rhoha23+q22*qd23*r21*rhoha23*r11*rho&
     &hd12-rhodot1*q13*rhoh11*r13*rhohd12*r21*rhoha22+rhodot2*q22*rhoh21&
     &*r22*rhoha23*r13*rhohd11+rhodot2*q22*rhoh23*r23*rhoha21*r12*rhohd1&
     &1-rhodot1*q13*rhoh12*r21*rhoha23*r12*rhohd11-q13*qd12*r21*rhoha23*&
     &r12*rhohd11+q13*qd11*r13*rhohd12*r22*rhoha21+q12*qd13*r13*rhohd11*&
     &r22*rhoha21-q11*qd13*r22*rhoha23*r12*rhohd11                      
      s8 = s6+rhodot1*q11*rhoh12*r22*rhoha23*r13*rhohd11-rhodot2*q22*rho&
     &h23*r23*rhoha21*r11*rhohd12-rhodot1*q12*rhoh11*r22*rhoha23*r13*rho&
     &hd11-q23*qd22*r23*rhoha21*r12*rhohd11-q23*qd21*r12*rhohd13*r21*rho&
     &ha22+rhodot1*q11*rhoh12*r13*rhohd12*r23*rhoha21+rhodot1*q11*rhoh12&
     &*r12*rhohd13*r21*rhoha23+rhodot1*q13*rhoh12*r23*rhoha21*r12*rhohd1&
     &1+rhodot2*q23*rhoh22*r21*rhoha23*r12*rhohd11-q21*qd22*r12*rhohd13*&
     &r21*rhoha23+q23*qd21*r12*rhohd13*r22*rhoha21                      
      s7 = s8+rhodot2*q21*rhoh23*r22*rhoha23*r12*rhohd11+q12*qd13*r21*rh&
     &oha23*r12*rhohd11+rhodot2*q21*rhoh22*r22*rhoha23*r11*rhohd13-rhodo&
     &t1*q12*rhoh11*r13*rhohd12*r23*rhoha21+q23*qd21*r22*rhoha23*r11*rho&
     &hd12+q22*qd21*r23*rhoha22*r11*rhohd13+q21*qd23*r13*rhohd12*r22*rho&
     &ha21-rhodot2*q22*rhoh21*r13*rhohd12*r21*rhoha23-rhodot2*q22*rhoh21&
     &*r22*rhoha23*r11*rhohd13-q21*qd23*r12*rhohd13*r22*rhoha21-q21*qd22&
     &*r23*rhoha22*r11*rhohd13+rhodot1*q13*rhoh11*r13*rhohd12*r22*rhoha2&
     &1                                                                 
      s8 = s7-q12*qd13*r13*rhohd11*r21*rhoha22+rhodot1*q11*rhoh12*r23*rh&
     &oha22*r11*rhohd13-q23*qd22*r21*rhoha23*r11*rhohd12+q23*qd21*r13*rh&
     &ohd12*r21*rhoha22-q13*qd11*r12*rhohd13*r22*rhoha21+q12*qd11*r13*rh&
     &ohd12*r21*rhoha23-rhodot2*q22*rhoh23*r11*rhohd13*r21*rhoha22-rhodo&
     &t2*q22*rhoh23*r13*rhohd11*r22*rhoha21+q22*qd21*r13*rhohd12*r23*rho&
     &ha21-rhodot1*q12*rhoh13*r11*rhohd13*r22*rhoha21-rhodot2*q23*rhoh21&
     &*r23*rhoha22*r11*rhohd12-rhodot2*q23*rhoh22*r13*rhohd11*r21*rhoha2&
     &2                                                                 
      s5 = s8-rhodot1*q13*rhoh12*r11*rhohd13*r21*rhoha22-rhodot2*q22*rho&
     &h23*r21*rhoha23*r12*rhohd11+rhodot2*q22*rhoh23*r11*rhohd13*r22*rho&
     &ha21+rhodot2*q23*rhoh22*r23*rhoha21*r11*rhohd12-rhodot1*q13*rhoh11&
     &*r22*rhoha23*r11*rhohd12-q22*qd21*r13*rhohd12*r21*rhoha23+rhodot1*&
     &q11*rhoh13*r12*rhohd13*r22*rhoha21+rhodot1*q11*rhoh13*r22*rhoha23*&
     &r11*rhohd12+rhodot2*q22*rhoh23*r21*rhoha23*r11*rhohd12-rhodot2*q21&
     &*rhoh22*r13*rhohd12*r23*rhoha21+q23*qd22*r21*rhoha23*r12*rhohd11-q&
     &21*qd23*r13*rhohd12*r21*rhoha22-q21*qd22*r13*rhohd12*r23*rhoha21  
      s8 = q13*qd11*r12*rhohd13*r21*rhoha22+q21*qd23*r22*rhoha23*r12*rho&
     &hd11+q21*qd23*r23*rhoha22*r11*rhohd12-rhodot1*q11*rhoh12*r23*rhoha&
     &22*r13*rhohd11-rhodot2*q21*rhoh22*r22*rhoha23*r13*rhohd11+rhodot2*&
     &q23*rhoh21*r12*rhohd13*r22*rhoha21-q11*qd13*r13*rhohd12*r22*rhoha2&
     &1-rhodot2*q22*rhoh21*r12*rhohd13*r23*rhoha21-q23*qd21*r22*rhoha23*&
     &r12*rhohd11-q11*qd13*r23*rhoha22*r11*rhohd12+q11*qd13*r13*rhohd12*&
     &r21*rhoha22+q13*qd12*r11*rhohd13*r22*rhoha21                      
      s7 = s8+rhodot2*q23*rhoh22*r13*rhohd11*r22*rhoha21+rhodot2*q21*rho&
     &h22*r12*rhohd13*r23*rhoha21-rhodot1*q11*rhoh13*r13*rhohd12*r22*rho&
     &ha21+rhodot1*q11*rhoh13*r13*rhohd12*r21*rhoha22-rhodot2*q21*rhoh22&
     &*r12*rhohd13*r21*rhoha23+rhodot1*q12*rhoh13*r21*rhoha23*r12*rhohd1&
     &1+rhodot2*q21*rhoh23*r13*rhohd12*r22*rhoha21-rhodot2*q21*rhoh23*r2&
     &3*rhoha22*r12*rhohd11-q12*qd13*r23*rhoha21*r12*rhohd11+rhodot2*q23&
     &*rhoh22*r11*rhohd13*r21*rhoha22-q13*qd12*r11*rhohd13*r21*rhoha22+q&
     &13*qd12*r13*rhohd11*r21*rhoha22                                   
      s8 = s7-q12*qd11*r12*rhohd13*r21*rhoha23+q22*qd23*r11*rhohd13*r22*&
     &rhoha21-q22*qd21*r12*rhohd13*r23*rhoha21+q12*qd13*r11*rhohd13*r21*&
     &rhoha22+q13*qd12*r21*rhoha23*r11*rhohd12+q21*qd22*r12*rhohd13*r23*&
     &rhoha21-rhodot1*q13*rhoh12*r13*rhohd11*r22*rhoha21-rhodot1*q13*rho&
     &h11*r23*rhoha22*r12*rhohd11+rhodot1*q12*rhoh11*r13*rhohd12*r21*rho&
     &ha23+rhodot1*q13*rhoh11*r22*rhoha23*r12*rhohd11+rhodot1*q13*rhoh11&
     &*r23*rhoha22*r11*rhohd12                                          
      s6 = s8+q21*qd22*r13*rhohd12*r21*rhoha23-q11*qd12*r23*rhoha22*r13*&
     &rhohd11+q11*qd13*r22*rhoha23*r11*rhohd12-q23*qd22*r13*rhohd11*r21*&
     &rhoha22-q11*qd12*r13*rhohd12*r21*rhoha23-rhodot2*q23*rhoh22*r11*rh&
     &ohd13*r22*rhoha21+q11*qd13*r23*rhoha22*r12*rhohd11-rhodot2*q21*rho&
     &h23*r12*rhohd13*r22*rhoha21-q22*qd23*r13*rhohd11*r22*rhoha21-rhodo&
     &t2*q23*rhoh22*r21*rhoha23*r11*rhohd12+q12*qd11*r12*rhohd13*r23*rho&
     &ha21-rhodot2*q23*rhoh21*r22*rhoha23*r12*rhohd11-rhodot1*q12*rhoh13&
     &*r23*rhoha21*r12*rhohd11                                          
      s8 = s6-rhodot2*q21*rhoh23*r13*rhohd12*r21*rhoha22-rhodot1*q13*rho&
     &h11*r12*rhohd13*r22*rhoha21+rhodot1*q11*rhoh13*r23*rhoha22*r12*rho&
     &hd11-rhodot1*q11*rhoh13*r23*rhoha22*r11*rhohd12-rhodot1*q11*rhoh13&
     &*r12*rhohd13*r21*rhoha22+rhodot2*q22*rhoh21*r12*rhohd13*r21*rhoha2&
     &3+rhodot2*q22*rhoh21*r13*rhohd12*r23*rhoha21+q13*qd11*r22*rhoha23*&
     &r12*rhohd11+q11*qd12*r13*rhohd12*r23*rhoha21-q21*qd22*r22*rhoha23*&
     &r13*rhohd11+q11*qd12*r23*rhoha22*r11*rhohd13                      
      s7 = s8+q22*qd23*r13*rhohd11*r21*rhoha22-q23*qd21*r13*rhohd12*r22*&
     &rhoha21-rhodot1*q11*rhoh12*r22*rhoha23*r11*rhohd13+q23*qd21*r23*rh&
     &oha22*r12*rhohd11+rhodot1*q12*rhoh13*r11*rhohd13*r21*rhoha22-q13*q&
     &d11*r23*rhoha22*r12*rhohd11-rhodot1*q12*rhoh13*r21*rhoha23*r11*rho&
     &hd12+rhodot2*q21*rhoh23*r23*rhoha22*r11*rhohd12-q13*qd12*r23*rhoha&
     &21*r11*rhohd12-q13*qd12*r13*rhohd11*r22*rhoha21-rhodot2*q21*rhoh22&
     &*r23*rhoha22*r11*rhohd13+s5-q12*qd11*r13*rhohd12*r23*rhoha21      
      s8 = s7-q22*qd21*r23*rhoha22*r13*rhohd11-q13*qd11*r22*rhoha23*r11*&
     &rhohd12-q22*qd23*r21*rhoha23*r12*rhohd11-rhodot1*q12*rhoh11*r23*rh&
     &oha22*r11*rhohd13+rhodot2*q22*rhoh23*r13*rhohd11*r21*rhoha22-q22*q&
     &d23*r23*rhoha21*r11*rhohd12+q21*qd22*r22*rhoha23*r11*rhohd13+q12*q&
     &d11*r23*rhoha22*r13*rhohd11+q13*qd11*r23*rhoha22*r11*rhohd12+q21*q&
     &d23*r12*rhohd13*r21*rhoha22-q21*qd23*r22*rhoha23*r11*rhohd12-q12*q&
     &d13*r21*rhoha23*r11*rhohd12                                       
      s4 = s8+q12*qd13*r23*rhoha21*r11*rhohd12-q22*qd23*r11*rhohd13*r21*&
     &rhoha22+q12*qd11*r22*rhoha23*r11*rhohd13+q23*qd22*r13*rhohd11*r22*&
     &rhoha21+q23*qd22*r11*rhohd13*r21*rhoha22-q12*qd11*r22*rhoha23*r13*&
     &rhohd11+q22*qd23*r23*rhoha21*r12*rhohd11+q21*qd22*r23*rhoha22*r13*&
     &rhohd11+q11*qd13*r12*rhohd13*r22*rhoha21-q12*qd11*r23*rhoha22*r11*&
     &rhohd13+q11*qd12*r12*rhohd13*r21*rhoha23+q23*qd22*r23*rhoha21*r11*&
     &rhohd12-q11*qd13*r12*rhohd13*r21*rhoha22                          
      s2 = s3*s4 
      detJ1 = s1*s2 
                                                                        
      detK1 = -rho1*rho2**2*(-rhohd22*r23*rhoha21-r21*rhohd23*rhoha22+r2&
     &1*rhohd22*rhoha23+r23*rhohd21*rhoha22-r22*rhohd21*rhoha23+rhohd23*&
     &r22*rhoha21)*(-r21*r12*rhohd13+r21*r13*rhohd12+r23*r12*rhohd11-r23&
     &*r11*rhohd12-r22*r13*rhohd11+r22*r11*rhohd13)                     
                                                                        
      s1 = -rho1 
      s3 = rho2 
      s8 = r13*rhoha11*q13*qd12*r22*rhoha21+r13*rhoha12*r21*rhoha23*q22*&
     &qd21+r12*rhoha11*q13*qd11*r23*rhoha22-r11*rhoha13*q23*qd22*r21*rho&
     &ha22+r13*rhoha11*q23*qd22*r21*rhoha22+r11*rhoha13*rhodot2*q21*rhoh&
     &22*r23*rhoha22+r11*rhoha13*rhodot1*q12*rhoh13*r22*rhoha21-r11*rhoh&
     &a12*rhodot1*q13*rhoh11*r23*rhoha22+r11*rhoha12*q13*qd12*r23*rhoha2&
     &1-r11*rhoha12*q22*qd23*r21*rhoha23-r12*rhoha13*rhodot2*q23*rhoh21*&
     &r22*rhoha21+r11*rhoha12*rhodot2*q23*rhoh22*r21*rhoha23            
      s7 = s8-r11*rhoha12*rhodot1*q12*rhoh13*r23*rhoha21-r11*rhoha12*rho&
     &dot1*q13*rhoh12*r21*rhoha23+r12*rhoha11*rhodot1*q12*rhoh13*r23*rho&
     &ha21+r11*rhoha12*rhodot1*q13*rhoh11*r22*rhoha23-r12*rhoha11*rhodot&
     &2*q22*rhoh23*r23*rhoha21+r12*rhoha11*rhodot1*q11*rhoh13*r22*rhoha2&
     &3-r11*rhoha13*rhodot1*q12*rhoh13*r21*rhoha22+r11*rhoha13*q22*qd21*&
     &r22*rhoha23+r11*rhoha13*rhodot2*q22*rhoh23*r21*rhoha22+r13*rhoha12&
     &*rhodot1*q13*rhoh11*r21*rhoha22-r13*rhoha12*r21*rhoha23*rhodot2*q2&
     &1*rhoh22-r13*rhoha12*r23*rhoha21*q22*qd21                         
      s8 = s7+r13*rhoha12*rhodot2*q23*rhoh21*r22*rhoha21+r13*rhoha12*rho&
     &dot2*q21*rhoh23*r21*rhoha22-r12*rhoha11*q23*qd21*r23*rhoha22-r12*r&
     &hoha11*q11*qd13*r23*rhoha22-r12*rhoha11*rhodot2*q23*rhoh21*r23*rho&
     &ha22-r12*rhoha13*r21*rhoha23*rhodot2*q22*rhoh21-r12*rhoha13*r21*rh&
     &oha23*rhodot1*q11*rhoh12+r12*rhoha11*rhodot2*q21*rhoh23*r23*rhoha2&
     &2+r12*rhoha13*rhodot1*q11*rhoh13*r21*rhoha22+r13*rhoha12*r23*rhoha&
     &21*q21*qd22-r13*rhoha11*rhodot1*q11*rhoh12*r22*rhoha23            
      s6 = s8-r11*rhoha12*rhodot1*q11*rhoh13*r22*rhoha23-r11*rhoha13*q21&
     &*qd22*r22*rhoha23-r12*rhoha13*r23*rhoha21*rhodot1*q12*rhoh11-r13*r&
     &hoha11*rhodot1*q12*rhoh13*r22*rhoha21-r13*rhoha11*rhodot1*q13*rhoh&
     &12*r21*rhoha22+r13*rhoha11*rhodot1*q13*rhoh12*r22*rhoha21-r13*rhoh&
     &a11*rhodot2*q21*rhoh22*r23*rhoha22+r13*rhoha12*q23*qd21*r22*rhoha2&
     &1+r12*rhoha13*q13*qd11*r22*rhoha21-r13*rhoha11*q23*qd22*r22*rhoha2&
     &1-r13*rhoha12*q21*qd23*r22*rhoha21-r13*rhoha11*q12*qd13*r22*rhoha2&
     &1+r13*rhoha11*rhodot1*q11*rhoh12*r23*rhoha22                      
      s8 = s6+r12*rhoha13*rhodot2*q23*rhoh21*r21*rhoha22+r13*rhoha11*q12&
     &*qd13*r21*rhoha22+r11*rhoha13*q23*qd22*r22*rhoha21+r13*rhoha12*q13&
     &*qd11*r21*rhoha22-r11*rhoha13*rhodot1*q11*rhoh12*r23*rhoha22+r13*r&
     &hoha12*r21*rhoha23*rhodot1*q11*rhoh12-r13*rhoha12*r21*rhoha23*rhod&
     &ot1*q12*rhoh11-r12*rhoha11*rhodot1*q12*rhoh13*r21*rhoha23-r12*rhoh&
     &a11*rhodot1*q13*rhoh12*r23*rhoha21-r12*rhoha13*q23*qd21*r22*rhoha2&
     &1+r12*rhoha13*r23*rhoha21*rhodot2*q22*rhoh21                      
      s7 = s8-r11*rhoha13*rhodot2*q22*rhoh23*r22*rhoha21+r13*rhoha12*r21&
     &*rhoha23*rhodot2*q22*rhoh21-r11*rhoha13*rhodot1*q12*rhoh11*r22*rho&
     &ha23-r12*rhoha13*r23*rhoha21*q12*qd11+r13*rhoha11*rhodot1*q12*rhoh&
     &13*r21*rhoha22+r12*rhoha13*q11*qd13*r21*rhoha22+r13*rhoha12*r21*rh&
     &oha23*q11*qd12-r13*rhoha12*r21*rhoha23*q12*qd11-r13*rhoha11*q13*qd&
     &12*r21*rhoha22-r12*rhoha13*r23*rhoha21*rhodot2*q21*rhoh22-r13*rhoh&
     &a12*rhodot2*q21*rhoh23*r22*rhoha21+r13*rhoha12*r23*rhoha21*q12*qd1&
     &1                                                                 
      s8 = s7+r11*rhoha13*rhodot1*q13*rhoh12*r21*rhoha22+r13*rhoha11*q11&
     &*qd12*r23*rhoha22+r12*rhoha11*rhodot2*q23*rhoh21*r22*rhoha23+r11*r&
     &hoha12*rhodot2*q21*rhoh23*r22*rhoha23-r13*rhoha11*q21*qd22*r23*rho&
     &ha22+r11*rhoha12*q23*qd21*r23*rhoha22+r13*rhoha11*q22*qd23*r22*rho&
     &ha21-r11*rhoha13*q13*qd12*r22*rhoha21-r11*rhoha12*rhodot2*q23*rhoh&
     &22*r23*rhoha21+r12*rhoha13*rhodot2*q21*rhoh23*r22*rhoha21+r13*rhoh&
     &a12*r23*rhoha21*rhodot2*q21*rhoh22+r12*rhoha13*q23*qd21*r21*rhoha2&
     &2                                                                 
      s5 = s8-r13*rhoha11*q22*qd23*r21*rhoha22+r12*rhoha13*r21*rhoha23*q&
     &21*qd22-r12*rhoha13*r21*rhoha23*q11*qd12-r11*rhoha13*rhodot2*q22*r&
     &hoh21*r23*rhoha22-r11*rhoha12*q12*qd13*r23*rhoha21+r11*rhoha13*q21&
     &*qd22*r23*rhoha22-r11*rhoha13*q11*qd12*r23*rhoha22-r11*rhoha13*q22&
     &*qd21*r23*rhoha22-r11*rhoha12*rhodot2*q23*rhoh21*r22*rhoha23+r12*r&
     &hoha11*q23*qd22*r23*rhoha21-r12*rhoha11*rhodot1*q11*rhoh13*r23*rho&
     &ha22-r12*rhoha11*q21*qd23*r22*rhoha23-r12*rhoha11*q13*qd11*r22*rho&
     &ha23                                                              
      s8 = r11*rhoha12*rhodot1*q11*rhoh13*r23*rhoha22-r11*rhoha13*q12*qd&
     &13*r21*rhoha22+r12*rhoha11*rhodot1*q13*rhoh12*r21*rhoha23+r12*rhoh&
     &a11*q11*qd13*r22*rhoha23-r13*rhoha12*q23*qd21*r21*rhoha22-r11*rhoh&
     &a12*q21*qd23*r23*rhoha22+r11*rhoha13*q13*qd12*r21*rhoha22+r11*rhoh&
     &a12*q21*qd23*r22*rhoha23+r12*rhoha11*rhodot2*q22*rhoh23*r21*rhoha2&
     &3+r13*rhoha11*rhodot1*q12*rhoh11*r22*rhoha23-r13*rhoha12*rhodot2*q&
     &23*rhoh21*r21*rhoha22-r12*rhoha11*q23*qd22*r21*rhoha23            
      s7 = s8+r11*rhoha12*q11*qd13*r23*rhoha22+r12*rhoha11*q12*qd13*r23*&
     &rhoha21+r12*rhoha13*r23*rhoha21*rhodot1*q11*rhoh12-r13*rhoha12*r23&
     &*rhoha21*rhodot2*q22*rhoh21-r11*rhoha13*rhodot1*q13*rhoh12*r22*rho&
     &ha21+r11*rhoha12*rhodot1*q12*rhoh13*r21*rhoha23-r12*rhoha11*q22*qd&
     &23*r23*rhoha21-r13*rhoha12*rhodot1*q13*rhoh11*r22*rhoha21+r12*rhoh&
     &a13*rhodot1*q13*rhoh11*r22*rhoha21-r12*rhoha11*rhodot2*q21*rhoh23*&
     &r22*rhoha23-r11*rhoha13*q22*qd23*r22*rhoha21+r11*rhoha12*q12*qd13*&
     &r21*rhoha23                                                       
      s8 = s7+r13*rhoha11*q12*qd11*r22*rhoha23+r11*rhoha12*rhodot2*q22*r&
     &hoh23*r23*rhoha21-r13*rhoha12*r23*rhoha21*rhodot1*q11*rhoh12+r13*r&
     &hoha11*rhodot2*q22*rhoh23*r22*rhoha21+r13*rhoha11*rhodot2*q23*rhoh&
     &22*r21*rhoha22-r12*rhoha13*rhodot1*q11*rhoh13*r22*rhoha21-r13*rhoh&
     &a11*rhodot2*q22*rhoh23*r21*rhoha22-r11*rhoha12*rhodot2*q22*rhoh23*&
     &r21*rhoha23-r12*rhoha13*rhodot1*q13*rhoh11*r21*rhoha22+r11*rhoha12&
     &*rhodot2*q23*rhoh21*r23*rhoha22+r13*rhoha11*rhodot2*q22*rhoh21*r23&
     &*rhoha22                                                          
      s6 = s8-r11*rhoha13*rhodot2*q23*rhoh22*r21*rhoha22+r11*rhoha13*rho&
     &dot2*q23*rhoh22*r22*rhoha21-r11*rhoha12*rhodot2*q21*rhoh23*r23*rho&
     &ha22-r12*rhoha13*q11*qd13*r22*rhoha21-r13*rhoha12*rhodot1*q11*rhoh&
     &13*r21*rhoha22+r13*rhoha12*rhodot1*q11*rhoh13*r22*rhoha21+r13*rhoh&
     &a11*rhodot2*q21*rhoh22*r22*rhoha23-r13*rhoha11*rhodot2*q22*rhoh21*&
     &r22*rhoha23+r11*rhoha12*rhodot1*q13*rhoh12*r23*rhoha21+r12*rhoha13&
     &*r21*rhoha23*rhodot1*q12*rhoh11+r11*rhoha12*q13*qd11*r22*rhoha23-r&
     &12*rhoha11*q13*qd12*r23*rhoha21+r11*rhoha13*rhodot1*q12*rhoh11*r23&
     &*rhoha22                                                          
      s8 = s6+r11*rhoha13*q11*qd12*r22*rhoha23-r11*rhoha13*rhodot2*q21*r&
     &hoh22*r22*rhoha23-r11*rhoha12*q23*qd21*r22*rhoha23-r13*rhoha11*rho&
     &dot1*q12*rhoh11*r23*rhoha22+r12*rhoha11*rhodot1*q13*rhoh11*r23*rho&
     &ha22-r12*rhoha11*rhodot2*q23*rhoh22*r21*rhoha23+r12*rhoha11*q13*qd&
     &12*r21*rhoha23+s5-r13*rhoha11*q11*qd12*r22*rhoha23-r12*rhoha13*rho&
     &dot2*q21*rhoh23*r21*rhoha22+r12*rhoha13*r21*rhoha23*rhodot2*q21*rh&
     &oh22                                                              
      s7 = s8-r13*rhoha11*rhodot2*q23*rhoh22*r22*rhoha21-r11*rhoha13*q12&
     &*qd11*r22*rhoha23+r13*rhoha12*r23*rhoha21*rhodot1*q12*rhoh11+r11*r&
     &hoha13*rhodot1*q11*rhoh12*r22*rhoha23+r12*rhoha13*r23*rhoha21*q11*&
     &qd12-r12*rhoha13*r21*rhoha23*q22*qd21+r12*rhoha11*q23*qd21*r22*rho&
     &ha23+r12*rhoha13*r21*rhoha23*q12*qd11-r12*rhoha13*q21*qd23*r21*rho&
     &ha22-r13*rhoha12*r21*rhoha23*q21*qd22-r11*rhoha12*q13*qd12*r21*rho&
     &ha23+r13*rhoha11*q22*qd21*r23*rhoha22+r11*rhoha13*q12*qd13*r22*rho&
     &ha21                                                              
      s8 = s7+r12*rhoha13*r23*rhoha21*q22*qd21+r13*rhoha12*q11*qd13*r22*&
     &rhoha21+r11*rhoha13*q22*qd23*r21*rhoha22-r12*rhoha11*q12*qd13*r21*&
     &rhoha23+r12*rhoha13*q21*qd23*r22*rhoha21-r13*rhoha11*q22*qd21*r22*&
     &rhoha23+r11*rhoha13*q12*qd11*r23*rhoha22+r11*rhoha12*q22*qd23*r23*&
     &rhoha21-r12*rhoha13*q13*qd11*r21*rhoha22-r12*rhoha13*r23*rhoha21*q&
     &21*qd22-r13*rhoha11*q12*qd11*r23*rhoha22+r12*rhoha11*q22*qd23*r21*&
     &rhoha23                                                           
      s4 = s8-r13*rhoha12*r23*rhoha21*q11*qd12-r13*rhoha12*q11*qd13*r21*&
     &rhoha22+r13*rhoha11*q21*qd22*r22*rhoha23+r11*rhoha12*q23*qd22*r21*&
     &rhoha23-r13*rhoha12*q13*qd11*r22*rhoha21-r11*rhoha12*q13*qd11*r23*&
     &rhoha22+r13*rhoha12*q21*qd23*r21*rhoha22+r12*rhoha11*q21*qd23*r23*&
     &rhoha22-r12*rhoha11*rhodot1*q13*rhoh11*r22*rhoha23-r11*rhoha12*q23&
     &*qd22*r23*rhoha21+r11*rhoha13*rhodot2*q22*rhoh21*r22*rhoha23-r11*r&
     &hoha12*q11*qd13*r22*rhoha23+r12*rhoha11*rhodot2*q23*rhoh22*r23*rho&
     &ha21                                                              
      s2 = s3*s4 
      detJ2 = s1*s2 
                                                                        
      s1 = -rho1 
      s3 = rho2 
      s8 = -r11*rhoha12*rhodot1*q11*rhoh13*r22*rhoha23+r11*rhoha12*q23*q&
     &d21*r23*rhoha22-r11*rhoha12*q21*qd23*r23*rhoha22-r12*rhoha13*q23*q&
     &d21*r22*rhoha21+r13*rhoha12*r23*rhoha21*rhodot2*q21*rhoh22-r13*rho&
     &ha12*r23*rhoha21*rhodot2*q22*rhoh21+r13*rhoha11*rhodot2*q21*rhoh22&
     &*r22*rhoha23-r11*rhoha13*rhodot1*q12*rhoh11*r22*rhoha23-r11*rhoha1&
     &3*q22*qd21*r23*rhoha22-r13*rhoha12*q23*qd21*r21*rhoha22+r11*rhoha1&
     &3*rhodot1*q13*rhoh12*r21*rhoha22+r12*rhoha13*r21*rhoha23*rhodot1*q&
     &12*rhoh11                                                         
      s7 = s8-r13*rhoha11*rhodot1*q12*rhoh13*r22*rhoha21-r13*rhoha11*rho&
     &dot1*q13*rhoh12*r21*rhoha22+r13*rhoha11*rhodot1*q13*rhoh12*r22*rho&
     &ha21-r11*rhoha12*q12*qd13*r23*rhoha21-r11*rhoha13*rhodot2*q21*rhoh&
     &22*r22*rhoha23+r11*rhoha13*rhodot2*q22*rhoh21*r22*rhoha23+r11*rhoh&
     &a13*q12*qd11*r23*rhoha22+r11*rhoha13*rhodot2*q21*rhoh22*r23*rhoha2&
     &2+r11*rhoha12*q13*qd12*r23*rhoha21-r12*rhoha13*r21*rhoha23*q22*qd2&
     &1+r13*rhoha11*rhodot1*q11*rhoh12*r23*rhoha22+r13*rhoha11*rhodot2*q&
     &22*rhoh21*r23*rhoha22                                             
      s8 = s7-r11*rhoha13*rhodot2*q23*rhoh22*r21*rhoha22+r11*rhoha12*q22&
     &*qd23*r23*rhoha21-r13*rhoha12*r21*rhoha23*rhodot1*q12*rhoh11-r12*r&
     &hoha11*rhodot1*q12*rhoh13*r21*rhoha23-r12*rhoha11*rhodot1*q13*rhoh&
     &12*r23*rhoha21-r13*rhoha11*rhodot2*q21*rhoh22*r23*rhoha22+r11*rhoh&
     &a13*q21*qd22*r23*rhoha22+r11*rhoha12*rhodot2*q22*rhoh23*r23*rhoha2&
     &1+r11*rhoha12*rhodot1*q13*rhoh11*r22*rhoha23-r11*rhoha13*rhodot1*q&
     &13*rhoh12*r22*rhoha21+r12*rhoha13*r21*rhoha23*q21*qd22            
      s6 = s8+r12*rhoha11*q22*qd23*r21*rhoha23-r12*rhoha13*rhodot1*q11*r&
     &hoh13*r22*rhoha21-r12*rhoha13*r23*rhoha21*q12*qd11+r13*rhoha11*q11&
     &*qd12*r23*rhoha22-r13*rhoha11*q22*qd23*r21*rhoha22-r11*rhoha12*q13&
     &*qd11*r23*rhoha22-r12*rhoha11*rhodot2*q23*rhoh21*r23*rhoha22-r11*r&
     &hoha13*q21*qd22*r22*rhoha23+r12*rhoha11*q12*qd13*r23*rhoha21-r12*r&
     &hoha13*rhodot2*q23*rhoh21*r22*rhoha21-r12*rhoha13*r21*rhoha23*q11*&
     &qd12+r12*rhoha13*r21*rhoha23*q12*qd11+r11*rhoha12*rhodot2*q21*rhoh&
     &23*r22*rhoha23                                                    
      s8 = s6-r12*rhoha13*r23*rhoha21*q21*qd22-r11*rhoha13*q11*qd12*r23*&
     &rhoha22-r13*rhoha12*r21*rhoha23*rhodot2*q21*rhoh22-r13*rhoha12*q11&
     &*qd13*r21*rhoha22+r13*rhoha12*q11*qd13*r22*rhoha21-r12*rhoha11*q13&
     &*qd11*r22*rhoha23-r13*rhoha12*q13*qd11*r22*rhoha21-r13*rhoha12*r23&
     &*rhoha21*rhodot1*q11*rhoh12+r13*rhoha11*rhodot2*q22*rhoh23*r22*rho&
     &ha21+r13*rhoha11*rhodot2*q23*rhoh22*r21*rhoha22-r11*rhoha12*q11*qd&
     &13*r22*rhoha23                                                    
      s7 = s8-r11*rhoha12*rhodot2*q23*rhoh21*r22*rhoha23+r12*rhoha11*q11&
     &*qd13*r22*rhoha23+r12*rhoha11*q23*qd22*r23*rhoha21-r12*rhoha11*q13&
     &*qd12*r23*rhoha21+r12*rhoha11*q13*qd12*r21*rhoha23+r11*rhoha12*q12&
     &*qd13*r21*rhoha23+r11*rhoha12*q13*qd11*r22*rhoha23+r13*rhoha12*r21&
     &*rhoha23*q22*qd21+r12*rhoha13*r23*rhoha21*q22*qd21+r11*rhoha13*rho&
     &dot2*q23*rhoh22*r22*rhoha21+r11*rhoha12*q11*qd13*r23*rhoha22+r12*r&
     &hoha13*rhodot1*q11*rhoh13*r21*rhoha22                             
      s8 = s7+r12*rhoha11*rhodot2*q22*rhoh23*r21*rhoha23+r12*rhoha11*rho&
     &dot2*q23*rhoh22*r23*rhoha21-r12*rhoha11*rhodot2*q23*rhoh22*r21*rho&
     &ha23-r13*rhoha12*rhodot2*q23*rhoh21*r21*rhoha22+r12*rhoha13*r23*rh&
     &oha21*q11*qd12+r13*rhoha12*r21*rhoha23*q11*qd12+r11*rhoha12*rhodot&
     &2*q23*rhoh22*r21*rhoha23-r11*rhoha12*rhodot1*q12*rhoh13*r23*rhoha2&
     &1+r12*rhoha13*rhodot2*q21*rhoh23*r22*rhoha21+r13*rhoha11*q22*qd21*&
     &r23*rhoha22+r13*rhoha12*r23*rhoha21*rhodot1*q12*rhoh11+r11*rhoha13&
     &*rhodot1*q11*rhoh12*r22*rhoha23                                   
      s5 = s8+r11*rhoha13*rhodot1*q12*rhoh11*r23*rhoha22+r11*rhoha13*q22&
     &*qd23*r21*rhoha22-r12*rhoha11*rhodot1*q11*rhoh13*r23*rhoha22+r12*r&
     &hoha11*q23*qd21*r22*rhoha23-r12*rhoha13*q13*qd11*r21*rhoha22+r13*r&
     &hoha11*q13*qd12*r22*rhoha21-r13*rhoha11*rhodot2*q22*rhoh21*r22*rho&
     &ha23-r11*rhoha12*rhodot2*q22*rhoh23*r21*rhoha23-r11*rhoha12*rhodot&
     &2*q23*rhoh22*r23*rhoha21+r13*rhoha12*r23*rhoha21*q12*qd11+r13*rhoh&
     &a12*rhodot1*q13*rhoh11*r21*rhoha22-r13*rhoha11*q11*qd12*r22*rhoha2&
     &3-r12*rhoha13*rhodot2*q21*rhoh23*r21*rhoha22                      
      s8 = r13*rhoha11*q23*qd22*r21*rhoha22-r13*rhoha11*q12*qd11*r23*rho&
     &ha22-r11*rhoha12*q22*qd23*r21*rhoha23-r11*rhoha12*q23*qd22*r23*rho&
     &ha21-r12*rhoha13*q21*qd23*r21*rhoha22+r12*rhoha13*r23*rhoha21*rhod&
     &ot2*q22*rhoh21-r12*rhoha11*q21*qd23*r22*rhoha23-r11*rhoha13*q13*qd&
     &12*r22*rhoha21+r12*rhoha11*rhodot2*q23*rhoh21*r22*rhoha23-r12*rhoh&
     &a11*rhodot2*q21*rhoh23*r22*rhoha23+r12*rhoha13*r21*rhoha23*rhodot2&
     &*q21*rhoh22-r13*rhoha11*rhodot2*q23*rhoh22*r22*rhoha21            
      s7 = s8-r11*rhoha12*q23*qd21*r22*rhoha23+r11*rhoha12*q21*qd23*r22*&
     &rhoha23+r13*rhoha12*q21*qd23*r21*rhoha22-r13*rhoha11*rhodot1*q12*r&
     &hoh11*r23*rhoha22-r13*rhoha11*q23*qd22*r22*rhoha21-r12*rhoha11*rho&
     &dot2*q22*rhoh23*r23*rhoha21-r12*rhoha11*rhodot1*q13*rhoh11*r22*rho&
     &ha23-r13*rhoha12*r21*rhoha23*q12*qd11-r13*rhoha12*r23*rhoha21*q22*&
     &qd21+r13*rhoha12*q23*qd21*r22*rhoha21+r11*rhoha13*rhodot2*q22*rhoh&
     &23*r21*rhoha22-r12*rhoha11*q23*qd22*r21*rhoha23                   
      s8 = s7+r11*rhoha12*rhodot1*q12*rhoh13*r21*rhoha23+r11*rhoha12*rho&
     &dot1*q11*rhoh13*r23*rhoha22+r12*rhoha11*rhodot2*q21*rhoh23*r23*rho&
     &ha22+r11*rhoha13*q23*qd22*r22*rhoha21-r11*rhoha13*q22*qd23*r22*rho&
     &ha21-r11*rhoha13*q23*qd22*r21*rhoha22-r13*rhoha12*rhodot1*q13*rhoh&
     &11*r22*rhoha21+r13*rhoha12*rhodot2*q23*rhoh21*r22*rhoha21+r13*rhoh&
     &a12*rhodot2*q21*rhoh23*r21*rhoha22-r13*rhoha12*r23*rhoha21*q11*qd1&
     &2+r12*rhoha11*rhodot1*q11*rhoh13*r22*rhoha23                      
      s6 = s8+r12*rhoha13*q21*qd23*r22*rhoha21-r11*rhoha12*rhodot1*q13*r&
     &hoh12*r21*rhoha23+r12*rhoha13*rhodot2*q23*rhoh21*r21*rhoha22-r13*r&
     &hoha11*q22*qd21*r22*rhoha23+r12*rhoha11*rhodot1*q12*rhoh13*r23*rho&
     &ha21-r11*rhoha13*rhodot1*q11*rhoh12*r23*rhoha22+r13*rhoha12*r21*rh&
     &oha23*rhodot1*q11*rhoh12-r11*rhoha13*q12*qd13*r21*rhoha22-r11*rhoh&
     &a13*rhodot2*q22*rhoh23*r22*rhoha21+r13*rhoha12*r21*rhoha23*rhodot2&
     &*q22*rhoh21+r11*rhoha13*rhodot1*q12*rhoh13*r22*rhoha21+r12*rhoha13&
     &*q13*qd11*r22*rhoha21+r12*rhoha13*r23*rhoha21*rhodot1*q11*rhoh12  
      s8 = s6-r12*rhoha13*r23*rhoha21*rhodot1*q12*rhoh11+r12*rhoha11*rho&
     &dot1*q13*rhoh12*r21*rhoha23-r12*rhoha11*q12*qd13*r21*rhoha23+r13*r&
     &hoha11*q21*qd22*r22*rhoha23-r11*rhoha13*q12*qd11*r22*rhoha23+r13*r&
     &hoha11*q12*qd13*r21*rhoha22-r12*rhoha13*r21*rhoha23*rhodot2*q22*rh&
     &oh21-r12*rhoha13*r21*rhoha23*rhodot1*q11*rhoh12+r11*rhoha12*rhodot&
     &2*q23*rhoh21*r23*rhoha22-r11*rhoha12*rhodot2*q21*rhoh23*r23*rhoha2&
     &2-r12*rhoha11*q23*qd21*r23*rhoha22                                
      s7 = s8+r12*rhoha11*q21*qd23*r23*rhoha22+r12*rhoha11*q13*qd11*r23*&
     &rhoha22-r12*rhoha11*q11*qd13*r23*rhoha22-r11*rhoha12*rhodot1*q13*r&
     &hoh11*r23*rhoha22+r13*rhoha11*rhodot1*q12*rhoh11*r22*rhoha23+r11*r&
     &hoha12*q23*qd22*r21*rhoha23-r11*rhoha12*q13*qd12*r21*rhoha23-r13*r&
     &hoha11*rhodot1*q11*rhoh12*r22*rhoha23+r11*rhoha13*q12*qd13*r22*rho&
     &ha21-r13*rhoha12*q21*qd23*r22*rhoha21+r13*rhoha12*q13*qd11*r21*rho&
     &ha22+r12*rhoha13*rhodot1*q13*rhoh11*r22*rhoha21-r11*rhoha13*rhodot&
     &2*q22*rhoh21*r23*rhoha22                                          
      s8 = s7+s5+r13*rhoha12*r23*rhoha21*q21*qd22-r12*rhoha13*q11*qd13*r&
     &22*rhoha21+r12*rhoha13*q23*qd21*r21*rhoha22+r11*rhoha13*q13*qd12*r&
     &21*rhoha22-r11*rhoha13*rhodot1*q12*rhoh13*r21*rhoha22+r11*rhoha12*&
     &rhodot1*q13*rhoh12*r23*rhoha21-r13*rhoha11*q12*qd13*r22*rhoha21-r1&
     &3*rhoha11*q13*qd12*r21*rhoha22+r13*rhoha11*rhodot1*q12*rhoh13*r21*&
     &rhoha22-r12*rhoha13*r23*rhoha21*rhodot2*q21*rhoh22-r13*rhoha12*r21&
     &*rhoha23*q21*qd22                                                 
      s4 = s8-r13*rhoha12*rhodot2*q21*rhoh23*r22*rhoha21+r12*rhoha11*rho&
     &dot1*q13*rhoh11*r23*rhoha22+r12*rhoha13*q11*qd13*r21*rhoha22-r12*r&
     &hoha11*q22*qd23*r23*rhoha21-r13*rhoha12*rhodot1*q11*rhoh13*r21*rho&
     &ha22+r13*rhoha12*rhodot1*q11*rhoh13*r22*rhoha21+r13*rhoha11*q22*qd&
     &23*r22*rhoha21-r13*rhoha11*rhodot2*q22*rhoh23*r21*rhoha22+r11*rhoh&
     &a13*q22*qd21*r22*rhoha23+r11*rhoha13*q11*qd12*r22*rhoha23+r13*rhoh&
     &a11*q12*qd11*r22*rhoha23-r12*rhoha13*rhodot1*q13*rhoh11*r21*rhoha2&
     &2-r13*rhoha11*q21*qd22*r23*rhoha22                                
      s2 = s3*s4 
      detJ2 = s1*s2 
                                                                        
      detK2 = rho1*rho2**2*(-rhohd22*r23*rhoha21-r21*rhohd23*rhoha22+r21&
     &*rhohd22*rhoha23+r23*rhohd21*rhoha22-r22*rhohd21*rhoha23+rhohd23*r&
     &22*rhoha21)*(-r21*r12*rhoha13+r21*r13*rhoha12+r12*rhoha11*r23+r11*&
     &rhoha13*r22-r11*rhoha12*r23-r13*rhoha11*r22)                      
                                                                        
      detJ3 = -rho1**2*(-rhoha13*r11*rhohd12+rhoha11*r13*rhohd12+r11*rho&
     &ha12*rhohd13-r13*rhoha12*rhohd11+r12*rhoha13*rhohd11-r12*rhoha11*r&
     &hohd13)*(-r11*q22*qd23+r11*q23*qd22+r11*q12*qd13-r11*q13*qd12-r11*&
     &rhodot2*q22*rhoh23+r11*rhodot2*q23*rhoh22+r11*rhodot1*q12*rhoh13-r&
     &11*rhodot1*q13*rhoh12-r12*q11*qd13-r12*rhodot2*q23*rhoh21-r13*q12*&
     &qd11-r13*rhodot2*q21*rhoh22+r13*q11*qd12+r12*q13*qd11+r13*q22*qd21&
     &+r13*rhodot2*q22*rhoh21+r13*rhodot1*q11*rhoh12+r12*rhodot1*q13*rho&
     &h11-r12*rhodot1*q11*rhoh13-r13*rhodot1*q12*rhoh11+r12*rhodot2*q21*&
     &rhoh23-r12*q23*qd21+r12*q21*qd23-r13*q21*qd22)                    
                                                                        
      detK3 = rho1**2*rho2*(-rhoha13*r11*rhohd12+rhoha11*r13*rhohd12+r11&
     &*rhoha12*rhohd13-r13*rhoha12*rhohd11+r12*rhoha13*rhohd11-r12*rhoha&
     &11*rhohd13)*(r11*r22*rhohd23-r11*r23*rhohd22+r12*r23*rhohd21+r13*r&
     &21*rhohd22-r12*r21*rhohd23-r13*r22*rhohd21)                       

    END SUBROUTINE determinants

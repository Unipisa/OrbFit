* Copyright (C) 2000 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: May 25, 2000
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                         E L E K D T                           *
*  *                                                               *
*  *        Time advance of a set of keplerian elements            *
*  *                                                               *
*  *****************************************************************
*
* INPUT:    ELEM      -  Orbital elements
*           TYPE      -  Element type
*           GM        -  G * ( M + m )
*           DT        -  Time interval
*
* OUTPUT:   ELEM      -  Updated orbital elements
*
      SUBROUTINE elekdt(elem,type,gm,dt)
      IMPLICIT NONE

      DOUBLE PRECISION elem(6),gm,dt
      CHARACTER*(*) type

      DOUBLE PRECISION enne,sma,anec,emme,sqe1,tanf,p,dcap,t0,delta
      DOUBLE PRECISION effe,esse,sinthe,costhe,u3,theta,sinf,cosf,b
      DOUBLE PRECISION beta2,beta,r,coshf,sinhf,expx,expmx

      DOUBLE PRECISION anecc,hypan
      EXTERNAL anecc,hypan

      IF(type.EQ.'KEP') THEN
          enne=SQRT(gm/elem(1)**3)
          elem(6)=elem(6)+enne*dt
      ELSEIF(type.EQ.'COM') THEN
          IF(elem(2).LT.1.D0) THEN
              sma=elem(1)/(1.d0-elem(2))
              enne=SQRT(gm/sma**3)
              sqe1=SQRT((1.d0-elem(2))/(1.d0+elem(2)))
              anec=sqe1*TAN(elem(6)/2.d0)
              anec=2.d0*ATAN(anec)
              emme=anec-elem(2)*SIN(anec)+enne*dt
              anec=anecc(emme,elem(2))
              tanf=TAN(anec/2)/sqe1
              elem(6)=2*ATAN(tanf)
          ELSEIF(elem(2).EQ.1.D0) THEN
              u3=1.d0/3.d0
              p=2.d0*elem(1)
              enne=SQRT(gm/p**3)
              dcap=TAN(elem(6)/2)
              t0=-dcap*(1+dcap**2/3)/(2*enne)
              delta=dt-t0
              IF(delta.EQ.0.d0) THEN
                  effe=0.d0
              ELSE
                  esse=3*enne*delta
                  esse=ATAN(1/esse)
                  sinthe=SIN(esse/2)
                  costhe=COS(esse/2)
                  IF(sinthe.GT.0.d0) THEN
                      sinthe=sinthe**u3
                  ELSEIF(sinthe.LT.0.d0) THEN
                      sinthe=-(ABS(sinthe)**u3)
                  END IF
                  IF(costhe.GT.0.d0) THEN
                      costhe=costhe**u3
                  ELSEIF(costhe.LT.0.d0) THEN
                      costhe=-(ABS(costhe)**u3)
                  END IF
                  theta=ATAN2(sinthe,costhe)
                  sinf=2.d0*COS(2*theta)
                  cosf=SIN(2*theta)
                  elem(6)=2*ATAN2(sinf,cosf)
              END IF
          ELSE
              p=(1.d0+elem(2))*elem(1)
              sma=elem(1)/(elem(2)-1.d0)
              beta2=elem(2)**2-1.d0
              beta=SQRT(beta2)
              enne=SQRT(gm/sma**3)
              b=SQRT(sma*p)
              IF(elem(6).EQ.0.d0) THEN
                  t0=0.d0
              ELSE
                  cosf=COS(elem(6))
                  sinf=SIN(elem(6))
                  r=p/(1.d0+elem(2)*cosf)
                  coshf=elem(2)-r*cosf/sma
                  sinhf=r*sinf/b
                  expx=coshf+sinhf
                  expmx=coshf-sinhf
                  IF(expx.GT.expmx) THEN
                      effe=LOG(expx)
                  ELSE
                      effe=-LOG(expmx)
                  END IF
                  t0=elem(2)*SINH(effe)-effe
                  t0=-t0/enne
              END IF
              emme=enne*(dt-t0)
              anec=hypan(emme,elem(2))
              tanf=SQRT((elem(2)+1)/(elem(2)-1))*TANH(anec/2)
              elem(6)=2*ATAN(tanf)
          END IF
      ELSE
          STOP '**** elekdt: internal error (01) ****'
      END IF

      END

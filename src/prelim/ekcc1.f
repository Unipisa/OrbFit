* Copyright (C) 2000 by Mario Carpino (carpino@brera.mi.astro.it)
* Version: May 29, 2000
* ---------------------------------------------------------------------
*
*  *****************************************************************
*  *                                                               *
*  *                          C C E K 1                            *
*  *                                                               *
*  *            Transformation from keplerian elements             *
*  *               to cartesian coordinates (pos/vel)              *
*  *                                                               *
*  *****************************************************************
*
*     see: A.E.Roy, Orbital Motion, Adam Hilger Ltd, Chapter 4
*
* INPUT:    ELEM(6)   -  Orbital elements
*                            ELEM(1) = a/q
*                            ELEM(2) = e
*                            ELEM(3) = i
*                            ELEM(4) = long. node
*                            ELEM(5) = arg.pericenter
*                            ELEM(6) = M/f
*           TYPE      -  Type of orbital elements (KEP/COM)
*           GM        -  G * ( M + m )
*           DT        -  Time interval between the epoch of orbital
*                        elements and the epoch at which cartesian
*                        coordinates are to be computed
*
* OUTPUT:   XV(6)     -  Position/velocity vector
*
      SUBROUTINE ekcc1(elem,type,xv,gm,dt)
      IMPLICIT NONE

      DOUBLE PRECISION elem(6),xv(6),gm,dt
      CHARACTER*(*) type

      DOUBLE PRECISION x1(3),v1(3),ra(3,3),rb(3,3),rc(3,3),rd(3,3)
      DOUBLE PRECISION rot(3,3),u3,ecc,ainc,argper,anod,sma,emme0
      DOUBLE PRECISION anec0,enne,emme,anec,b,cose,sine,danec,q,p
      DOUBLE PRECISION effe0,dcap,t0,delta,effe,esse,sinthe,costhe
      DOUBLE PRECISION theta,sinf,cosf,r,vel,fnorm,beta2,beta
      DOUBLE PRECISION coshf,sinhf,expx,expmx

      INCLUDE 'trig.h'

      DOUBLE PRECISION anecc,hypan
      EXTERNAL anecc,hypan

      u3=1.d0/3.d0

      ecc=elem(2)
      ainc=elem(3)
      argper=elem(5)
      anod=elem(4)

* Coordinates in the orbital plane: elliptic orbit
      IF(ecc.LT.1.d0) THEN
          IF(type.EQ.'KEP') THEN
              sma=elem(1)
              emme0=elem(6)
          ELSEIF(type.EQ.'COM') THEN
              sma=elem(1)/(1.d0-ecc)
              anec0=SQRT((1.d0-ecc)/(1.d0+ecc))*TAN(elem(6)/2.d0)
              anec0=2.d0*ATAN(anec0)
              emme0=anec0-ecc*SIN(anec0)
          ELSE
              STOP '**** ekcc1: TYPE = ??? ****'
          END IF
          enne=SQRT(gm/sma**3)
          emme=emme0+enne*dt
          emme=dmod(emme,dpig)
          anec=anecc(emme,ecc)
          b=SQRT(1.d0-ecc**2)*sma
          cose=COS(anec)
          sine=SIN(anec)
          x1(1)=sma*(cose-ecc)
          x1(2)=b*sine
          x1(3)=0.d0
          danec=enne/(1.d0-ecc*cose)
          v1(1)=-sma*sine*danec
          v1(2)=b*cose*danec
          v1(3)=0.d0

* Coordinates in the orbital plane: parabolic orbit
      ELSEIF(ecc.EQ.1.d0) THEN
          IF(type.EQ.'COM') THEN
              q=elem(1)
              effe0=elem(6)
          ELSEIF(type.EQ.'KEP') THEN
              STOP '**** ekcc1: TYPE=''KEP'' not valid if ECC=1 ****'
          ELSE
              STOP '**** ekcc1: TYPE = ??? ****'
          END IF

* T0 (epoch of perihelion passage)
          p=2.d0*q
          enne=SQRT(gm/p**3)
          dcap=TAN(effe0/2.d0)
          t0=-dcap*(1.d0+dcap**2/3.d0)/(2.d0*enne)

* True anomaly
          delta=dt-t0
          IF(delta.EQ.0.d0) THEN
              effe=0.d0
          ELSE
              esse=3.d0*enne*delta
              esse=ATAN(1.d0/esse)
              sinthe=SIN(esse/2.d0)
              costhe=COS(esse/2.d0)
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
              sinf=2.d0*COS(2.d0*theta)
              cosf=SIN(2.d0*theta)
              effe=2.d0*ATAN2(sinf,cosf)
          END IF

* Coordinates in the orbital plane
          cosf=COS(effe)
          r=p/(1.d0+cosf)
          x1(1)=r*cosf
          x1(2)=r*SIN(effe)
          x1(3)=0.d0
          vel=SQRT(2.d0*gm/r)
          v1(1)=-x1(2)
          v1(2)=p
          fnorm=SQRT(v1(1)**2+v1(2)**2)
          fnorm=vel/fnorm
          v1(1)=fnorm*v1(1)
          v1(2)=fnorm*v1(2)
          v1(3)=0.d0

* Coordinates in the orbital plane: hyperbolic orbit
      ELSE
          IF(type.EQ.'COM') THEN
              q=elem(1)
              effe0=elem(6)
          ELSEIF(type.EQ.'KEP') THEN
              STOP '**** ekcc1: TYPE=''KEP'' not valid if ECC>1 ****'
          ELSE
              STOP '**** ekcc1: TYPE = ??? ****'
          END IF

* T0 (epoch of perihelion passage)
          p=(1.d0+ecc)*q
          sma=q/(ecc-1.d0)
          beta2=ecc**2-1.d0
          beta=SQRT(beta2)
          enne=SQRT(gm/sma**3)
          b=SQRT(sma*p)
          IF(effe0.EQ.0.d0) THEN
              t0=0.d0
          ELSE
              cosf=COS(effe0)
              sinf=SIN(effe0)
              r=p/(1.d0+ecc*cosf)
              coshf=ecc-r*cosf/sma
              sinhf=r*sinf/b
              expx=coshf+sinhf
              expmx=coshf-sinhf
              IF(expx.GT.expmx) THEN
                  effe=LOG(expx)
              ELSE
                  effe=-LOG(expmx)
              END IF
              t0=ecc*SINH(effe)-effe
              t0=-t0/enne
          END IF

* Hyperbolic eccentric anomaly
          emme=enne*(dt-t0)
          anec=hypan(emme,ecc)

* Coordinates in the orbital plane
          coshf=COSH(anec)
          sinhf=SINH(anec)
          r=sma*(ecc*coshf-1.d0)
          x1(1)=sma*(ecc-coshf)
          x1(2)=b*sinhf
          x1(3)=0.d0
          vel=2.d0/r+1.d0/sma
          vel=SQRT(gm*vel)
          v1(1)=-sinhf
          v1(2)=beta*coshf
          fnorm=SQRT(v1(1)**2+v1(2)**2)
          fnorm=vel/fnorm
          v1(1)=fnorm*v1(1)
          v1(2)=fnorm*v1(2)
          v1(3)=0.d0
      END IF

* Rotation of coordinates in the inertial frame
      CALL rotmt(-argper,ra,3)
      CALL rotmt(-ainc,rb,1)
      CALL rotmt(-anod,rc,3)
      CALL prodmm(rd,rb,ra)
      CALL prodmm(rot,rc,rd)
      CALL prodmv(xv(1),rot,x1)
      CALL prodmv(xv(4),rot,v1)

      END

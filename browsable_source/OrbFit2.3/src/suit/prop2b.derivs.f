c
c  PROP2B
c
c  Task:       Solves the two body problem in equinoctal coordinates
c
c  Features:   the Newton's method is exploited in the interval
c              [W, 2PI + W] where W is the peri-planet longitude:
c              W = w_p + RAAN
c
c  Input arguments:
c              E         =     equinoctal orbital elements 
c                               e=(a,h,k,p,q,lambda)
c              T0        =     epoch time
c              T1        =     prediction time
c              GM        =     gravitational constant of the system GM
c              IDER       =     flag for derivatives option :
c                                0 = only position and velocites
c                                1 = first derivatives
c                                2 = first and second partial derivatives 
c  Output arguments:
c              X         =     position and velocity components
c                               in absolute cartesian coordinates
c              DXDE      =     first derivatives of position vector x
c                               with respect to elements
c              DDXDE     =     second derivatives of position vector x   
c                               with respect to elements 
c
c****************
c   static memory not required
c****************
      subroutine prop2b(t0,e,t1,x,gm,ider,dxde,ddxde)
      implicit double precision (a-h,o-z)
      parameter (iter=25)
      dimension e(6),x(6)
      dimension f(3),g(3),w(3)
      dimension dxde(6,6),ddxde(3,6,6)
c  Trgonometric constants 
      include 'trig.h'
c  Mean motion 
      enne=sqrt(gm/e(1)**3)
c  Mean longitude PML at time T1
      pml=e(6)+enne*(t1-t0)
c  Absolute peri-planet longitude POL at epoch T0
      ecc2=e(2)**2+e(3)**2
      errm=roff(n)
      eps=errm*1.d2
      if(ecc2.lt.eps)then
         pol=0.d0
      elseif(ecc2.ge.1.d0)then
         write(*,*)' dcc.ge.1, ecc**2=',ecc2
         stop
      else
         pol=atan2(e(2),e(3))
         pol=princ(pol)
      endif
c  Mean longitude restriction to [POL, POL + 2*PIGR]
c  (Newton's method continuity conditions)
      pml = princ(pml)
      if (pml.lt.pol)then
         pml = pml + 2.d0*pig
      endif
c  Newton's method for motion equation solution:
c  R(F,lambda) = F - ksinF + hcosF - lambda = 0
c  search for F for a given lambda within [POL, POL + 2PIGR]
      el = pig+pol
      do 70 j=1,iter
             sinel=sin(el)
             cosel=cos(el)
             rf = el - e(3)*sinel + e(2)*cosel - pml
             rdf = 1.d0 - e(3)*cosel - e(2)*sinel
             del = - rf/rdf
             el = el + del
             if (abs(del).lt.eps) goto 100
 70   continue
      write(*,*)' Too many iter. in newton, iter=',iter,' del=',del
      write(*,*) ' eq,eps ',e, eps
      stop
c  Computation of position and velocity on the orbit plane
c  in equinoctal cartesian coordinates (f,g,w)
 100  beta=1.d0/(1.d0+sqrt(1.d0-ecc2))
      xe=e(1)*((1.d0-beta*e(2)**2)*cosel+e(2)*e(3)*beta*sinel-e(3))
      ye=e(1)*((1.d0-beta*e(3)**2)*sinel+e(2)*e(3)*beta*cosel-e(2))
c  Equinoctal reference frame
      upq=1.d0+e(4)**2+e(5)**2
      f(1)=(1.d0-e(4)**2+e(5)**2)/upq
      f(2)=2.d0*e(4)*e(5)/upq
      f(3)=-2.d0*e(4)/upq
      g(1)=2.d0*e(4)*e(5)/upq
      g(2)=(1.d0+e(4)**2-e(5)**2)/upq
      g(3)=2.d0*e(5)/upq
c  Conversion from equinoctal to absolute coordinates
      call lincom(f,xe,g,ye,x)
c  Computation of velocities 
      coe=enne*e(1)**2/sqrt(xe**2+ye**2)
      xpe=coe*(e(2)*e(3)*beta*cosel-(1.d0-beta*e(2)**2)*sinel)
      ype=coe*((1.d0-beta*e(3)**2)*cosel-e(2)*e(3)*beta*sinel)
      call lincom(f,xpe,g,ype,x(4))
c  Computation of partials if required
      if(ider.lt.1)return
c  Equinoctal reference frame, third vector 
      w(1)=2.d0*e(4)/upq
      w(2)=-2.d0*e(5)/upq
      w(3)=(1.d0-e(4)**2-e(5)**2)/upq
c  Computation of same temporary variables
      r=sqrt(xe**2+ye**2)
      tmp1=pml-el
      tmp2=beta+e(2)**2*beta**3/(1.d0-beta)
      tmp3=e(2)*e(3)*beta**3/(1.d0-beta)
      tmp4=beta*e(2)-sinel
      tmp5=beta*e(3)-cosel
      tmp6=beta+e(3)**2*beta**3/(1.d0-beta)
      tmp7=1.d0-r/e(1)
      tmp8=sinel-e(2)
      tmp9=cosel-e(3)
      tmp10=e(1)*cosel/r
      tmp11=e(1)*sinel/r
      tmp12=enne*e(1)**2/r
c  Computation of derivatives of position vector w. r. to the elements 
      dxde(1,1)=(x(1)-3.d0*x(4)*(t1-t0)/2.d0)/e(1)
      dxde(2,1)=(x(2)-3.d0*x(5)*(t1-t0)/2.d0)/e(1)
      dxde(3,1)=(x(3)-3.d0*x(6)*(t1-t0)/2.d0)/e(1)
      dx1de2=-e(1)*(tmp1*tmp2+e(1)*cosel*tmp4/r)          
      dx2de2=e(1)*(tmp1*tmp3-1.d0+e(1)*cosel*tmp5/r)
      call lincom(f,dx1de2,g,dx2de2,dxde(1,2)) 
      dx1de3=-e(1)*(tmp1*tmp3+1.d0-e(1)*sinel*tmp4/r)
      dx2de3=e(1)*(tmp1*tmp6-e(1)*sinel*tmp5/r)
      call lincom(f,dx1de3,g,dx2de3,dxde(1,3))
      dxde(1,4)=2.d0*(e(5)*(ye*f(1)-xe*g(1))-xe*w(1))/upq
      dxde(2,4)=2.d0*(e(5)*(ye*f(2)-xe*g(2))-xe*w(2))/upq
      dxde(3,4)=2.d0*(e(5)*(ye*f(3)-xe*g(3))-xe*w(3))/upq
      dxde(1,5)=2.d0*(e(4)*(-ye*f(1)+xe*g(1))+ye*w(1))/upq          
      dxde(2,5)=2.d0*(e(4)*(-ye*f(2)+xe*g(2))+ye*w(2))/upq
      dxde(3,5)=2.d0*(e(4)*(-ye*f(3)+xe*g(3))+ye*w(3))/upq
      dxde(1,6)=x(4)/enne
      dxde(2,6)=x(5)/enne
      dxde(3,6)=x(6)/enne
c  Computation of derivatives of velocity vector w. r. to the elements
      dxde(4,1)=-(x(4)-3.d0*gm*x(1)*(t1-t0)/r**3)/(2.d0*e(1))
      dxde(5,1)=-(x(5)-3.d0*gm*x(2)*(t1-t0)/r**3)/(2.d0*e(1))
      dxde(6,1)=-(x(6)-3.d0*gm*x(3)*(t1-t0)/r**3)/(2.d0*e(1))
      dx4de2=tmp12*(tmp7*tmp2+e(1)**2*tmp8*tmp4/r**2+tmp10*cosel)          
      dx5de2=-tmp12*(tmp7*tmp3+e(1)**2*tmp8*tmp5/r**2-tmp10*sinel)
      call lincom(f,dx4de2,g,dx5de2,dxde(4,2))
      dx4de3=tmp12*(tmp7*tmp3+e(1)**2*tmp9*tmp4/r**2-tmp11*cosel)
      dx5de3=-tmp12*(tmp7*tmp6+e(1)**2*tmp9*tmp5/r**2+tmp11*sinel)
      call lincom(f,dx4de3,g,dx5de3,dxde(4,3))
      dxde(4,4)=2.d0*(e(5)*(ype*f(1)-xpe*g(1))-xpe*w(1))/upq
      dxde(5,4)=2.d0*(e(5)*(ype*f(2)-xpe*g(2))-xpe*w(2))/upq
      dxde(6,4)=2.d0*(e(5)*(ype*f(3)-xpe*g(3))-xpe*w(3))/upq
      dxde(4,5)=2.d0*(e(4)*(-ype*f(1)+xpe*g(1))+ype*w(1))/upq          
      dxde(5,5)=2.d0*(e(4)*(-ype*f(2)+xpe*g(2))+ype*w(2))/upq
      dxde(6,5)=2.d0*(e(4)*(-ype*f(3)+xpe*g(3))+ype*w(3))/upq
      dxde(4,6)=-enne*e(1)**3*x(1)/r**3
      dxde(5,6)=-enne*e(1)**3*x(2)/r**3
      dxde(6,6)=-enne*e(1)**3*x(3)/r**3
c  Computation of second derivatives if required
      if(ider.lt.2)return
      e1=e(1)
      e2=e(2)
      e3=e(3)
      e4=e(4)
      e5=e(5)
      e6=e(6)
c  Begins maple computation
      t1872 = e2**2
      t1873 = e3**2
      t1875 = sqrt(1.d0-t1872-t1873)
      t1876 = 1.d0+t1875
      t1877 = 1/t1876
      t1879 = 1.d0-t1877*t1872
      t1880 = cos(el)
      t1881 = t1879*t1880
      t1882 = e2*e3
      t1883 = sin(el)
      t1885 = t1882*t1877*t1883
      t1886 = t1881+t1885-e3
      t1887 = e1*t1886
      t1888 = e4**2
      t1889 = e5**2
      t1890 = 1.d0+t1888+t1889
      t1891 = 1/t1890
      t1892 = e4*t1891
      t1893 = t1887*t1892
      t1894 = 1.d0-t1888+t1889
      t1895 = t1890**2
      t1896 = 1/t1895
      t1897 = t1894*t1896
      t1898 = t1897*e4
      t1901 = 1.d0-t1877*t1873
      t1902 = t1901*t1883
      t1904 = t1882*t1877*t1880
      t1905 = t1902+t1904-e2
      t1906 = e1*t1905
      t1907 = e5*t1891
      t1908 = t1906*t1907
      t1909 = t1888*e5
      t1910 = t1909*t1896
      t1911 = t1906*t1910
      t1912 = sqrt(gm)
      t1913 = sqrt(e1)
      t1914 = t1912*t1913
      t1915 = e1**2
      t1916 = t1886**2
      t1918 = t1905**2
      t1920 = t1915*t1916+t1915*t1918
      t1921 = sqrt(t1920)
      t1922 = 1/t1921
      t1923 = t1914*t1922
      t1925 = t1904-t1879*t1883
      t1926 = t1925*e4
      t1927 = t1926*t1891
      t1928 = t1923*t1927
      t1929 = t1925*t1894
      t1930 = t1896*e4
      t1934 = t1901*t1880-t1885
      t1935 = t1934*e5
      t1936 = t1935*t1891
      t1937 = t1923*t1936
      t1939 = e5*t1896
      t1942 = -2*t1928-2*t1923*t1929*t1930+2.d0*t1937
     +        -4.d0*t1923*t1934*t1888*t1939
      t1943 = t1-t0
      t1946 = 1/e1
      dxd14 = (-2*t1893-2*t1887*t1898+2.d0*t1908-4.d0*t1911
     +         -0.15d1*t1942*t1943)*t1946
      t1947 = t1887*t1907
      t1948 = t1897*e5
      t1950 = t1906*t1892
      t1951 = e4*t1889
      t1952 = t1951*t1896
      t1953 = t1906*t1952
      t1956 = t1923*t1925*e5*t1891
      t1959 = t1934*e4
      t1961 = t1923*t1959*t1891
      t1962 = t1889*t1896
      t1965 = 2*t1956-2*t1923*t1929*t1939+2.d0*t1961
     +      -4.d0*t1923*t1959*t1962
      dxd15 = (2*t1947-2*t1887*t1948+2.d0*t1950-4.d0*t1953
     +       -0.15d1*t1965*t1943)*t1946
      t1968 = t1887*t1910
      t1969 = 1.d0+t1888-t1889
      t1970 = t1969*t1896
      t1971 = t1970*e4
      t1973 = t1925*t1888
      t1976 = t1934*t1969
      t1979 = 2.d0*t1956-4.d0*t1923*t1973*t1939+2*t1961
     +       -2*t1923*t1976*t1930
      dyd14 = (2.d0*t1947-4.d0*t1968+2*t1950-2*t1906*t1971
     +       -0.15d1*t1979*t1943)*t1946
      t1982 = t1887*t1952
      t1983 = t1970*e5
      t1989 = 2.d0*t1928-4.d0*t1923*t1926*t1962-2*t1937
     +        -2*t1923*t1976*t1939
      dyd15 = (2.d0*t1893-4.d0*t1982-2*t1908-2*t1906*t1983
     +        -0.15d1*t1989*t1943)*t1946
      t1992 = t1887*t1891
      t1993 = t1888*t1896
      t1994 = t1887*t1993
      t1995 = t1939*e4
      t1996 = t1906*t1995
      t2004 = -2.d0*t1914*t1922*t1925*t1891+4.d0*t1923*t1973*t1896
     +        -4.d0*t1923*t1935*t1930
      dzd14 = (-2.d0*t1992+4.d0*t1994-4.d0*t1996
     +        -0.15d1*t2004*t1943)*t1946
      t2007 = t1887*t1995
      t2008 = t1906*t1891
      t2009 = t1906*t1962
      t2018 = 4.d0*t1923*t1926*t1939+2.d0*t1914*t1922*t1934*t1891
     +        -4.d0*t1923*t1934*t1889*t1896
      dzd15 = (4.d0*t2007+2.d0*t2008-4.d0*t2009
     +        -0.15d1*t2018*t1943)*t1946
      t2021 = 1/t1915
      t2025 = e6+t1912*t1913*t2021*t1943-el
      t2026 = t1876**2
      t2028 = 1/t2026/t1876
      t2030 = 1.d0-t1877
      t2031 = 1/t2030
      t2032 = t1872*t2028*t2031
      t2033 = t1877+t2032
      t2035 = e1*t1880
      t2036 = t1877*e2
      t2037 = t2036-t1883
      t2038 = t2037*t1922
      t2039 = t2035*t2038
      t2040 = t2025*t2033+t2039
      t2041 = e1*t2040
      t2042 = t1892*t2041
      t2045 = t2025*e2
      t2047 = e3*t2028*t2031
      t2048 = t2045*t2047
      t2049 = t1877*e3
      t2050 = t2049-t1880
      t2051 = t2050*t1922
      t2052 = t2035*t2051
      t2053 = t2048-1.d0+t2052
      t2054 = e1*t2053
      t2055 = t1907*t2054
      t2056 = t1896*e1
      t2057 = t2056*t2053
      dxd24 = 2*t2042+2*t1897*t2041*e4+2.d0*t2055-4.d0*t1909*t2057
      t2059 = t1907*t2041
      t2060 = t2041*e5
      t2062 = t1892*t2054
      dxd25 = -2*t2059+2*t1897*t2060+2.d0*t2062-4.d0*t1951*t2057
      t2064 = t2056*t2040
      t2066 = t2054*e4
      dyd24 = -2.d0*t2059+4.d0*t1909*t2064+2*t2062-2*t1970*t2066
      dyd25 = -2.d0*t2042+4.d0*t1951*t2064-2*t2055-2*t1970*t2054*e5
      t2071 = t1891*e1
      dzd24 = 2.d0*t2071*t2040-4.d0*t1993*t2041-4.d0*t1939*t2066
      dzd25 = -4.d0*t1930*t2060+2.d0*t2071*t2053-4.d0*t1962*t2054
      t2078 = e1*t1883
      t2079 = t2078*t2038
      t2080 = t2048+1.d0-t2079
      t2081 = e1*t2080
      t2082 = t1892*t2081
      t2086 = t1873*t2028*t2031
      t2087 = t1877+t2086
      t2089 = t2078*t2051
      t2090 = t2025*t2087-t2089
      t2091 = e1*t2090
      t2092 = t1907*t2091
      t2093 = t2056*t2090
      dxd34 = 2*t2082+2*t1897*t2081*e4+2.d0*t2092-4.d0*t1909*t2093
      t2095 = t1907*t2081
      t2096 = t2081*e5
      t2098 = t1892*t2091
      dxd35 = -2*t2095+2*t1897*t2096+2.d0*t2098-4.d0*t1951*t2093
      t2100 = t2056*t2080
      t2102 = t2091*e4
      dyd34 = -2.d0*t2095+4.d0*t1909*t2100+2*t2098-2*t1970*t2102
      dyd35 = -2.d0*t2082+4.d0*t1951*t2100-2*t2092-2*t1970*t2091*e5
      dzd34 = 2.d0*t2071*t2080-4.d0*t1993*t2081-4.d0*t1939*t2102
      dzd35 = -4.d0*t1930*t2096+2.d0*t2071*t2090-4.d0*t1962*t2091
      t2118 = t1894*t1891
      t2119 = t1906*t2118
      t2120 = e4*e5
      t2121 = t2120*t1891
      t2122 = t1887*t2121
      t2126 = (e5*(t2119-2.d0*t2122)-2.d0*t1893)*t1896
      dxd44 = 2.d0*(e5*(-2*t1950-2*t1906*t1898-2.d0*t1947+4.d0*t1968)
     +      -2.d0*t1992+4.d0*t1994)*t1891-4.d0*t2126*e4
      t2128 = t1906*t1948
      dxd45 = 2.d0*(t2119-2.d0*t2122+e5*(2*t1908-2*t2128-2.d0*t1893
     +      +4.d0*t1982)+4.d0*t2007)*t1891-4.d0*t2126*e5
      t2139 = t1906*t2121
      t2140 = t1969*t1891
      t2141 = t1887*t2140
      t2145 = (e5*(2.d0*t2139-t2141)+2.d0*t1947)*t1896
      dyd44 = 2.d0*(e5*(2.d0*t1908-4.d0*t1911-2*t1893+2*t1887*t1971)
     +       -4.d0*t2007)*t1891-4.d0*t2145*e4
      t2147 = t1887*t1983
      t2150 = t1887*t1962
      dyd45 = 2.d0*(2.d0*t2139-t2141+e5*(2.d0*t1950
     +       -4.d0*t1953+2*t1947+2*t2147)+2.d0*t1992
     +       -4.d0*t2150)*t1891-4.d0*t2145*e5
      t2157 = 1.d0-t1888-t1889
      t2158 = t2157*t1896
      t2165 = t2157*t1891
      t2168 = (e5*(-2.d0*t1950-2.d0*t1947)-t1887*t2165)*t1896
      dzd44 = 2.d0*(e5*(-2.d0*t2008+4.d0*t1906*t1993+4.d0*t2007)
     +       +2*t1893+2*t1887*t2158*e4)*t1891-4.d0*t2168*e4
      t2172 = t2158*e5
      dzd45 = 2.d0*(-2.d0*t1950+e5*(4.d0*t1996-2.d0*t1992+4.d0*t2150)
     +       +2*t1887*t2172)*t1891-4.d0*t2168*e5
      dxd55 = 2.d0*(e4*(-2*t1908+2*t2128+2.d0*t1893-4.d0*t1982)
     +       -4.d0*t1996)*t1891-4.d0*(e4*(-t2119+2.d0*t2122)
     +       +2.d0*t1950)*t1896*e5
      dyd55 = 2.d0*(e4*(-2.d0*t1950+4.d0*t1953-2*t1947-2*t2147)
     +       -2.d0*t2008+4.d0*t2009)*t1891-4.d0*(e4*(-2.d0*t2139+t2141)
     +       -2.d0*t1908)*t1896*e5
      dzd55 = 2.d0*(e4*(-4.d0*t1996+2.d0*t1992-4.d0*t2150)-2*t1908
     +       -2*t1906*t2172)*t1891-4.d0*(e4*(2.d0*t1950+2.d0*t1947)
     +       +t1906*t2165)*t1896*e5
      t2206 = 1/t1912
      t2208 = t1913*e1
      dxd46 = t1942*t2206*t2208
      dxd56 = t1965*t2206*t2208
      dyd46 = t1979*t2206*t2208
      dyd56 = t1989*t2206*t2208
      dzd46 = t2004*t2206*t2208
      dzd56 = t2018*t2206*t2208
      t2214 = 1/t2026
      t2215 = t1872*e2
      t2217 = 1/t1875
      t2219 = -t2214*t2215*t2217-2*t2036
      t2221 = t2049*t1883
      t2222 = t1872*e3
      t2225 = t2222*t2214*t1883*t2217
      t2226 = t2219*t1880+t2221+0.1d1*t2225
      t2227 = e1*t2226
      t2229 = t2214*t1873
      t2232 = t2229*t2217*e2*t1883
      t2233 = t2049*t1880
      t2235 = t2214*t1880*t2217
      t2236 = t2222*t2235
      t2237 = -t2232+t2233+0.1d1*t2236-1
      t2238 = e1*t2237
      t2240 = t1920**2
      t2242 = t1921/t2240
      t2243 = t1914*t2242
      t2244 = t1915*t1886
      t2246 = t1915*t1905
      t2248 = 2*t2244*t2226+2*t2246*t2237
      t2249 = t1891*t2248
      t2253 = t2233+0.1d1*t2236-t2219*t1883
      t2258 = t1914*t2242*t1934
      t2259 = t2120*t2249
      t2262 = e2*t1873*t2235
      t2263 = -t2262-t2221-t2225
      t2271 = e1*t1925
      t2273 = e1*t1934
      t2277 = 2*t2244*t1925+2*t2246*t1934
      t2278 = t1891*t2277
      t2281 = -t1885-t1881
      t2285 = t2120*t2278
      t2287 = -t1902-t1904
      t2291 = -t2243*t1929*t2278/2+t1923*t2281*t1894*t1891
     +       -t2258*t2285+2.d0*t1923*t2287*e4*t1907
      t2293 = t2271*t2118+2.d0*t2273*t2121-0.15d1*t2291*t1943
      dxd12 = (t2227*t2118+2.d0*t2238*t2121-0.15d1*
     +      (-t2243*t1929*t2249/2+t1923*t2253*t1894*t1891
     +      -t2258*t2259+2.d0*t1923*t2263*e4*t1907)*t1943)*t1946
     +      -t2293*t1880*t1922
      t2299 = t1914*t2242*t1925
      t2301 = t2253*e4
      t2316 = t2281*e4
      t2324 = -t2299*t2285+2.d0*t1923*t2316*t1907-t2243*t1976*t2278/2
     +      +t1923*t2287*t1969*t1891
      t2326 = 2.d0*t2271*t2121+t2273*t2140-0.15d1*t2324*t1943
      dyd12 = (2.d0*t2227*t2121+t2238*t2140-0.15d1*(-t2299*t2259+
     +        2.d0*t1923*t2301*t1907-t2243*t1976*t2249/2
     +       +t1923*t2263*t1969*t1891)*t1943)*t1946-t2326*t1880*t1922
      t2355 = 0.1d1*t2243*t1926*t2278-2.d0*t1923*t2316*t1891
     +        -t2243*t1935*t2278+2.d0*t1923*t2287*e5*t1891
      t2357 = -2.d0*t2271*t1892+2.d0*t2273*t1907-0.15d1*t2355*t1943
      dzd12 = (-2.d0*t2227*t1892+2.d0*t2238*t1907-0.15d1*
     +        (0.1d1*t2243*t1926*t2249-2.d0*t1923*t2301*t1891
     +        -t2243*t1935*t2249+2.d0*t1923*t2263*e5*t1891)*t1943)
     +        *t1946-t2357*t1880*t1922
      t2360 = t2036*t1883
      t2361 = -t2236+t2360+0.1d1*t2232-1
      t2362 = e1*t2361
      t2364 = t1873*e3
      t2367 = -t2214*t2364*t2217-2*t2049
      t2369 = t2036*t1880
      t2370 = t2367*t1883+t2369+0.1d1*t2262
      t2371 = e1*t2370
      t2375 = 2*t2244*t2361+2*t2246*t2370
      t2376 = t1891*t2375
      t2379 = t2369+0.1d1*t2262+0.1d1*t2225
      t2383 = t2120*t2376
      t2386 = t2367*t1880-t2360-t2232
      dxd13 = (t2362*t2118+2.d0*t2371*t2121-0.15d1*(-t2243*t1929*t2376/2
     +   +t1923*t2379*t1894*t1891-t2258*t2383+2.d0*t1923*t2386*e4*t1907)
     +   *t1943)*t1946+t2293*t1883*t1922
      t2399 = t2379*e4
      dyd13 = (2.d0*t2362*t2121+t2371*t2140-0.15d1*(-t2299*t2383
     +   +2.d0*t1923*t2399*t1907-t2243*t1976*t2376/2+
     +    t1923*t2386*t1969*t1891)*t1943)*t1946+t2326*t1883*t1922
      dzd13 = (-2.d0*t2362*t1892+2.d0*t2371*t1907-0.15d1*
     +     (0.1d1*t2243*t1926*t2376-2.d0*t1923*t2399*t1891
     +     -t2243*t1935*t2376+2.d0*t1923*t2386*e5*t1891)*t1943)
     +      *t1946+t2357*t1883*t1922
      dxd16 = t2293*t1922
      dyd16 = t2326*t1922
      dzd16 = t2357*t1922
      t2430 = t2214*t2217
      t2431 = t2430*e3
      t2432 = t2026**2
      t2433 = 1/t2432
      t2435 = t2031*t2217
      t2439 = 1/t2432/t1876
      t2441 = t2030**2
      t2442 = 1/t2441
      t2443 = t2442*t2217
      t2451 = e3*t2217*e2*t1922
      t2452 = t2035*t2214*t2451
      t2453 = t2037*t2242
      t2454 = t2453*t2375
      t2456 = t2025*(0.1d1*t2431+0.3d1*t1872*t2433*t2435*e3
     +    +0.1d1*t1872*t2439*t2443*e3)+0.1d1*t2452-t2035*t2454/2
      t2457 = e1*t2456
      t2459 = t2028*t2031
      t2460 = t2045*t2459
      t2461 = t2045*t1873
      t2463 = t2433*t2031*t2217
      t2464 = t2461*t2463
      t2466 = t2439*t2442*t2217
      t2467 = t2461*t2466
      t2470 = (0.1d1*t2229*t2217+t1877)*t1922
      t2472 = t2050*t2242
      t2473 = t2472*t2375
      t2475 = t2460+0.3d1*t2464+0.1d1*t2467+t2035*t2470-t2035*t2473/2
      t2478 = t1880**2
      t2481 = t2453*t2277
      t2483 = -t1877-t2032-t2079-e1*t2478*t1922-t2035*t2481/2
      t2484 = e1*t2483
      t2486 = t1882*t2459
      t2487 = t1883*t1922
      t2488 = t2035*t2487
      t2489 = t2472*t2277
      t2491 = -t2486-t2089+t2488-t2035*t2489/2
      t2495 = (-t2118*t2484+2.d0*t2120*t2071*t2491)*e1
      dxd23 = -t2118*t2457+2.d0*t2120*t2071*t2475+t2495*t2487
      t2499 = e1*t2475
      t2503 = e1*t2491
      t2506 = (-2.d0*t2120*t2071*t2483+t2140*t2503)*e1
      dyd23 = -2.d0*t2120*t2071*t2456+t2140*t2499+t2506*t2487
      t2513 = (2.d0*t1892*t2484+2.d0*t1907*t2503)*e1
      dzd23 = 2.d0*t1892*t2457+2.d0*t1907*t2499+t2513*t2487
      t2515 = e1*t2033
      t2518 = t2120*t2071*t2486
      dxd26 = -t2118*t2515+2.d0*t2518+t2495*t1922
      dyd26 = -2.d0*t2120*t2071*t2033+t2140*e1*t2486+t2506*t1922
      dzd26 = 2.d0*t1892*t2515+2.d0*t1907*e1*t2486+t2513*t1922
      t2534 = -t2486-t2039+t2488+t2078*t2481/2
      t2535 = e1*t2534
      t2537 = t1883**2
      t2541 = -t1877-t2086-t2052-e1*t2537*t1922+t2078*t2489/2
      t2545 = (-t2118*t2535+2.d0*t2120*t2071*t2541)*e1
      dxd36 = -t2118*e1*t2486+2.d0*t2120*t2071*t2087+t2545*t1922
      t2547 = e1*t2087
      t2551 = e1*t2541
      t2554 = (-2.d0*t2120*t2071*t2534+t2140*t2551)*e1
      dyd36 = -2.d0*t2518+t2140*t2547+t2554*t1922
      t2562 = (2.d0*t1892*t2535+2.d0*t1907*t2551)*e1
      dzd36 = 2.d0*t1892*e1*t2486+2.d0*t1907*t2547+t2562*t1922
      t2570 = t1912/t1913*t1922
      t2571 = t1929*t1891
      t2576 = t1891*(2*e1*t1916+2*e1*t1918)
      t2579 = t1959*t1907
      t2581 = t2120*t2576
      t2596 = t1913/t1915/e1
      t2599 = t1912*t1943*t1922
      dxd11 = (t1886*t1894*t1891+2.d0*t1905*e4*t1907-0.15d1*(t2570*t2571
     +     /2-t2243*t1929*t2576/2+0.1d1*t2570*t2579
     +     -t2258*t2581)*t1943)*t1946
     +    -(t1887*t2118+2.d0*t2139-0.15d1*(t1923*t2571
     +     +2.d0*t1923*t2579)*t1943)*t2021-0.15d1*t2293*t2596*t2599
      t2601 = t1886*e4
      t2605 = t1926*t1907
      t2608 = t1976*t1891
      dyd11 = (2.d0*t2601*t1907+t1905*t1969*t1891-0.15d1*
     +   (0.1d1*t2570*t2605-t2299*t2581+t2570*t2608/2
     +   -t2243*t1976*t2576/2)*t1943)*t1946-(2.d0*t2122
     +   +t1906*t2140-0.15d1*(2.d0*t1923*t2605+t1923*t2608)*t1943)
     +   *t2021-0.15d1*t2326*t2596*t2599
      dzd11 = (-2.d0*t2601*t1891+2.d0*t1905*e5*t1891-0.15d1*
     +   (-t2570*t1927+0.1d1*t2243*t1926*t2576+0.1d1*t2570*t1936
     +    -t2243*t1935*t2576)*t1943)*t1946-(-2.d0*t1893+2.d0*t1908
     +    -0.15d1*(-2.d0*t1928+2.d0*t1937)*t1943)*t2021
     +    -0.15d1*t2357*t2596*t2599
      t2660 = t2025*(0.1d1*t2430*e2+2*e2*t2028*t2031
     +    +0.3d1*t2215*t2433*t2435+0.1d1*t2215*t2439*t2443)
     +    +t2035*(0.1d1*t2214*t1872*t2217+t1877)*t1922
     +    -t2035*t2453*t2248/2
      t2661 = e1*t2660
      t2666 = t2025*t1872*e3
      t2671 = t2025*e3*t2459+0.3d1*t2666*t2463+0.1d1*t2666*t2466
     +    +0.1d1*t2452-t2035*t2472*t2248/2
      t2674 = t1880*t1922
      dxd22 = -t2118*t2661+2.d0*t2120*t2071*t2671-t2495*t2674
      t2678 = e1*t2671
      dyd22 = -2.d0*t2120*t2071*t2660+t2140*t2678-t2506*t2674
      dzd22 = 2.d0*t1892*t2661+2.d0*t1907*t2678-t2513*t2674
      t2687 = t2460+0.3d1*t2464+0.1d1*t2467-t2078*t2214*t2451
     +    +t2078*t2454/2
      t2688 = e1*t2687
      t2698 = t2025*(0.1d1*t2431+2*t2047+0.3d1*t2364*t2433*t2435
     +    +0.1d1*t2364*t2439*t2443)-t2078*t2470+t2078*t2473/2
      dxd33 = -t2118*t2688+2.d0*t2120*t2071*t2698+t2545*t2487
      t2704 = e1*t2698
      dyd33 = -2.d0*t2120*t2071*t2687+t2140*t2704+t2554*t2487
      dzd33 = 2.d0*t1892*t2688+2.d0*t1907*t2704+t2562*t2487
      t2712 = t1913*t1915*t1922
      dxd66 = t2291*t2206*t2712
      dyd66 = t2324*t2206*t2712
      dzd66 = t2355*t2206*t2712
c  Give names with indexes
      ddxde(1,1,4)=dxd14
      ddxde(1,1,5)=dxd15
      ddxde(1,2,4)=dxd24
      ddxde(1,2,5)=dxd25
      ddxde(1,3,4)=dxd34
      ddxde(1,3,5)=dxd35
      ddxde(1,4,4)=dxd44
      ddxde(1,4,5)=dxd45
      ddxde(1,4,6)=dxd46
      ddxde(1,5,5)=dxd55
      ddxde(1,5,6)=dxd56
      ddxde(2,1,4)=dyd14
      ddxde(2,1,5)=dyd15
      ddxde(2,2,4)=dyd24
      ddxde(2,2,5)=dyd25
      ddxde(2,3,4)=dyd34
      ddxde(2,3,5)=dyd35
      ddxde(2,4,4)=dyd44
      ddxde(2,4,5)=dyd45
      ddxde(2,4,6)=dyd46
      ddxde(2,5,5)=dyd55
      ddxde(2,5,6)=dyd56
      ddxde(3,1,4)=dzd14
      ddxde(3,1,5)=dzd15
      ddxde(3,2,4)=dzd24
      ddxde(3,2,5)=dzd25
      ddxde(3,3,4)=dzd34
      ddxde(3,3,5)=dzd35
      ddxde(3,4,4)=dzd44
      ddxde(3,4,5)=dzd45
      ddxde(3,4,6)=dzd46
      ddxde(3,5,5)=dzd55
      ddxde(3,5,6)=dzd56
c  Second derivatives with respect to a,h,k,lambda
      ddxde(1,1,1)=dxd11
      ddxde(1,1,2)=dxd12
      ddxde(1,1,3)=dxd13
      ddxde(1,1,6)=dxd16
      ddxde(1,2,2)=dxd22
      ddxde(1,2,3)=dxd23
      ddxde(1,2,6)=dxd26
      ddxde(1,3,3)=dxd33
      ddxde(1,3,6)=dxd36
      ddxde(1,6,6)=dxd66
      ddxde(2,1,1)=dyd11
      ddxde(2,1,2)=dyd12
      ddxde(2,1,3)=dyd13
      ddxde(2,1,6)=dyd16
      ddxde(2,2,2)=dyd22
      ddxde(2,2,3)=dyd23
      ddxde(2,2,6)=dyd26
      ddxde(2,3,3)=dyd33
      ddxde(2,3,6)=dyd36
      ddxde(2,6,6)=dyd66
      ddxde(3,1,1)=dzd11
      ddxde(3,1,2)=dzd12
      ddxde(3,1,3)=dzd13
      ddxde(3,1,6)=dzd16
      ddxde(3,2,2)=dzd22
      ddxde(3,2,3)=dzd23
      ddxde(3,2,6)=dzd26
      ddxde(3,3,3)=dzd33
      ddxde(3,3,6)=dzd36
      ddxde(3,6,6)=dzd66
c  Schwartz
      do 88 i=1,3
        do 88 k=1,6
        do 88 j=k+1,6
 88       ddxde(i,j,k)=ddxde(i,k,j)
      return
      end






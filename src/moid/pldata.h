c planet data

c ap = semimajor axis of the planet, ky = Gauss constant 
c gm = (masse pianeti)/(massa Sole)

       INTEGER nplax,nplax2,npl,inpl,ioupl,ndum
       parameter (nplax=10)
       parameter (nplax2=20)
       DOUBLE PRECISION ap(nplax),gm(nplax)
       DOUBLE PRECISION ky,bigg 

       COMMON/pldata/gm,ap,ky,ndum,npl,inpl,ioupl,bigg



      real function airy(x)
c  calculate airy function (ai(-x))
      pi=4*atan(1.)
c
      x3=x*x*x
c
      if(x-2.2990176) 10,10,20
   10 fx=1.-x3*(1-x3*(1-x3*(1-x3/132)/72)/30)/6
      gx=x*(1-x3*(1-x3*(1-x3*(1-x3/156)/90)/42)/12)
      airy=0.355028053887817*fx+0.258819403792807*gx
      return
c
   20 x2=0.125/x3
      x32=0.66666666667*sqrt(x3)
      tht=pi/4-x32*(1-x2*(1.25-x2*(11.510417-x2*647.07031)))
      am=(1-1.25*x2*(1-28.875*x2*(1-92.08333333*x2)))
     *  /(pi*sqrt(x))
      if(am.le.0) am=-am
      am=sqrt(am)
      airy=am*cos(tht)
      return
      end
      real function dairy*8(x)
      implicit real*8 (a-h,o-z)
c  calculate airy function (ai(-x))
      pi=4*datan(1.d0)
c
      x3=x*x*x
c
      if(x-2.2990176d0) 10,10,20
   10 fx=1.d0-x3*(1-x3*(1-x3*(1-x3/132)/72)/30)/6
      gx=x*(1-x3*(1-x3*(1-x3*(1-x3/156)/90)/42)/12)
      airy=0.355028053887817d0*fx+0.258819403792807d0*gx
      return
c
   20 x2=0.125d0/x3
      x32=0.66666666667d0*dsqrt(x3)
      tht=pi/4-x32*(1-x2*(1.25d0-x2*(11.510417d0-x2*647.07031d0)))
      am=(1-1.25d0*x2*(1-28.875d0*x2*(1-92.08333333d0*x2)))
     *  /(pi*dsqrt(x))
      if(am.le.0) am=-am
      am=dsqrt(am)
      dairy=am*dcos(tht)
      return
      end

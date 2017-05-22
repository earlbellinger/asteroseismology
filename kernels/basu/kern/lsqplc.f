      subroutine lsqplc(x,nn,ix,kk,ac,iac,idiag)
c
c  calculates expansion coefficients ac which determine coefficients
c  least-squares polynomial fit of degree kk. thus the coefficients
c  in the fit
c
c       y = a(1) + a(2)*x + ... + a(kk+1)*x**kk
c
c  to data points (x(1,i), y(1,i)), i = 1, ..., nn,
c  are given by
c
c      a(k)=sum(ac(k,i)*y(1,i)), k = 1, ..., kk + 1
c            i
c
c  ix and iac are first dimensions of x and ac.
c
c  ......................................................................
c
      dimension x(ix,1),ac(iac,1),w(1000),aa(20,20),w1(1000)
c
      do 10 n=1,nn
   10 w(n)=1
c
c  to avoid over and underflow problems rescale x to
c  be between -1 and 1.
c
      ixmax=isamax(nn,x,ix)
      xscl=x(1,ixmax)
      xscli=1./xscl
c
      do 12 n=1,nn
   12 w1(n)=xscli*x(1,n)
c
      kk1=kk+1
      kk2=kk+kk1
c
      do 30 k=1,kk2
      sk=ssum(nn,w,1)
      if(k.le.kk1) call scopy(nn,w,1,ac(k,1),iac)
      i1=max0(1,k-kk1+1)
      i2=min0(k,kk1)
      do 15 i=i1,i2
   15 aa(i,k+1-i)=sk
      if(k.eq.kk2) go to 30
      do 20 n=1,nn
   20 w(n)=w1(n)*w(n)
   30 continue
c
c  test for output of equations
c
      if(idiag.lt.2) go to 36
c
      write(6,30090)
      do 35 i=1,kk1
   35 write(6,30091) (aa(i,j),j=1,kk1),(ac(i,j),j=1,nn)
30090 format(///' aa:'/)
30091 format(1p10e13.5)
c
   36 call leq(aa,ac,kk1,nn,20,iac,det)
c
c
c  rescale coefficients
c
      scl=1
      do 38 k=1,kk1
      call sscal(nn,scl,ac(k,1),iac)
   38 scl=scl*xscli
c
      return
      end

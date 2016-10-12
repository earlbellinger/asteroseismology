      subroutine derivk(x,y,dydx,ii,nn,iy,idy,karg)  
c 
c      2k-th order differentiation routine
c      ***********************************
c 
c  sets derivative of y(i,n) with respect to x into dydx(i,n),i=1,ii,n=1  
c  id and idy are first dimensions of y and dydx respectively.
c  if x is found to be non-monotonic ii is set to 0.  
c 
c  if derivk is called with karg .le. 0, it is assumed that a table
c  of differentiation coefficients has been set up.
c
c  derivk uses as work space 
c      wrk((2k+1)*(2k+1)) 
c  the size of wrk set in derivk is 100, which is sufficient, say, for   
c  8-th order differentiation in 20 dependent variables. if more 
c  work space is needed it must be set in common/wrklir/ 
c 
c  note that this version of derivk also uses storage of size
c
c         dc((2*k+1)*nn)
c
c  in common /cderst/ to store differentiation coefficients. the
c  default size of 500 probably in general has to be increased in
c  the calling programme.
c
c  original version: 18/3/1986
c
      dimension x(1),y(iy,1),dydx(idy,1)
      common/wrklir/ wrk(100)
      common/cderst/ dc(500)
c
      save /cderst/,k,k2
c
c  test for setting up coefficients
c
      if(karg.le.0) go to 60
c
      k=karg
c
c  check order   
c
      k2=2*k+1   
      if(nn.lt.k2) then
        write(6,100) k2,nn  
        k=(nn-1)/2 
        k2=2*k+1   
      end if
c
c  first dimension of wrk, storage information   
c
      iw=k2  
      mst=iw*(iw-1)+2
c
c  total range   
c
      diff=x(nn)-x(1)
      if(diff.eq.0) go to 95
c
c  start loop for setting coefficients
c
   20 jdc=0
c
      do 50 n=1,nn   
      xtr=x(n)   
c
c  determine i   
c
      i=min0(max0(0,n-1-k),nn-k2)
c  length of interval
      dlx=x(i+k2)-x(i+1) 
      if(diff*dlx.le.0) go to 97
c
c  set equations 
c
   32 dlxi=1./dlx
c
      j1=jdc
      m=0
      do 40 l=1,k2   
      l1=i+l 
      j1=j1+1
      xl=(x(l1)-xtr)*dlxi
      aa=1.  
      m=m+1
      wrk(m)=aa  
      do 35 ir=2,k2  
      m=m+1
      aa=xl*aa   
   35 wrk(m)=aa  
   40 dc(j1)=0
c
      dc(jdc+2)=dlxi
c
c  solve equations   
c
      call leq(wrk,dc(jdc+1),k2,1,iw,1,det)   
c
      if(det.eq.0) go to 98
c
   50 jdc=jdc+k2
c
c  end setting coefficients
c
      write(6,110) jdc
c
c            *****************************************
c
c  set derivatives
c
   60 jdc=0
c
      do 70 n=1,nn
c
c  determine i   
c
      i=min0(max0(0,n-1-k),nn-k2)
c
      do 65 is=1,ii
      j1=jdc
      sum=0
      do 62 l=1,k2
      j1=j1+1
   62 sum=sum+dc(j1)*y(is,i+l)
c
   65 dydx(is,n)=sum
c
   70 jdc=jdc+k2
c
      return 
c
c  diagnostics   
c
   95 write(6,120) x(1)  
      ii=0   
      return 
   97 i1=i+1 
      ik2=i+k2   
      write(6,130) x(1),x(nn),i1,x(i1),ik2,x(ik2)
      ii=0   
      return 
   98 i1=i+1 
      ik2=i+k2   
      write(6,140) (l,x(l),l=i1,ik2) 
      ii=0   
      return 
  100 format(1x,10('*'),' 2k+1 =',i4,' is greater than nn =',i4, 
     *  ' in derivk.'/11x,'k has been reset to (nn-1)/2')
  110 format(//' storage needed in common/cderst/ in derivk:',
     *  i6)
  120 format(//1x,10('*'),' range is zero in derivk. x(1)=', 
     *  1pe15.5,' = x(nn)'//)
  130 format(//1x,10('*'),' independent variable is not monotonic',  
     *  ' in derivk'//' x(1)=',1pe15.5,' x(nn)=',e15.5,  
     *  2(' x(',i4,')=',e15.5)//)
  140 format(//1x,10('*'),' points coincide in derivk'/  
     *  (5(' x(',i4,')=',1pe13.4)))  
      end

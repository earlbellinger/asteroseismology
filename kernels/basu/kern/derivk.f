      subroutine derivk(x,y,dydx,ii,nn,iy,idy,k)  
c 
c      2k-th order differentiation routine
c      ***********************************
c 
c  sets derivative of y(i,n) with respect to x into dydx(i,n),i=1,ii,n=1  
c  id and idy are first dimensions of y and dydx respectively.
c  if x is found to be non-monotonic ii is set to 0.  
c 
c  derivk uses as work space 
c      wrk((2k+1)*(2k+1+ii)) 
c  the size of wrk set in derivk is 300, which is sufficient, say, for   
c  8-th order differentiation in 20 dependent variables. if more 
c  work space is needed it must be set in common/wrklir/ 
c
      dimension x(1),y(iy,1),dydx(idy,1)
      common/wrklir/ wrk(300)
    5 iwr=9  
c  check order   
      k2=2*k+1   
      if(nn-k2) 10,15,15 
   10  write(iwr,100) k2,nn  
      k=(nn-1)/2 
      k2=2*k+1   
c  first dimension of wrk, storage information   
   15 iw=k2  
      istrh=iw*iw+1  
      mst=iw*(iw-1)+2
c  total range   
      diff=x(nn)-x(1)
      if(diff) 20,95,20  
   20 do 70 n=1,nn   
      xtr=x(n)   
c  determine i   
   30 i=min0(max0(0,n-1-k),nn-k2)
c  length of interval
      dlx=x(i+k2)-x(i+1) 
      if(diff*dlx) 97,97,32  
c  set equations 
   32 do 40 l=1,k2   
      l1=i+l 
      xl=(x(l1)-xtr)/dlx 
      aa=1.  
      wrk(l)=aa  
      m=l
      do 35 ir=2,k2  
      m=m+iw 
      aa=xl*aa   
   35 wrk(m)=aa  
      do 40 is=1,ii  
      m=m+iw 
   40 wrk(m)=y(is,l1)
c  solve equations   
      call leq(wrk,wrk(istrh),k2,ii,iw,iw,det)   
      if(det) 50,98,50   
c  set derivatives   
   50 m=mst  
      do 63 is=1,ii  
      m=m+iw 
   63 dydx(is,n)=wrk(m)/dlx  
   70 continue   
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
  120 format(//1x,10('*'),' range is zero in derivk. x(1)=', 
     *  1pe15.5,' = x(nn)'//)
  130 format(//1x,10('*'),' independent variable is not monotonic',  
     *  ' in derivk'//' x(1)=',1pe15.5,' x(nn)=',e15.5,  
     *  2(' x(',i4,')=',e15.5)//)
  140 format(//1x,10('*'),' points coincide in derivk'/  
     *  (5(' x(',i4,')=',1pe13.4)))  
      end

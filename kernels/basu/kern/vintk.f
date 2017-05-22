      subroutine vintk(x,a,y,ii,nn,ia,id,k) 
c   
c      2k-1-th order integration routine
c      *********************************
c   
c  sets integral of a(i,n) with respect to x into y(i,n),i=1,ii,n=1,nn. 
c  the values of y(i,1), i=1,ii, must be supplied.  
c  ia and id are first dimensions of a and y respectively.  
c  if x is found to be non-monotonic ii is set to 0.
c   
c  vintk uses as work space 
c      wrk(2k*(2k+ii))  
c  the size of wrk set in vintk is 300, which is sufficient, say, for   
c  9-th order integration in 20 dependent variables. if more
c  work space is needed it must be set in common/wrklir/
c   
      dimension x(1),a(ia,1),y(id,1),wrk(300),  
     *  di(50)  
      common/wrklir/ wrk
      call store(y(1,1),di(1),ii)   
    5 iwr=6 
c  check order  
      k2=2*k
      if(nn-k2) 10,15,15
   10  write(iwr,100) k2,nn 
      k=nn/2
      k2=2*k
c  first dimension of wrk, storage information  
   15 iw=k2 
      istrh=iw*iw+1 
      mst=iw*(iw-1) 
c  total range  
      diff=x(nn)-x(1)   
      if(diff) 20,95,20 
   20 do 70 n=2,nn  
      n1=n-1
      xtr=x(n1) 
c  determine i  
   30 i=min0(max0(0,n1-k),nn-k2)
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
   40 wrk(m)=a(is,l1)   
c  solve equations  
      call leq(wrk,wrk(istrh),k2,ii,iw,iw,det)  
      if(det) 50,98,50  
c  set integrals
   50 xl=(x(n)-xtr)/dlx 
      aa=dlx
      do 60 l=1,k2  
      aa=aa*xl  
      ab=aa/l   
      m=mst+l   
      do 60 is=1,ii 
      m=m+iw
   60 di(is)=di(is)+wrk(m)*ab   
      call store(di(1),y(1,n),ii)   
c   
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
  100 format(1x,10('*'),' 2k =',i4,' is greater than nn =',i4,  
     *  ' in vintk.'/11x,'k has been reset to nn/2')
  120 format(//1x,10('*'),' range is zero in vintk. x(1)=', 
     *  1pe15.5,' = x(nn)'//)   
  130 format(//1x,10('*'),' independent variable is not monotonic', 
     *  ' in vintk'//' x(1)=',1pe15.5,' x(nn)=',e15.5,  
     *  2(' x(',i4,')=',e15.5)//)   
  140 format(//1x,10('*'),' points coincide in vintk'/  
     *  (5(' x(',i4,')=',1pe13.4))) 
      end   

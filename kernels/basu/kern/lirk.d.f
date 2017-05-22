      subroutine lirk(x,xi,y,yi,ii,id,nt,il,inter,k) 
c
c      interpolation routine 
c      ********************* 
c
c  sets y(i) = yi(i,a), i=1,ii for a such that x=xi(a), using 2k-1-th
c   order interpolation. 
c  note that extrapolation is not allowed.   
c  --------------------------------------
c  xi(n),yi(i,n) must be supplied for n=1,nt, i=1,ii.
c  id is first dimension of yi.  
c  if il.le.1 scan to find the xi(n) which immediately bounds x 
c   starts at n=1.  
c  if il.gt.1 scan starts at the value of n from previous call of lirk. 
c  inter is set to 1 if interpolation is successful. if x is outside ran
c   of xi inter is set to 0. if xi is found to be non-monotonic inter   
c   is set to -1.   
c   
c  lirk uses as work space  
c      wrk(2k*(2k+ii))  
c  the sixe of wrk set in lirk is 300, which is sufficient, say, for
c  10-th order interpolation in 20 dependent variables. if more 
c  work space is needed it must be set in common/wrklir/
c   
c  Double precision version.
c  ++++++++++++++++++++++++
c
c  Dated 10/3/90.
c
c  Note: this double precision version of the routine has same name
c  as single precision version, but is distinguished by its file name
c
      implicit double precision(a-h,o-z)
      dimension xi(1),yi(id,1),y(1),wrk(300),yst(30)
      common/wrklir/ wrk
      save
c
    5 iwr=6 
      inter=1   
c  check order  
      k2=2*k
      if(nt-k2) 10,15,15
   10  write(iwr,100) k2,nt 
      k=nt/2
      k2=2*k
c  find xi(j),xi(j+1) bracketing x  
   15 if(il.gt.1) go to 17  
      diff=xi(nt)-xi(1) 
      j=1   
      iold=-1   
   17 if(diff) 18,95,19 
   18 if(x-xi(j)) 20,80,25  
   19 if(x-xi(j)) 25,80,20  
   20 if(nt-j) 90,90,21 
   21 j=j+1 
      go to 17  
   25 if(j-1) 90,90,26  
   26 j=j-1 
c  determine i  
   30 i=min0(max0(0,j-k),nt-k2) 
      if(i.eq.iold) go to 50
c  set  equations   
      k21=k2-1  
c  length of interval, translation  
      dlx=xi(i+k2)-xi(i+1)  
      if(diff*dlx) 97,97,32 
   32 xtr=x 
c  first dimension of wrk   
      iw=k2 
      do 40 l=1,k2  
      l1=i+l
      xl=(xi(l1)-xtr)/dlx   
      aa=1. 
      wrk(l)=aa 
      m=l   
      do 35 ir=2,k2 
      m=m+iw
      aa=xl*aa  
   35 wrk(m)=aa 
      do 40 is=1,ii 
      m=m+iw
   40 wrk(m)=yi(is,l1)  
c  solve equations  
      m=1+iw*iw 
      call leq(wrk,wrk(m),k2,ii,iw,iw,det)  
      if(det) 50,98,50  
c  set interpolated values  
   50 aa=1. 
      m=m-iw
      do 55 is=1,ii 
      m=m+iw
   55 yst(is)=wrk(m)
      if(x-xtr) 57,70,57
   57 xx=(x-xtr)/dlx
      m1=k2*iw  
      do 60 l=2,k2  
      m1=m1+1   
      m=m1  
      aa=aa*xx  
      do 60 is=1,ii 
      m=m+iw
   60 yst(is)=yst(is)+wrk(m)*aa 
   70 continue  
      call store(yst,y,ii)  
      return
c  x falls on a meshpoint   
   80 call store(yi(1,j),y,ii)  
      return
c  we do not allow extrapolation
   90  write(iwr,110) x,xi(1),xi(nt)
      inter=0   
      return
   95 write(6,120) xi(1)
      write(6,125) (i,(xi(i+j),j=0,min(4,nt-i)),i=1,nt,5)
      inter=-1  
      return
   97 i1=i+1
      ik2=i+k2  
      write(6,130) xi(1),xi(nt),i1,xi(i1),ik2,xi(ik2)   
      inter=-1  
      return
   98 i1=i+1
      ik2=i+k2  
      write(6,140) (l,xi(l),l=i1,ik2)   
      inter=-1  
      return
  100 format(1x,10('*'),' 2k =',i4,' is greater than nt =',i4,  
     *  ' in lirk.'/11x,'k has been reset to nt/2') 
  110 format(///1x,10('*'),' x =',1pe15.5,' is outside range (',
     *  2e14.5,') of xi in lirk'//) 
  120 format(//1x,10('*'),' range is zero in lirk. xi(1)=', 
     *  1pe15.5,' = xi(nt)'//)  
  125 format(//' xi:'/(i5,1p5e13.5))
  130 format(//1x,10('*'),' independent variable is not monotonic', 
     *  ' in lirk'//' xi(1)=',1pe15.5,' xi(nt)=',e15.5, 
     *  2(' xi(',i4,')=',e15.5)//)  
  140 format(//1x,10('*'),' points coincide in lirk'/   
     *  (5(' xi(',i4,')=',1pe13.4))//)  
      end   

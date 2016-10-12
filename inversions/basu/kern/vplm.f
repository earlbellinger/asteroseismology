      subroutine vplm(l,m,x,plm,dplm,nx,iplm,idplm,inorm)   
c   
c  calculates the associated legendre function Plm and its derivative   
c  for all elements in the array x.
c  Hence plm(1,n) and dplm(1,n) return the value and the derivative of
c  Plm at x(n).
c  abs(m) must be .le. l. if not, plm and dplm are set to 0. 
c   
c  nb: sign convention is as in abramowitz and stegun   
c  **************************************************   
c  for inorm = 1 normalization is as in abramowitz and stegun   
c  for inorm = 2 normalization is such that 
c      integral from -1 to 1 (plm**2)dx = 4/(1+delta(m,0))
c   
c  Uses recursion relations, starting from l=m  and going towards 
c  increasing l. Hence, for most efficient use the routine should
c  be called in that order (e.g., m = 0, l = 0, 1, 2, ..., 
c  m=1, l = 1, 2, ... etc); however, the routine may be called
c  with any l, m.
c   
c  x(n) is restricted to have abs(x) .le. 1 - epsx, where epsx  
c  is initialized to 1.e-6. in addition m is reset to -m, if   
c  m .lt. 0.
c   
      dimension x(1),plm(iplm,1),dplm(idplm,1)  
      common/wrkplm/ wrk(6,401)
c
      save
c
      data lp,mp,inormp,nxp /-1,0,-1,-1/
      data epsx /1.e-6/
c   
      iwrk=6
c   
c  restrict m to be .ge. 0  
c   
      m=iabs(m) 
c   
c  test that abs(m).le.l
c   
      if(m.le.l) go to 10   
c   
c  otherwise set plm and dplm to zero   
c   
      do 5 i=1,nx   
      plm(1,i)=0
    5 dplm(1,i)=0   
      return
c   
c  restrict x to have abs(x) .lt. 1 
c   
   10 xlim=1.-epsx  
      do 12 n=1,nx  
   12 x(n)=amin1(xlim,amax1(x(n),-xlim))
c   
c  test for same x as in previous call (wrk(5,n) contains previous x)   
c   
      if(nx.ne.nxp) go to 20
c   
      do 14 n=1,nx  
   14 wrk(6,n)=x(n)-wrk(5,n)
      idx=isamax(nx,wrk(6,1),iwrk)  
      if(abs(wrk(6,idx)).gt.1.e-8) go to 20 
c   
c  test for initialization  
c   
      if(inorm.ne.inormp.or.m.ne.mp.or.l.lt.lp) go to 20
c   
c  test for identical case  
c   
      if(l.gt.lp) go to 30  
c   
c  identical case. store previous plm and dplm  
c   
      call scopy(nx,wrk(iw3,1),iwrk,plm,iplm)   
      call scopy(nx,wrk(4,1),iwrk,dplm,idplm)   
      return
c   
c  initialize at l = m  
c   
   20 do 22 n=1,nx  
      wrk(4,n)=1-x(n)*x(n)  
      wrk(3,n)=sqrt(wrk(4,n))   
   22 wrk(3,n)=wrk(3,n)**m  
c   
c  set normalization factor, depending on inorm 
c   
      if(inorm.eq.2) go to 25   
      fct=1 
      if(m.le.1) go to 24   
      do 23 k=2,m   
   23 fct=fct*(2*k-1)   
   24 if(mod(m,2).eq.1) fct=-fct
      go to 28  
c   
   25 fct=2*m+1 
      if(m.le.1) go to 27   
      do 26 k=2,m   
   26 fct=fct*(2*k-1)/float(2*k)
   27 fct=sqrt(fct)
c  
c  multiply by normalization factor
c  
   28 call sscal(nx,fct,wrk(3,1),iwrk) 
c  
c  initialize iw3 and lp   
c  
      iw3=3
      lp=m 
c  
c  if l = m set derivative and go to finish
c  
      if(l.gt.m) go to 30  
      do 29 n=1,nx 
   29 wrk(4,n)=-m*x(n)*wrk(3,n)/wrk(4,n)   
      go to 70 
c  
c  use recursion to get plm from previous values   
c  
   30 lp1=lp+1 
      do 50 k=lp1,l
c  
c  increment storage indices   
c  
      iw1=iw2  
      iw2=iw3  
      iw3=1+mod(iw3,3) 
c  
c  set factors, depending on inorm 
c  
      if(inorm.eq.2) go to 35  
      fc1=(2*k-1)/float(k-m)   
      fc2=(k+m-1)/float(k-m)   
      go to 40 
c  
   35 fc1=(2*k+1)/float((k+m)*(k-m))   
      fc2=0
      if(k.gt.1) fc2=sqrt(fc1*(k-1+m)*(k-1-m)/float(2*k-3))
      fc1=sqrt(fc1*(2*k-1))
c  
c  set first term  
c  
   40 do 42 n=1,nx 
   42 wrk(iw3,n)=fc1*x(n)*wrk(iw2,n)   
c  
c  when k = m+1 this is all
c  
      if(k.eq.m+1) go to 50
c  
c  add second term 
c  
      do 44 n=1,nx 
   44 wrk(iw3,n)=wrk(iw3,n)-fc2*wrk(iw1,n) 
   50 continue 
c  
c  finally set derivative  
c  
c  set factor depending on inorm   
c  
      if(inorm.eq.2) go to 55  
      fcd=l+m  
      go to 60 
c  
   55 fcd=sqrt((2*l+1)*(l-m)*(l+m)/float(2*l-1))   
c  
   60 do 62 n=1,nx 
   62 wrk(4,n)=(fcd*wrk(iw2,n)-l*x(n)*wrk(iw3,n))/(1-x(n)*x(n))
c  
c  finish. store in plm and dplm   
c  
   70 call scopy(nx,wrk(iw3,1),iwrk,plm,iplm)  
      call scopy(nx,wrk(4,1),iwrk,dplm,idplm)  
c  
c  set current values and x
c  
      inormp=inorm 
      lp=l 
      mp=m 
      nxp=nx   
      call scopy(nx,x,1,wrk(5,1),iwrk) 
c  
      return   
      end  

      subroutine armax(x,a,ia,nn,xext,aext,d2aext)
c
c  subroutine armax finds the maximum aext of the function a = a(x),
c  and the point xext where the maximum occurs. x(n) contains x and
c  a(n) contains a, for n = 1,...,nn. ia is the first dimension of a.
c  entry armin similarly finds the minimum of a(x), and
c  entry aramax finds the maximum of abs(a(x)).
c  in each case the extremum is found by fitting a parabola through
c  the functional values at the three points (x(n),a(n)) closest to
c  the extremum, if this is not at n = 1 or n = nn.
c  if the extremum is in the interior of the interval considered
c  d2aext is set to the second derivative of a wrt x at the
c  extremum. otherwise d2aext is set to 0.
c
c  note that x, a, xext, aext and d2aext are real.
c
c  calls leq.
c  **********
c
c  ......................................................................
c
      logical max,abm
      dimension x(nn),a(ia,nn),wrk(3,3),b(3)
      max=.true.
      go to 10
c
      entry armin(x,a,ia,nn,xext,aext,d2aext)
      max=.false.
      go to 10
c
      entry aramax(x,a,ia,nn,xext,aext,d2aext)
      max=.true.
      abm=.true.
      go to 12
c
   10 abm=.false.
c  find extremal value of a
   12 axt=a(1,1)
      if(abm) axt=abs(axt)
      nxt=1
      do 20 n=2,nn
      aa=a(1,n)
      if(abm) aa=abs(aa)
      if(max) go to 15
      if(aa-axt) 17,17,20
   15 if(axt-aa) 17,17,20
   17 nxt=n
      axt=aa
   20 continue
      if(nxt.ne.1.and.nxt.ne.nn) go to 30
      xext=x(nxt)
      aext=axt
      d2aext=0
      return
c  fit parabola around x(nxt)
   30 xtr=x(nxt)
      scl=x(nxt+1)-x(nxt-1)
      if(scl) 35,90,35
   35 do 40 i=1,3
      t=1
      wrk(i,1)=t
      i1=nxt-2+i
      xx=(x(i1)-xtr)/scl
      do 38 k=2,3
      t=t*xx
   38 wrk(i,k)=t
   40 b(i)=a(1,i1)
c
      call leq(wrk,b,3,1,3,1,det)
c
c  find extremum
      if(det) 42,90,42
   42 if(b(3)) 45,95,45
   45 xext=-b(2)/(2*b(3))
      aext=b(1)+xext*(b(2)+xext*b(3))
      if(abm) aext=abs(aext)
      xext=xtr+scl*xext
      d2aext=2*b(3)/(scl*scl)
      return
   90 write(6,100) x(nxt-1),x(nxt),x(nxt+1)
      d2aext=0
      return
   95 write(6,110) b
      d2aext=0
      return
  100 format(///1x,10('*'),' points coincide in armax'/
     *  ' x(nxt-1),x(nxt),x(nxt+1):',1p3e15.5)
  110 format(///1x,10('*'),' coefficient to x**2 is zero in armax'/
     *  ' b =',1p3e15.5)
      end

      integer function intgps(x)
c  finds integer part of real*4 x, i.e. largest integer  .le. x
c  (note: differs from standard fortran function int for negative x)
c
c     integer function intgpd*4(x):  same for real*8 x
      real*8 xd
      x1=x
      ipr=1
      go to 10
c
      entry intgpd(xd)
      x1=xd
      ipr=2
c
   10 ix=int(x1)
      if(x1.lt.0.and.ix.ne.x1) ix=ix-1
c
      go to (15,20), ipr
   15 intgps=ix
      return
c
   20 intgpd=ix
      return
      end

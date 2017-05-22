      subroutine binave(x,y,nn,ix,iy,nmean,xm,ym,nnbin)
c
c  make binned averages of x(1,n) and y(1,n), in nmean uniform
c  bins in x. Averages are set to xm(k), ym(k), k = 1, ..., nmean.
c  Number of points in each bin are returned in nnbin(k).
c
c  Note: x need not be sorted.
c
c  Original version: 21/111986
c
      dimension x(ix,nn),y(iy,nn),xm(nmean),ym(nmean),nnbin(nmean)
c
c  set range in x.
c
      xmin=x(1,1)
      xmax=xmin
c
      do 20 n=2,nn
      xmin=amin1(xmin,x(1,n))
   20 xmax=amax1(xmax,x(1,n))
c
      dx=(xmax-xmin)/(nmean-1)
c
c  step through points, setting averages
c
   30 do 32 n=1,nmean
      xm(n)=0
      ym(n)=0
   32 nnbin(n)=0
c
      do 35 n=1,nn
      nb=(x(1,n)-xmin)/dx+1
      nb=min0(nb,nmean)
      nnbin(nb)=nnbin(nb)+1
      xm(nb)=xm(nb)+x(1,n)
   35 ym(nb)=ym(nb)+y(1,n)
c
c  set averages
c
      nb1=0
      do 40 nb=1,nmean
      nbin=nnbin(nb)
      if(nbin.gt.0) then
        nb1=nb1+1
        xm(nb1)=xm(nb)/nbin
        ym(nb1)=ym(nb)/nbin
        nnbin(nb1)=nbin
      end if
c
   40 continue
c
c  reset nmean to actual number of bins set
c
      nmean=nb1
      return
      end

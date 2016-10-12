      subroutine linint(x,y,rhs,ii,iy,ig,isn,nd1,nn)
c
c  linear differential equation solver. differs from original version
c  of linint by incorporating gaussian elimination without pivoting
c  in subroutine.
c
c  original version: 17/3/1986
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
      dimension x(nn),y(iy,nn),fd(20,20),w(20,20),y1(20)
      common/modfac/ r21,ntld
      external rhs
c
c  explicitly disable check for linear dependence
c
      ifd=20
c
      iw=20
c  right hand sides at first point
    5 n=1
      nc=1
      nd=nd1
      ntst=nn-2
      iig=ii*ig
      x2=x(1)
      call rhs(x2,fd,ifd,nd)
      go to 20
c  right hand sides at next point
   10 nd=nd+isn
      call rhs(x2,fd,ifd,nd)
c  set coefficient matrix for linear equations
      do 15 i=1,ii
      do 13 j=1,ii
   13 w(i,j)=dx*fd(i,j)
   15 w(i,i)=1+w(i,i)
c
c  solve linear equations by gaussian elimination without pivoting
c
      ii1=ii-1
      if(ii1.gt.0) then
c
c  triangularize matrix
c
        do 17 i=1,ii1
        r=w(i,i)
c
c  test for non-zero diagonal element
c
        if(r.eq.0) then
          write(6,110) i,n
          return
        end if
c
        r=-1./r
c
        i1=i+1
c
        do 17 j=i1,ii
        rj=r*w(j,i)
        do 16 k=i1,ii
   16   w(j,k)=w(j,k)+rj*w(i,k)
c
        js=j
        is=i
c
        do 17 l=1,ig
        y(js,n)=y(js,n)+rj*y(is,n)
        js=js+ii
   17   is=is+ii
c
      end if
c
c  now matrix is on triangular form
c
c  start solution
c
      i=ii+1
      do 18 icnt=1,ii
      i=i-1
c
      r=w(i,i)
c
c  test for non-zero diagonal element
c
      if(r.eq.0) then
        write(6,110) i,n
        return
      end if
c
      r=1./r
      is=i
      i1=i+1
c
      do 18 l=1,ig
      sum=y(is,n)
c
      if(i.lt.ii) then
c
        j1=is
        do 17100 j=i1,ii
        j1=j1+1
17100   sum=sum-w(i,j)*y(j1,n)
c
      end if
c
      y(is,n)=r*sum
c
   18 is=is+ii
      if(nc.eq.nn) return
c
c  set up right hand side of linear equation at next point
c
   20 n1=n
      n=n+isn
      nc=nc+1
      x1=x2
      x2=x(n)
      dx=0.5*(x1-x2)
      kg=-ii
      do 25 l=1,ig
      kg=kg+ii
      do 25 i=1,ii
      k=kg+i
      sum=0
      do 22 j=1,ii
   22 sum=sum+fd(i,j)*y(kg+j,n1)
   25 y(k,n)=y(k,n1)-dx*sum
      if(ntld.ne.1) go to 10
      if(nc.ne.ntst.or.ig.ne.2) go to 10
c  test for linear dependence (only implemented for ig = 2)
      do 30 i=1,iig
   30 y1(i)=y(i,n1)
      aym=0
      im=0
      aym2=0
      do 32 i=1,ii
      ay=abs(y1(i))
      aym2=max(aym2,abs(y1(i+ii)))
      if(ay.le.aym) go to 32
      aym=ay
      im=i
   32 continue
c
      r21=y1(ii+im)/y1(im)
      s1=0
      s2=0
      j=ii
      do 34 i=1,ii
      j=j+1
      yi=y1(i)
      if(yi.ne.0) go to 33
      if(abs(y1(j)).gt.1.e-5*aym2) go to 10
      go to 34
   33 ri=y1(j)/y1(i)
      s1=s1+abs(r21-ri)
      s2=s2+abs(ri)
   34 continue
c  test
      if(s1/s2.gt.1.e-7) go to 10
c  nearly linear dependence. modify initial values and start again
      j=0
      aym=0
      do 36 i=1,ii
      y1(i)=y(i+ii,1)-r21*y(i,1)
   36 aym=max(aym,abs(y1(i)))
c  normalize initial y's
      aym=1./aym
      do 38 i=1,ii
   38 y(ii+i,1)=aym*y1(i)
c
      write(6,100) (y(i,1),i=1,iig)
      go to 5
c
  100 format(//' initial solution has been modified. new y(n=1):'/
     *  (1p10e13.5))
c
  110 format(//'  ***** in s/r linint zero diagonal element at i =',
     *  i5,'   n =',i5)
      end

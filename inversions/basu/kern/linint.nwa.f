      subroutine linint(x,y,rhs,ii,iy,ig,isn,nd1,nn)
c
c  linear differential equation solver. differs from original version
c  of linint by incorporating gaussian elimination without pivoting
c  in subroutine.
c
c  as linint.new. but using special Alliant routines for gaussian
c  elimination
c
c  original version: 17/3/1986
c
      dimension x(nn),y(iy,nn),fd(400),w(400),y1(20),w1(400)
      common/modfac/ r21,ntld
      external rhs
c
c  initialize w1 to delta(i,j)
c
      ii2=ii*ii
      do 2 i=1,ii2
    2 w1(i)=0
      do 3 i=1,ii
    3 w1(i+ii*(i-1))=1
c  right hand sides at first point
    5 n=1
      nc=1
      nd=nd1
      ntst=nn-2
      iig=ii*ig
      x2=x(1)
      call rhs(x2,fd,ii,nd)
      go to 20
c  right hand sides at next point
   10 nd=nd+isn
      call rhs(x2,fd,ii,nd)
c  set coefficient matrix for linear equations
c..      do 15 i=1,ii
c..      do 13 j=1,ii
c..   13 w(i,j)=dx*fd(i,j)
c..   15 w(i,i)=1+w(i,i)
      call scopy(ii2,w1,1,w,1)
      call saxpy(ii2,dx,fd,1,w,1)
c
c  solve linear equations by gaussian elimination without pivoting
c
      ii1=ii-1
      if(ii1.gt.0) then
c
c  triangularize matrix
c
        do 17 i=1,ii1
        idiag=i+ii*(i-1)
        r=w(idiag)
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
        icols=ii-i
        i1=i+1
        ji=idiag
c
cvd$ cncall
        do 17 j=i1,ii
        ji=ji+1
        rj=r*w(ji)
c..        do 16 k=i1,ii
c..   16   w(j,k)=w(j,k)+rj*w(i,k)
        call saxpy(icols,rj,w(idiag+ii),ii,w(ji+ii),ii)
c
c..        js=j
c..        is=i
c..c
c..        do 17 l=1,ig
c..        y(js,n)=y(js,n)+rj*y(is,n)
c..        js=js+ii
c..   17   is=is+ii
        call saxpy(ig,rj,y(i,n),ii,y(j,n),ii)
   17   continue
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
      idiag=i+ii*(i-1)
c
      r=w(idiag)
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
cvd$ nocncall
      do 18 l=1,ig
      sum=y(is,n)
c
      if(i.lt.ii) then
c
        j1=is
c..        do 17100 j=i1,ii
c..        j1=j1+1
c..17100   sum=sum-w(i,j)*y(j1,n)
c
        sum=sum-sdot(ii-i,w(idiag+ii),ii,y(is+1,n),1)
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
c..      sum=0
c..      do 22 j=1,ii
c..   22 sum=sum+fd(i,j)*y(kg+j,n1)
c..   25 y(k,n)=y(k,n1)-dx*sum
c
   25 y(k,n)=y(k,n1)-dx*sdot(ii,fd(i),ii,y(kg+1,n1),1)
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
      aym2=amax1(aym2,abs(y1(i+ii)))
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
   36 aym=amax1(aym,abs(y1(i)))
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

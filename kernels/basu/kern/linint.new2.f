      subroutine linint(x,y,rhs,ii,iy,ig,isn,nd1,nn)
c
c  linear differential equation solver. differs from original version
c  of linint by incorporating gaussian elimination without pivoting
c  in subroutine.
c
c  original version: 17/3/1986
c
      dimension x(nn),y(1),fd(400),w(400),y1(20)
      common/modfac/ r21,ntld
      external rhs
c
      isny=isn*iy
c
c  right hand sides at first point
c
    5 n=1
      ny=0
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
      id=1
      do 15 i=1,ii
      ij=i
      do 13 j=1,ii
      w(ij)=dx*fd(ij)
   13 ij=ij+ii
      w(id)=1+w(id)
   15 id=id+ii+1
c
c  solve linear equations by gaussian elimination without pivoting
c
      ii1=ii-1
      if(ii1.gt.0) then
c
c  triangularize matrix
c
        id=1
        do 17 i=1,ii1
        r=w(id)
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
        ji=id+1
        do 16700 j=i1,ii
        rj=r*w(ji)
        jk=ji+ii
        ik=id+ii
        do 16 k=i1,ii
        w(jk)=w(jk)+rj*w(ik)
        jk=jk+ii
   16   ik=ik+ii
c
        js=j+ny
        is=i+ny
c
        do 16500 l=1,ig
        y(js)=y(js)+rj*y(is)
        js=js+ii
16500   is=is+ii
c
16700   ji=ji+1
c
   17   id=id+ii+1
c
      end if
c
c  now matrix is on triangular form
c
c  start solution
c
      i=ii
      id=ii*ii
      do 18 icnt=1,ii
c
      r=w(id)
      i1=i+1
c
c  test for non-zero diagonal element
c
      if(r.eq.0) then
        write(6,120) i,n
        return
      end if
c
      r=1./r
      is=i+ny
c
      do 17200 l=1,ig
      sum=y(is)
c
      if(i.lt.ii) then
c
        j1=is+1
        ij=id+ii
        do 17100 j=i1,ii
        sum=sum-w(ij)*y(j1)
        ij=ij+ii
17100   j1=j1+1
c
      end if
c
      y(is)=r*sum
c
17200 is=is+ii
c
      id=id-ii-1
   18 i=i-1
      if(nc.eq.nn) return
c
c  set up right hand side of linear equation at next point
c
   20 n1=n
      ny1=ny
      n=n+isn
      ny=ny+isny
      nc=nc+1
      x1=x2
      x2=x(n)
      dx=0.5*(x1-x2)
      kg=ny1
c
      do 25 l=1,ig
      do 24 i=1,ii
      ij=i
      k=kg+i
      sum=0
      do 22 j=1,ii
      sum=sum+fd(ij)*y(kg+j)
   22 ij=ij+ii
   24 y(k+isny)=y(k)-dx*sum
   25 kg=kg+ii
      if(ntld.ne.1) go to 10
      if(nc.ne.ntst.or.ig.ne.2) go to 10
c
c  test for linear dependence (only implemented for ig = 2)
c
      do 30 i=1,iig
   30 y1(i)=y(i+ny1)
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
      y1(i)=y(i+ii)-r21*y(i)
   36 aym=amax1(aym,abs(y1(i)))
c  normalize initial y's
      aym=1./aym
      do 38 i=1,ii
   38 y(ii+i)=aym*y1(i)
c
      write(6,100) (y(i),i=1,iig)
      go to 5
c
  100 format(//' initial solution has been modified. new y(n=1):'/
     *  (1p10e13.5))
c
  110 format(//'  ***** in s/r linint zero diagonal element at i =',
     *  i5,'   n =',i5,'  in triangularization')
  120 format(//'  ***** in s/r linint zero diagonal element at i =',
     *  i5,'   n =',i5,'  in back substitution')
      end

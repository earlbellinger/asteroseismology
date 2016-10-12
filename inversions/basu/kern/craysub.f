      function sdot(n,a,i1,b,i2)
      dimension a(1),b(1)
      sdot=0.
      j1=1
      j2=1
      do 1 i=1,n
      sdot=sdot+a(j1)*b(j2)
      j1=j1+i1
      j2=j2+i2
    1 continue
      return
      end 
      function ssum(n,a,i1)
      dimension a(1)
      ssum=0.
      j1=1
      do 1 i=1,n
      ssum=ssum+a(j1)
      j1=j1+i1
    1 continue
      return
      end 
      subroutine sscal(n,fak,a,i1)
      dimension a(1)
      j1=1
      do 1 i=1,n
      a(j1)=fak*a(j1)
    1 j1=j1+i1
      return
      end 
      integer function ismin(n,a,nstep) 
      dimension a(1)
      ndim=n*nstep
      ismin=1
      k=1 
      x=a(1)
      do 1 i=1,ndim,nstep
      if(a(i).ge.x)go to 1
      ismin=k
      x=a(i)
    1 k=k+1
      return
      end 
      subroutine scopy(n,a,na,b,nb)
      dimension a(1),b(1)
      ia=1
      ib=1
      do 1 i=1,n
      b(ib)=a(ia)
      ia=ia+na
      ib=ib+nb
    1 continue
      return
      end 
      subroutine ccopy(n,a,na,b,nb)
      complex a(1),b(1)
      ia=1
      ib=1
      do 1 i=1,n
      b(ib)=a(ia)
      ia=ia+na
      ib=ib+nb
    1 continue
      return
      end 
      integer function ismax(n,a,nstep) 
      dimension a(1)
c
      ndim=n*nstep
      x=a(1)
      ismax=1
      k=1 
      do 1 i=1,ndim,nstep
      if(a(i).le.x)go to 1
      ismax=k
      x=a(i)
    1 k=k+1
      return
      end 
      integer function isamax(n,a,nstep) 
      dimension a(1)
c
      ndim=n*nstep
      x=abs(a(1))
      isamax=1
      k=1 
      do 1 i=1,ndim,nstep
      aa=abs(a(i))
      if(aa.le.x)go to 1
      isamax=k
      x=aa
    1 k=k+1
      return
      end 
      subroutine saxpy(n,c,a,na,b,nb)
      dimension a(1),b(1)
      ia=1
      ib=1
      do 1 i=1,n
      b(ib)=c*a(ia)+b(ib)
      ia=ia+na
    1 ib=ib+nb
      return
      end 
      subroutine saxpyd(n,c,a,na,b,nb)
c
c  double precision version of SCILIB emulator saxpy
c
      double precision c,a,b
      dimension a(1),b(1)
      ia=1
      ib=1
      do 1 i=1,n
      b(ib)=c*a(ia)+b(ib)
      ia=ia+na
    1 ib=ib+nb
      return
      end 
      subroutine cscal(n,cfak,a,na)
      complex cfak,a(1)
      ndim=n*na
      do 1 i=1,ndim,na
    1 a(i)=cfak*a(i)
      return
      end 
      subroutine caxpy(n,cx,a,na,b,nb)
      complex cx,a(1),b(1)
      ia=1
      ib=1
      do 1 i=1,n
      b(ib)=cx*a(ia)+b(ib)
      ia=ia+na
      ib=ib+nb
    1 continue
      return
      end 
      complex function cdotu(n,a,na,b,nb)
      complex a(1),b(1)
      ia=1
      ib=1
      cdotu=0.
      do 1 i=1,n
      cdotu=cdotu+a(ia)*b(ib) 
      ia=ia+na
    1 ib=ib+nb
      return
      end 
      function cvmgp(x1,x2,x3)
      if(x3)1,2,2
    1 cvmgp=x2
      return
    2 cvmgp=x1
      return
      end 
      function cvmgz(x1,x2,x3)
      if(x3)1,2,1
    1 cvmgz=x2
      return
    2 cvmgz=x1
      return
      end 
      function cvmgn(x1,x2,x3)
      if(x3)2,1,2
    1 cvmgn=x2
      return
    2 cvmgn=x1
      return
      end 
      function cvmgt(x1,x2,x3)
      logical x3
      if(x3) then
        cvmgt=x1
      else
        cvmgt=x2
      end if
      return
      end 
      function cvmgm(x1,x2,x3)
      if(x3)2,1,1
    1 cvmgm=x2
      return
    2 cvmgm=x1
      return
      end 
      integer function idamax(n,a,nstep) 
      implicit double precision(a-h,o-z)
      dimension a(1)
      data initpr /0/
c
      if(initpr.eq.0) then
        initpr=1
        write(6,99999) 
99999   format(//' ***** emulate idamax *****'//)
      end if
c
      ndim=n*nstep
      x=abs(a(1))
      idamax=1
      k=1 
      do 1 i=1,ndim,nstep
      aa=abs(a(i))
      if(aa.le.x)go to 1
      idamax=k
      x=aa
    1 k=k+1
      return
      end 
      subroutine daxpy(n,c,a,na,b,nb)
      implicit double precision(a-h,o-z)
      dimension a(1),b(1)
      data initpr /0/
c
      if(initpr.eq.0) then
        initpr=1
        write(6,99999) 
99999   format(//' ***** emulate daxpy *****'//)
      end if
      ia=1
      ib=1
      do 1 i=1,n
      b(ib)=c*a(ia)+b(ib)
      ia=ia+na
    1 ib=ib+nb
      return
      end 
      subroutine dscal(n,fak,a,i1)
      implicit double precision(a-h,o-z)
      dimension a(1)
      data initpr /0/
c
      if(initpr.eq.0) then
        initpr=1
        write(6,99999) 
99999   format(//' ***** emulate dscal *****'//)
      end if
      j1=1
      do 1 i=1,n
      a(j1)=fak*a(j1)
    1 j1=j1+i1
      return
      end 
      function ddot(n,a,i1,b,i2)
c
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
      dimension a(1),b(1)
      data initpr /0/
c
      if(initpr.eq.0) then
        initpr=1
        write(6,99999) 
99999   format(//' ***** emulate ddot *****'//)
      end if
      ddot=0.
      j1=1
      j2=1
      do 1 i=1,n
      ddot=ddot+a(j1)*b(j2)
      j1=j1+i1
      j2=j2+i2
    1 continue
      return
      end 

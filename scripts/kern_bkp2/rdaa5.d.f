      subroutine rdaa5(ichaa,data,x,aa,nn,iend)
      implicit double precision (a-h,o-z)
      logical diag
      parameter(nnmax=10000)
      dimension x(nnmax),aa(5,nnmax),data(8)
      common/work/aa1(10,nnmax)
      diag = .true.
      iend=0
      read(ichaa,end=90) nmod,nn,data
      close(ichaa)
      idata8=nint(data(8))
      if(idata8.lt.10) then
	icase = 0
	ivar = 5
      else if(idata8.lt.100) then
	icase = 1
	ivar = 6
      else
	icase = 2
	ivar = 8
      end if
      write(6,*) 'In rdaa5  nn, ivar =', nn, ivar
c
      call openf(ichaa,'o','u')
      ivar1=ivar+1
      read(ichaa,end=90) nmod,nn,data,((aa1(i,n),i=1,ivar1),n=1,nn)
         do 50 l1=1,nn
         x(l1) = aa1(1,l1)
            do 40 l2=1,5
            aa(l2,l1)=aa1(l2+1,l1)
   40       continue
   50    continue
      if (diag) write(6,1001) (data(i),i=1,8)
      write(6,*) 'After read of aa, nn, ivar =',nn, ivar
c
c     following diagnostic inserted because of confusion over
c     contents of aa-file.  aa(3,.) should contain gamma1.
c
      diag = .true.
      if (diag) write (6,1000) nn,aa(3,1),aa(3,nn)
      if (diag) write (6,1000) nn,aa(1,1),aa(1,nn)
      diag = .false.
      return
   90 iend=999
      return
 1000 format(//i5,' mesh points in aa-5 model'/
     +  ' aa(3,.) at innermost/outermost meshpoint = ',1p2e15.7,0p//)
 1001 format(//'data: '/(4e14.5))
      end

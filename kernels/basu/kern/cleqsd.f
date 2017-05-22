      subroutine cleqsd(a,b,nn,mm,ia,ib,err,iad,ibd)
      implicit complex (a-h,o-z)
      real r,r1,ri
      logical norhs
      dimension a(iad),b(ibd)
      data one,zero/(1.,0.),(0.,0.)/
c
c     this routine will find the inverse of a matrix by the method of
c     partial pivoting and gaussian elimination if b is set equal to
c     the identity matrix
c
c     n - dimension of segment of a to be used
c     m - number of right hand columns of b to be used
c     ia - the total number of rows in large array a
c     ib - the total number of rows in large array b
c     the matrix equation is    ax=b
c     err = det(a)
c
c     if mm = 0 cleqsd calculates err = det(a)
c
c      iad  -  total size of a
c
c      ibd  -  total size of b
c
      n=nn
      m=mm
c
      err=zero
      norhs=m.le.0
      detsc=one
c
c
      n1=ia*n
      ia1=ia+1
      m1=ib*m
c
c     n.le.1
c
      if(n-1) 280,281,285
c     no equation
  280 write(6,295) n
      err=zero
      return
c     n = 1
  281 err=a(1)
      ai=one/err
      do 282 j=1,m1,ib
  282 b(j)=b(j)*ai
      return
c
c     find maximum element in each row and divide the row by it
c
  285 do 1 i=1,n
      r=cabs(a(i))
      do 2 j=  1,n1,ia
      ij=j+i-1
    2 r=amax1(r,cabs(a(ij)))
      if(r)31,30,31
   30 write(6,298)i
      return
   31 ri=1./r
      do 3 j=1,n1,ia
      ij=j+i-1
    3 a(ij)=a(ij)*ri
      if(norhs) go to 1
      do 4 j=1,m1,ib
      ij=j+i-1
    4 b(ij)=b(ij)*ri
    1 detsc=detsc*r
c
c
c     find maximum element in the i'th column
      n2=n-1
      do 5 i=1,n
      ialow=(i-1)*ia+i
      iaup=(i-1)*ia+n
      r=cabs(a(ialow))
      imax=ialow
      if(i.eq.n) go to 81
      ialow1=ialow+1
      do 6 j=ialow1,iaup
      r1=cabs(a(j))
      if(r-r1)7,6,6
    7 imax=j
      r=r1
    6 continue
c     test that maximum element is non-zero
   81 if(r) 32,32,33
   32 write(6,299) i,i
      return
   33 if(imax-ialow)8,8,9
c     replace the i'th row with the row that has the maximum element in
c         the respective column and put the i'th row in its place
    9 im=imax
   72 if(im-ia)70,70,71
   71 im=im-ia
      go to 72
   70 do 10 j=1,n1,ia
      jj=i+j-1
      ji=im+j-1
      rr=a(jj)
      a(jj)=a(ji)
   10 a(ji)=rr
c     change sign of determinant
      detsc=-detsc
c
      if(norhs) go to 8
c
      do 11 j=1,m1,ib
      jj=i+j-1
      ji=im+j-1
      rr=b(jj)
      b(jj)=b(ji)
   11 b(ji)=rr
c     multiply the i'th row by (the negative of each i'th column element
c       below the diagonal divided by the diabonal element) and add the
c     resulting row to the respective row of the element used
    8 iaup1=iaup-1
c
c
c     contribution to determinant
      detsc=detsc*a(ialow)
c     test new size of determinant (it may have become zero due to
c     underflow)
      r=cabs(detsc)
      if(r) 32,32,83
c     invert diagonal element
   83 a(ialow)=one/a(ialow)
c
      if(i.eq.n) go to 5
c
      do 12 j=ialow,iaup1
      j1=j+1
      rr=-a(j1)*a(ialow)
      i1=ialow
      do 13 k=j1,n1,ia
      a(k)=a(k)+a(i1)*rr
   13 i1=i1+ia
      if(norhs) go to 12
      i1=i
      i2=j-ialow+i+1
      do 14 k=i2,m1,ib
      b(k)=b(k)+b(i1)*rr
   14 i1=i1+ib
   12 continue
c
c
    5 continue
c
c
c     the matrix is now in triangular form
c
c     set determinant
      err=detsc
c
      if(norhs) return
c     find solution to ax=b
      do 18 k=1,m
      ka=(n-1)*ia+n
      kb=(k-1)*ib+n
      b(kb)=b(kb)*a(ka)
      do 19 l=1,n2
      i=n-l
      rr=zero
      imax=i+1
      do 20 j=imax,n
      jj=i+n+1-j
      ja=(jj-1)*ia+i
      jb=(k-1)*ib+jj
   20 rr=rr+a(ja)*b(jb)
      la=(i-1)*ia+i
      lb=(k-1)*ib+i
   19 b(lb)=(b(lb)-rr)*a(la)
   18 continue
c
c
      return
  295 format(///4h n =,i10,13h. no equation)
  298 format(26h1  all the elements in row,i5,20h  are zero therefore,
     &55h the solution, matrix x, cannot be found by this method )
  299 format(50h1  the solution, matrix x, cannot be found by this  ,
     &57h method because there is a zero array element in the main ,
     &9h diagonal  / 30x,2ha(,i4,1h,,i4,8h) = zero  )
      end

      subroutine istore(ia,ib,n)
c
c     stores first n elements of integer a into integer b
c
      dimension ia(n), ib(n)
c
   10 do 11 i=1,n
   11 ib(i)=ia(i)
      return
c
      end

      subroutine leqsd(a,b,n,m,ia,ib,det,isa,isb)
      implicit double precision (a-h,o-z)
      dimension a(ia,1),b(ib,1)
c
      call leq(a,b,n,m,ia,ib,det)
      return
      end

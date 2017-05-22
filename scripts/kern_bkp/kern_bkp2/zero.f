      subroutine  zero(a,nn)
c
c
c     sets a(n)=0., n=1,nn
c
c
      dimension a(nn)
c
c
      do 1 n=1,nn
    1 a(n)=0.
      return
c
      end

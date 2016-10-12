      subroutine store(a,b,n)
c
c     stores first n elements of single precision a 
c        into single precision b
c
      dimension a(n),b(n)
c
   10 do 11 i=1,n
   11 b(i)=a(i)
      return
c
      end
c
      subroutine stsd(a,db,n)
c     stores first n elements of single precision a 
c        in double precision db
c
      double precision db
      dimension a(n),db(n)
c
   10 do 11 i=1,n
   11 db(i)=a(i)
      return
c
      end
c
      subroutine stds(da,b,n)
c     stores first n elements of double precision da 
c        in single precision b
      double precision da
      dimension b(n),da(n)
c
   10 do 11 i=1,n
   11 b(i)=da(i)
      return
c
      end
c
      subroutine stdd(da,db,n)
c     stores first n elements of double precision da 
c        in double precision db
      double precision da,db
      dimension da(n),db(n)
c
   10 do 11 i=1,n
   11 db(i)=da(i)
      return
c
      end

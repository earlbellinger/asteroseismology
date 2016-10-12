      subroutine lexsrt(a,ia,nn,in,ida)
      dimension a(ida,nn),in(nn)
      dimension w(1000),fn(1000)
c  find index permutation so that (a(i,in(n)),i=1,ia) is sorted
c  lexicographically.
      do 10 n=1,nn
      w(n)=a(1,n)
   10 fn(n)=n
      call sort(w,fn,nn)
c
      do 12 n=1,nn
   12 in(n)=fn(n)
      if(ia.eq.1) return
c  step up in i
      do 50 i=2,ia
c  find points with equal values of a(k,.), k=1,...,i-1
      j=i-1
      n1=1
   15 m=0
      do 20 n=n1,nn
      mp=m
      m=in(n)
      if(n.eq.n1) go to 18
      do 17 k=1,j
   17 if(a(k,m).ne.a(k,mp)) go to 25
   18 w(n)=a(i,m)
      fn(n)=m
   20 n2=n
c  sort after a(i,.) for n between n1 and n2
   25 if(n2.eq.n1) go to 30
      nns=n2-n1+1
      call sort(w(n1),fn(n1),nns)
      do 27 n=n1,n2
   27 in(n)=fn(n)
c
   30 if(n2.eq.nn) go to 50
      n1=n2+1
      go to 15
   50 continue
      return
      end

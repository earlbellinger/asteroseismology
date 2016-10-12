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

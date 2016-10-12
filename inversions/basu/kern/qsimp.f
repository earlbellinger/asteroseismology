      function qsimp(func,a,b,eps)
c
c  evaluates integral of function func(x) from a to b with 
c  accuracy eps
c
      parameter (jmax=20)
      ost=-1.e30
      os= -1.e30
      do 11 j=1,jmax
        call trapzd(func,a,b,st,j)
        s=(4.*st-ost)/3.
        if (abs(s-os).lt.eps*abs(os)) then
	  qsimp=s
	  return
	end if
        os=s
        ost=st
11    continue
      write(6,110)
      stop 
  110 format(/' *** too many steps in qsimp.')
      end
      subroutine trapzd(func,a,b,s,n)
      external func
      if (n.eq.1) then
        s=0.5*(b-a)*(func(a)+func(b))
        it=1
      else
        tnm=it
        del=(b-a)/tnm
        x=a+0.5*del
        sum=0.
        do 11 j=1,it
          sum=sum+func(x)
          x=x+del
11      continue
        s=0.5*(s+(b-a)*sum/tnm)
        it=2*it
      endif
      return
      end



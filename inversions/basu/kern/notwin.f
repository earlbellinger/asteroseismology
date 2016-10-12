      logical function notwin(xw1,xw2,x) 
c
c  windowing function
c  returns true, if x is not in window defined by xw1 and xw2
c
      if(xw1.gt.xw2) then
        notwin = .false.
      else
        notwin = x.lt.xw1.or.x.gt.xw2
      end if
      return
      end

      function loctim
      character*(*) loctim
      character*24 timest,fdate
c
c  returns date and time as a string. Note this routine
c  may depend on installation.
c
      timest= fdate()
      loctim=timest
      return
      end

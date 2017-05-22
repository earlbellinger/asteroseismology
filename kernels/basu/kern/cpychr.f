      subroutine cpychr(idsin,idsout)
c
c  copies file from unit idsin to unit idsout, both assumed
c  to be opened, truncating blankt at end of lines
c  Lines have to be less than 256 characters
c
c  Original version: 7/5/92
c
      character s*256
c
   10 read(idsin,100,end=20) s
      ls=length(s)
      write(idsout,100) s(1:ls)
      go to 10
c
   20 continue
      return
  100 format(a)
      end

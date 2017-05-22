      subroutine cdate(idat)
c
c  emulated Cray routine date on Alliant
c
      character*(*) idat
      dimension iadat(3)
c
      call idate(iadat)
      iadat(3)=mod(iadat(3),100)
      write(idat,100) iadat(2), iadat(1), iadat(3)
      return
  100 format(i2.2,'/',i2.2,'/',i2.2)
      end
      subroutine clock(itim)
c
c  emulates Cray routine clock on Alliant
c
      character*(*) itim
      dimension iatim(3)
c
      call itime(iatim)
      write(itim,100) iatim
      return
  100 format(i2.2,':',i2.2,':',i2.2)
      end

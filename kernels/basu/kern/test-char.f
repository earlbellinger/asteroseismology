      program main
c
      character*40 a, b
      data a /'1234567891234567891234567891234567890'/
c
      do 10 i=10,50,5
      b = a
      b(i:i) = char(0)
      write(6,*) b, i
   10 continue
      stop
      end

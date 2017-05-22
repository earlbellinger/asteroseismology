c
c     reads kernel from unformatted file and produces formatted output
c
      parameter(nnmax=1500)
      implicit double precision (a-h,o-z)
      dimension x(nnmax),rker(nnmax)
      character*1 yesno
      character*32 filnam,filnam1

      write(*,*) 'Which unformatted file contains kernel?'
      read(*,1001) filnam
      open(8,file=filnam,status='old',form='unformatted')

      write(*,*) 'Where to put formatted kernel output?'
      read(*,1001) filnam1
      open(9,file=filnam1,status='unknown',form='formatted')

      read(8,err=900,end=950) nn,(x(i),i=1,nn)
      if (nn.gt.nnmax) goto 920
      write(*,*) ''
      write(*,*) 'No. of mesh points is ',nn
      write(*,*) ''
  100 read(8,err=900,end=950) iel,iord,rnu,(rker(i),i=1,nn)
      write(*,*) 'l,n,freq = ',iel,iord,rnu
  110 write(*,*) 'Do you want this kernel in ',filnam1,'? (y/n)'
      read(*,1002) yesno
      if (yesno.eq.'Y'.or.yesno.eq.'y') goto 150
      if (yesno.eq.'N'.or.yesno.eq.'n') then
         goto 100
      else 
         goto 110
      endif 
  150 write(9,*) '# Kernel from file ',filnam
      write(9,*) '#  l,n,freq = ',iel,iord,rnu
         do 200 i=1,nn
         write(9,*) x(i),rker(i)
  200    continue
      write(*,*) 'OK'
      goto 100

  900 write(*,*) 'Error reading from unformatted kernel file'
      stop
  920 write(*,*) 'Error! nnmax too small  --  nn,nnmax = ',nn,nnmax
      stop
  950 stop
 1001 format(1a32)
 1002 format(1a1)
      end 

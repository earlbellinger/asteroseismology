      implicit double precision (a-h,o-z)
      logical diag
      character*1 noyes
      dimension delta(2,2000), rk(2,2000)
      dimension x(2000),xgg(2000),deltgg(2,2000)
      real rrnu, rx(2000)
      common/consts/ fourpi,gconst
      common/work/dummy(10000)
      diag = .false.
      call ofiles
      call openf(9,'o','f')
      call openf(11,'o','u')
      call openf(12,'o','u')
      call openf(14,'u','f')
      open(20,file='diagc.out',status='unknown')
      open(21,file='diagrho.out',status='unknown')
      call skpcom(9)
   80 write(6,*) '22, 19 or 7 column format differences? (22/19/7)'
      read(5,*,err=80) ncol
      write(*,*) ncol
      if (ncol.ne.22.and.ncol.ne.19.and.ncol.ne.7) goto 80
 1005 format(1a1)
      write(*,*) 'Input 0 for ln(c),     ln(rho)'
      write(*,*) '      1 for ln(c^2),   ln(rho)'
      write(*,*) '      2 for ln(Gamma1),ln(rho)'
      write(*,*) '      3 for ln(c^2),   ln(Gamma1)'
      write(*,*) '      4 for ln(u),   ln(Gamma1)'
      write(*,*) '      5 for ln(u),   Y'
      write(*,*) '      6 for ln(rho),   Y'
      write(*,*) '      7 for ln(c),   ln(Gamma1/c)'
      read(*,*) icase
      write(*,*) icase
      ichdgg = 9
      ichkr1= 11
      ichkr2= 12
      gconst = 6.6732e-8
      fourpi = 8.*acos(0.0)
c
      write(14,*) '# Convolutions of kernels with model differences'
      if (icase.eq.0) then
         write(14,*) '# for variables ln(c) and ln(rho)'
      elseif (icase.eq.1) then
         write(14,*) '# for variables ln(c^2) and ln(rho)'
      elseif (icase.eq.2) then
         write(14,*) '# for variables ln(Gamma1) and ln(rho)'
      elseif (icase.eq.3) then
         write(14,*) '# for variables ln(c^2) and ln(Gamma1)'
      elseif (icase.eq.4) then
         write(14,*) '# for variables ln(u) and ln(Gamma1)'
      elseif (icase.eq.5) then
         write(14,*) '# for variables ln(u) and y'
      elseif (icase.eq.6) then
         write(14,*) '# for variables ln(rho) and y'
      elseif (icase.eq.7) then
         write(14,*) '# for variables ln(c) and ln(Gamma1/c)'
      endif
      write(14,*) '#    l     n     freq      ',
     +        ' From Ker1        Ker2           Total' 
      read(ichkr1) nn,(rx(l1),l1=1,nn)
         do 800 l1=1,nn
         x(l1)=rx(l1)
  800    continue
      nngg = 0
      if (ncol.eq.7) goto 92
      ioff = 2
      if (ncol.eq.19) ioff = 1
   90 read(ichdgg,*,end=95) (dummy(l1),l1=1,19)
      nngg = nngg+1
      xgg(nngg) = dummy(1)
      if (icase.eq.0) then
        deltgg(1,nngg) = dummy(15+ioff)
        deltgg(2,nngg) = dummy(7+ioff)
      endif
      if (icase.eq.1) then
        deltgg(1,nngg) = dummy(15+ioff)*2.
        deltgg(2,nngg) = dummy(7+ioff)
      endif
      if (icase.eq.2) then
        deltgg(1,nngg) = dummy(10+ioff)
        deltgg(2,nngg) = dummy(7+ioff)
      endif
      if (icase.eq.3) then
        deltgg(1,nngg) = dummy(15+ioff)*2.
        deltgg(2,nngg) = dummy(10+ioff)
      endif
      if (icase.eq.4) then
        deltgg(1,nngg) = dummy(15+ioff)*2.-dummy(10+ioff)
        deltgg(2,nngg) = dummy(10+ioff)
      endif
      if (icase.eq.5) then
        deltgg(1,nngg) = dummy(15+ioff)*2.-dummy(10+ioff)
        deltgg(2,nngg) = -dummy(5+ioff)
      endif
      if (icase.eq.6) then
        deltgg(1,nngg) = dummy(7+ioff)
        deltgg(2,nngg) = -dummy(5+ioff)
      endif
      if (icase.eq.7) then
        deltgg(1,nngg) = dummy(15+ioff)
        deltgg(2,nngg) = dummy(10+ioff)-dummy(15+ioff)
      endif
      goto 90

   92 read(ichdgg,*,end=95) (dummy(l1),l1=1,7)
      nngg = nngg+1
      xgg(nngg) = dummy(1)
      if (icase.eq.0) then 
        deltgg(1,nngg) = dummy(7)
        deltgg(2,nngg) = dummy(5)
      endif
      if (icase.eq.1) then 
        deltgg(1,nngg) = dummy(7)*2.
        deltgg(2,nngg) = dummy(5)
      endif
      if (icase.eq.2) then 
        deltgg(1,nngg) = dummy(7)*2.+dummy(5)-dummy(4)
        deltgg(2,nngg) = dummy(5)
      endif
      if (icase.eq.3) then 
        deltgg(1,nngg) = dummy(7)*2.
        deltgg(2,nngg) = dummy(7)*2.+dummy(5)-dummy(4)
      endif
      if (icase.eq.4) then 
        deltgg(1,nngg) = dummy(4)-dummy(5)
        deltgg(2,nngg) = dummy(7)*2.+dummy(5)-dummy(4)
      endif
      if (icase.eq.5.or.icase.eq.6) then
         write(*,*) 'Sorry, not programmed to test ncol=7 Y kernels'
         stop
      endif
      goto 92

   95 continue
      
      write(*,*) 'Does second kernel file also contain mesh & nl-s',
     +  '? (y/n)'
      read(*,*) noyes
      write(*,*) noyes
      if (noyes.eq.'Y') noyes = 'y'
      if (noyes.eq.'y') then
         read(ichkr2) nn,(rx(l1),l1=1,nn)
      endif
      ii = 2
      id = 2
         do 100 l1=1,nn
         z = x(l1)
         call lir(z,xgg,delta(1,l1),deltgg,ii,id,nngg,l1,inter)
  100    continue
  105 read(ichkr1,end=900) iel,iord,rrnu,(rx(l1),l1=1,nn)
      rnu = rrnu
         do 801 l1=1,nn
         rk(1,l1)=rx(l1)
  801    continue
      if (noyes.eq.'y') then
         read(ichkr2) iel,iord,rrnu,(rx(l1),l1=1,nn)
      else
         read(ichkr2) (rx(l1),l1=1,nn)
      endif
         do 802 l1=1,nn
         rk(2,l1)=rx(l1)
  802    continue

      sum1 = 0.
      sum2 = 0.
      nn1 = nn-1
	 do 200 l1=1,nn
c        write(*,*) x(l1),rk(1,l1),delta(1,l1),delta(2,l1)
	 dummy(l1) = rk(1,l1)*delta(1,l1)
	 dummy(l1+2000) = rk(2,l1)*delta(2,l1)
  200    continue
      call vinta(x,dummy,dummy(4001),nn,1,1)
      call vinta(x,dummy(2001),dummy(6001),nn,1,1)

      do 350 l1=1,nn
      write(20,1001) x(l1),delta(1,l1),rk(1,l1),dummy(4000+l1)
      write(21,1001) x(l1),delta(2,l1),rk(2,l1),dummy(6000+l1)
      write(*,1001) x(l1),delta(1,l1),delta(2,l1),
     +    dummy(4000+l1),dummy(6000+l1)
  350 continue

      
      sum1 = dummy(4000+nn)
      sum2 = dummy(6000+nn)
      ijk = 1
      if (ijk.eq.1) goto 250
         do 110 l1=1,nn1
         d = (x(l1+1)-x(l1))/2.
         sum1 = sum1+rk(1,l1)*delta(1,l1)*d
     +             +rk(1,l1+1)*delta(1,l1+1)*d
         sum2 = sum2+rk(2,l1)*delta(2,l1)*d
     +             +rk(2,l1+1)*delta(2,l1+1)*d
  110    continue
  250 write(6,*) iel,iord,sum1,sum2,sum1+sum2
      write(14,1004) iel,iord,rnu,sum1,sum2,sum1+sum2
      goto 105
  900 stop
 1000 format(i5)
 1001 format(1p5e13.5)
 1003 format(i5,f10.3,e15.7)
 1004 format(2i6,f11.4,2x,1p2e15.7,1x,1pe15.7)
      end
      
      subroutine stop1
      ijk = 1
      if (ijk.eq.1) stop
      return
      end

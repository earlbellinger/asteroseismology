      subroutine kercg(x,nn,aa,aax,rkcro,rkgro,ik)
c
c     calculates kernels for ln(c^2) and ln(gamma1) 
c       given equilibrium state and kernels for ln(c^2) and ln(rho)
c
      implicit double precision(a-h,o-z)
      integer*4 v
      logical zero1
      parameter (nnmax=1500)
      dimension x(nn), aa(5,nn), aax(3,nn)
      dimension rkcro(ik,nn),rkcg(ik,nn)
      dimension ea(2,3),v(2)
      common/consts/twopi,fourpi,gdimen,gconst,tiny
      common/centre/zero1,inonz
      common/solde/a11(nnmax),a12(nnmax),a21(nnmax),a22(nnmax),
     +   b1(nnmax),b2(nnmax),g11,g12,g13,g21,g22,g23
      external rhs,bc

      nn1 = nn-1

      call zero(a11,nn)
      call zero(b1,nn)
      call zero(b2,nn)
         do 100 l1=2,nn
	 a12(l1) = 1.
	 a21(l1) = -4.*aa(2,l1)*aa(3,l1)/x(l1)**2
	 a22(l1) = -a21(l1)*x(l1)/2. + 2./x(l1) + aax(3,l1)/aax(1,l1)
  100    continue
c     put (rho r^2/p) d(Kroc)/dr into b2
      call derivk(x,rkcro(2,1),b2,1,nn,2,1,2)
         do 110 l1=2,nn
         b2(l1) = b2(l1)*aa(3,l1)*aa(2,l1)/aa(1,l1)
  110    continue

      ii = 2
      kk = 0
      ka = 1
      kb = 1
      ki = 0
      id = 2
      v(1) = 1
      v(2) = 2
      ucy = 1.

  150 continue 
      
      call nrk(x(2),rkcg(1,2),zk,ap,aq,rhs,bc,ii,kk,ka,kb,ki,nn1,id,ucy,ea,
     +         det,v)
      write(*,*) 'ea:'
      write(*,*) ea
      write(*,*) 'Continue? (0/1)'
      read(*,*) iyesno
      if (iyesno.eq.999) goto 200
      if (iyesno.eq.1) then
         goto 100
      else
         stop
      endif

  200 return
      end

     
      subroutine rhs(x,y,zk,ap,aq,f,fd,h,hd,ifd,n)
      implicit real*8 (a-h,o-z)
      parameter (nnmax=1500)
      dimension y(1),f(1),fd(ifd,1),h(1),hd(ifd,1),zk(1),ap(1),aq(1)
      common/solde/a11(nnmax),a12(nnmax),a21(nnmax),a22(nnmax),
     +   b1(nnmax),b2(nnmax),g11,g12,g13,g21,g22,g23
      fd(1,1) = 0.
      fd(1,2) = 1.
      f(1) = fd(1,1)*y(1) + fd(1,2)*y(2)
      fd(2,1) = a21(n+1)
      fd(2,2) = a22(n+1)
      f(2) = fd(2,1)*y(1) + fd(2,2)*y(2) + b2(n+1)
      return
      end
 
      subroutine bc(x1,x2,y1,y2,zk,ap,aq,g,gd,ig,id,n)
      implicit real*8(a-h,o-z)
      parameter (nnmax=1500)
      dimension y1(1),y2(1),g(1),gd(ig,1),zk(1),ap(1),aq(1)
      common/solde/a11(nnmax),a12(nnmax),a21(nnmax),a22(nnmax),
     +   b1(nnmax),b2(nnmax),g11,g12,g13,g21,g22,g23
      gd(1,1) = -3.
      gd(1,2) = x1
      g(1) = gd(1,1)*y1(1) + gd(1,2)*y1(2)
      gd(2,1) = 1.
      gd(2,2) = 0.
      g(2) = gd(2,1)*y2(1) + gd(2,2)*y2(2) - 1.
      return
      end


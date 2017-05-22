      implicit real*8(a-h,o-z)
      integer*4 v
      parameter(nnmax=100000)
      dimension x(nnmax),y(2,nnmax),ea(2,3),v(2)
      common/work/dummy(1000000)
      external rhs, bc

      write(*,*) 'Input nn'
      read(*,*) nn
         do 90 l1=1,nn
         x(l1) = float(l1-1)/float(nn-1)*1.0
         y(1,l1) = 0.
         y(2,l1) = 1.
   90    continue

      s10 = sin(10.)
      c10 = cos(10.)
      s5  = sin(5.)
      c5  = cos(5.)
      t10 = s10/c10
      bb  = 1./(1.+t10)
      aa  = 1. - bb

      ii = 2
      kk = 0
      ka = 1
      kb = 1
      ki = 0
      id = 2
      v(1) = 1
      v(2) = 2
      ucy = 1.

  100 continue
      call nrk(x,y,zk,ap,aq,rhs,bc,ii,kk,ka,kb,ki,nn,id,ucy,ea,
     +         det,v)
      write(*,*) 'det = ',det
      write(*,*) 'ea:'
      write(*,*) ea
      write(*,*) ' '
      write(*,*) 'v: ',v
      nnm = (nn+1)/2.
      write(*,*) 'At x=0:  ',y(1,1),y(2,1)
      write(*,*) 'x(nn) =  ',x(nn)
      write(*,*) 'Midpoint:',x(nnm),y(1,nnm),y(1,nnm)-(aa*s5+bb*c5)
      write(*,*) 'Continue? (0/1)'
      read(*,*) i
      if (i.eq.1) goto 100
      stop
      end
      
     
      subroutine rhs(x,y,zk,ap,aq,f,fd,h,hd,ifd,n)
      implicit real*8 (a-h,o-z)
      dimension y(1),f(1),fd(ifd,1),h(1),hd(ifd,1),zk(1),ap(1),aq(1)
      fd(1,1) = 0.
      fd(1,2) = 1.
      f(1) = fd(1,1)*y(1) + fd(1,2)*y(2)
      fd(2,1) = 0.
      fd(2,2) = 0.
      f(2) = fd(2,1)*y(1) + fd(2,2)*y(2) + x
c     write(*,*) 'In RHS, n,y = ',n,y(1),y(2)
c     write(*,*) '          f = ',f(1),f(2)
      return
      end

      subroutine bc(x1,x2,y1,y2,zk,ap,aq,g,gd,ig,id,n)
      implicit real*8(a-h,o-z)
      dimension y1(1),y2(1),g(1),gd(ig,1),zk(1),ap(1),aq(1)
      gd(1,1) = 1.
      gd(1,2) = 0.
      g(1) = gd(1,1)*y1(1) + gd(1,2)*y1(2)
      gd(2,1) = 1.
      gd(2,2) = 0.
      g(2) = gd(2,1)*y2(1) + gd(2,2)*y2(2) - 1.
c     write(*,*) 'In BC, y1,y2,g = ',y1(1),y1(2)
c     write(*,*) '                 ',y2(1),y2(2)
      write(*,*) '                 ',g(1),g(2)
      return
      end
